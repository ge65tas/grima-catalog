use std::{
    fs::File,
    io::{BufWriter, Read, Seek, SeekFrom::Start, Write},
    path::Path,
};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum CatalogError {
    #[error("File has incomplete header data")]
    IncompleteHeader,
}

const BLOCK_SIZE: usize = 4096;
pub struct BinCatReader {
    file: File,
    band_size: f64,
    offsets: Vec<Vec<usize>>,
}

impl BinCatReader {
    pub fn new(file_path: impl AsRef<Path>) -> Result<BinCatReader, CatalogError> {
        let mut buffer = [0; BLOCK_SIZE];
        let mut byte_buffer = [0; 8];
        let mut file = File::open(file_path).unwrap();
        let mut bytes_read = file.read(&mut buffer).unwrap();
        if bytes_read < 16 {
            return Err(CatalogError::IncompleteHeader);
        }
        byte_buffer.copy_from_slice(&buffer[0..8]);
        let band_size = f64::from_le_bytes(byte_buffer);
        byte_buffer.copy_from_slice(&buffer[8..16]);
        let header_size = usize::from_le_bytes(byte_buffer);
        let mut offsets: Vec<usize> = Vec::with_capacity(header_size / 8 - 1);
        offsets.push(header_size);
        for i in 2..(header_size / 8) {
            let j = i % (BLOCK_SIZE / 8);
            if j == 0 {
                bytes_read = file.read(&mut buffer).unwrap();
                if bytes_read != BLOCK_SIZE && bytes_read < header_size % BLOCK_SIZE {
                    return Err(CatalogError::IncompleteHeader);
                }
            }
            byte_buffer.copy_from_slice(&buffer[(j * 8)..((j + 1) * 8)]);
            let offset = usize::from_le_bytes(byte_buffer);
            offsets.push(offset);
        }
        let mut suboffsets = Vec::with_capacity(offsets.len());
        for band in offsets.windows(2) {
            let start = band[0] / 8;
            let end = band[1] / 8;
            let mut current = Vec::with_capacity((end - start) / 8);
            for i in start..end {
                let j = i % (BLOCK_SIZE / 8);
                if j == 0 {
                    bytes_read = file.read(&mut buffer).unwrap();
                    if bytes_read != BLOCK_SIZE && bytes_read < end % BLOCK_SIZE {
                        return Err(CatalogError::IncompleteHeader);
                    }
                }
                byte_buffer.copy_from_slice(&buffer[(j * 8)..((j + 1) * 8)]);
                let suboffset = usize::from_le_bytes(byte_buffer);
                current.push(suboffset);
            }
            suboffsets.push(current);
        }
        println!("offsets: {:?}", offsets);
        for i in 0..(suboffsets.len() - 1) {
            println!("{}/{}", i + 1, suboffsets.len());
            let band_end = *suboffsets[i + 1].first().unwrap();
            suboffsets[i].push(band_end);
        }
        suboffsets
            .last_mut()
            .unwrap()
            .push(file.metadata().unwrap().len() as usize);
        println!(
            "first offsets: {:?}, \nlast offsets: {:?}",
            suboffsets.first().unwrap(),
            suboffsets.iter().last().unwrap()
        );
        println!("number of offsets: {}", offsets.len());
        Ok(Self {
            file,
            offsets: suboffsets,
            band_size,
        })
    }

    pub fn read_region(
        &mut self,
        ra_lower: f64,
        ra_upper: f64,
        dec_lower: f64,
        dec_upper: f64,
    ) -> Vec<[f64; 3]> {
        let lower_index_dec = ((dec_lower + 90.0) / (self.band_size)) as usize;
        let upper_index_dec =
            (((dec_upper + 90.0) / (self.band_size)) as usize).min(self.offsets.len() - 1);
        let mut out = Vec::new();
        for i in lower_index_dec..(upper_index_dec + 1) {
            let current_band = &self.offsets[i];
            let current_band_size = 360.0 / (current_band.len() as f64);
            let lower_index_ra = (ra_lower / current_band_size) as usize;
            let upper_index_ra =
                ((ra_upper / current_band_size) as usize + 1).min(current_band.len() - 1);
            println!(
                "{}, {}, {}, {}",
                lower_index_ra, upper_index_ra, lower_index_dec, upper_index_dec
            );
            let lower_offset = current_band[lower_index_ra];
            let upper_offset = current_band[upper_index_ra];
            let mut buffer: Vec<u8> = vec![0; upper_offset - lower_offset];
            self.file.seek(Start(lower_offset as u64)).unwrap();
            self.file.read_exact(&mut buffer).unwrap();
            out.extend(buffer.chunks(24).map(|x| {
                [
                    Self::bytes_to_float(&x[0..8]),
                    Self::bytes_to_float(&x[8..16]),
                    Self::bytes_to_float(&x[16..24]),
                ]
            }));
        }
        out
    }

    fn bytes_to_float(bytes: &[u8]) -> f64 {
        let mut byte_buffer = [0; 8];
        byte_buffer.copy_from_slice(bytes);
        f64::from_le_bytes(byte_buffer)
    }
}

pub fn write_bin_cat(file: impl AsRef<Path>, mut sources: Vec<[f64; 3]>, band_size: f64) {
    const BYTES_PER_ENTRY: usize = 8 * 3;
    let mut file = BufWriter::new(File::create(file).unwrap());
    let mut decs: Vec<f64> = Vec::new();
    sources.sort_unstable_by(|x, y| y[1].total_cmp(&x[1]));
    let mut current_pos = -90.0 + band_size;
    while current_pos < 90.0 {
        decs.push(current_pos);
        current_pos += band_size;
    }

    let mut bands = Vec::new();
    let mut suboffsets: Vec<Vec<usize>> = Vec::new();
    let mut band_start = 0;
    for (index, dec) in decs.iter().enumerate() {
        println!("{}/{}", index, decs.len());
        let lower_index = match sources
            .iter()
            .enumerate()
            .rfind(|x| x.1[1] > *dec)
            .map(|x| x.0)
        {
            Some(x) => x + 1,
            None => {
                if sources.iter().all(|x| x[1] < *dec) && !sources.is_empty() {
                    println!("last band!");
                    println!("{} sources remaining", sources.len());
                    0
                } else {
                    println!("empty band!");
                    suboffsets.push(vec![band_start]);
                    bands.push(Vec::new());
                    continue;
                }
            }
        };
        let mut current: Vec<[f64; 3]> = sources
            .iter()
            .rev()
            .take(sources.len() - lower_index)
            .cloned()
            .collect();
        sources.truncate(lower_index);
        current.sort_unstable_by(|x, y| x[0].total_cmp(&y[0]));
        let band_size_ra = band_size
            / dec
                .to_radians()
                .cos()
                .max((dec + band_size).to_radians().cos());
        let chunks = ((360.0 / band_size_ra) as usize).max(1);
        let current_offsets: Vec<usize> = (0..chunks)
            .map(|i| {
                let current_ra = i as f64 * band_size_ra;
                current.iter().filter(|&x| x[0] < current_ra).count() * BYTES_PER_ENTRY + band_start
            })
            .collect();
        band_start += current.len() * BYTES_PER_ENTRY;
        bands.push(current);
        suboffsets.push(current_offsets);
    }
    // start at size of header to give file-global offsets. band_size + header_size/file_body_start + band_end for each band + chunk_end for each chunk within the bands
    let subheader_size = (2 + bands.len()) * 8;
    let header_size = suboffsets.iter().map(|x| x.len()).sum::<usize>() * 8 + subheader_size;
    let offsets: Vec<usize> = suboffsets
        .iter()
        .scan(subheader_size, |state, x| {
            *state += x.len() * 8;
            Some(*state)
        })
        .collect();
    file.write_all(&band_size.to_le_bytes()).unwrap();
    file.write_all(&subheader_size.to_le_bytes()).unwrap();
    for offset in offsets {
        file.write_all(&offset.to_le_bytes()).unwrap();
    }
    for offset in suboffsets.into_iter().flatten() {
        let true_offset = offset + header_size;
        file.write_all(&true_offset.to_le_bytes()).unwrap();
    }
    for band in bands {
        for entry in band {
            for number in entry {
                file.write_all(&number.to_le_bytes()).unwrap();
            }
        }
    }
}
