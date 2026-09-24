use std::fs::File;
use std::path::Path;

use memmap2::Mmap;
use numpy::PyArray1;
use pyo3::prelude::*;
use rayon::prelude::*;

const HEADER_SIZE: usize = 48;
const PACKED_MISSING: u8 = 3;

pub struct TGenoMmap {
    mmap: Mmap,
    nsnp: usize,
    bytes_per_ind: usize,
}

impl TGenoMmap {
    pub fn open(path: &str) -> std::io::Result<Self> {
        let file = File::open(Path::new(path))?;
        let mmap = unsafe { Mmap::map(&file)? };
        if mmap.len() < HEADER_SIZE {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "TGENO header too short",
            ));
        }
        let header = &mmap[..HEADER_SIZE];
        if !header.starts_with(b"TGENO") {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "not a TGENO file",
            ));
        }
        let parts: Vec<&[u8]> = header.split(|b| *b == b' ').filter(|p| !p.is_empty()).collect();
        if parts.len() < 3 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "invalid TGENO header",
            ));
        }
        let nind = parse_usize(parts[1])?;
        let nsnp = parse_usize(parts[2])?;
        let bytes_per_ind = (nsnp + 3) / 4;
        let expected = HEADER_SIZE + bytes_per_ind * nind;
        if mmap.len() != expected {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "TGENO size mismatch",
            ));
        }
        Ok(Self {
            mmap,
            nsnp,
            bytes_per_ind,
        })
    }

    fn decode_individual(&self, ind_index: usize, out: &mut [f32]) {
        let offset = HEADER_SIZE + ind_index * self.bytes_per_ind;
        let packed = &self.mmap[offset..offset + self.bytes_per_ind];
        let mut idx = 0usize;
        for &byte in packed {
            for shift in [6u8, 4, 2, 0] {
                if idx >= self.nsnp {
                    break;
                }
                let g = (byte >> shift) & 3;
                out[idx] = if g == PACKED_MISSING {
                    f32::NAN
                } else {
                    g as f32
                };
                idx += 1;
            }
        }
    }
}

fn parse_usize(bytes: &[u8]) -> std::io::Result<usize> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidData, "bad header int"))?;
    text.parse::<usize>()
        .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidData, "bad header int"))
}

fn open_tgeno(geno_path: &str, nsnp: usize) -> Result<TGenoMmap, String> {
    let tgeno = TGenoMmap::open(geno_path).map_err(|e| format!("TGENO open failed: {e}"))?;
    if tgeno.nsnp != nsnp {
        return Err(format!(
            "nsnp mismatch: expected {nsnp}, got {}",
            tgeno.nsnp
        ));
    }
    Ok(tgeno)
}

fn freq_from_partials(
    partials: Vec<(Vec<f64>, Vec<u32>)>,
    nsnp: usize,
) -> Vec<f32> {
    let mut sums = vec![0.0f64; nsnp];
    let mut counts = vec![0u32; nsnp];
    for (local_sums, local_counts) in partials {
        for i in 0..nsnp {
            sums[i] += local_sums[i];
            counts[i] += local_counts[i];
        }
    }
    let mut freq = vec![f32::NAN; nsnp];
    for i in 0..nsnp {
        if counts[i] > 0 {
            freq[i] = (sums[i] / counts[i] as f64) as f32;
        }
    }
    freq
}

fn accumulate_individual(
    tgeno: &TGenoMmap,
    ind_index: usize,
    nsnp: usize,
) -> (Vec<f64>, Vec<u32>) {
    let mut dosages = vec![0.0f32; nsnp];
    tgeno.decode_individual(ind_index, &mut dosages);
    let mut local_sums = vec![0.0f64; nsnp];
    let mut local_counts = vec![0u32; nsnp];
    for (i, d) in dosages.iter().enumerate() {
        if d.is_finite() {
            local_sums[i] += (*d as f64) / 2.0;
            local_counts[i] += 1;
        }
    }
    (local_sums, local_counts)
}

#[pyfunction]
pub fn population_allele1_freq(
    py: Python<'_>,
    geno_path: &str,
    sample_indices: Vec<usize>,
    nsnp: usize,
) -> PyResult<Py<PyArray1<f32>>> {
    if sample_indices.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err("no sample indices"));
    }
    let freq = py
        .allow_threads(|| -> Result<Vec<f32>, String> {
            let tgeno = open_tgeno(geno_path, nsnp)?;
            let partials: Vec<(Vec<f64>, Vec<u32>)> = sample_indices
                .par_iter()
                .map(|&ind_index| accumulate_individual(&tgeno, ind_index, nsnp))
                .collect();
            Ok(freq_from_partials(partials, nsnp))
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
    Ok(PyArray1::from_vec(py, freq).into())
}

#[pyfunction]
pub fn build_population_freqs_single_pass(
    py: Python<'_>,
    geno_path: &str,
    nsnp: usize,
    ind_indices: Vec<usize>,
    ind_group_offsets: Vec<usize>,
    ind_group_ids: Vec<usize>,
    n_groups: usize,
) -> PyResult<Vec<Py<PyArray1<f32>>>> {
    if n_groups == 0 {
        return Ok(Vec::new());
    }
    if ind_indices.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err("no sample indices"));
    }
    if ind_group_offsets.len() != ind_indices.len() + 1 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "ind_group_offsets length mismatch",
        ));
    }
    let freqs = py
        .allow_threads(|| -> Result<Vec<Vec<f32>>, String> {
            let tgeno = open_tgeno(geno_path, nsnp)?;
            let partials: Vec<(Vec<Vec<f64>>, Vec<Vec<u32>>)> = ind_indices
                .par_iter()
                .enumerate()
                .map(|(ii, &ind_index)| {
                    let mut sums = vec![vec![0.0f64; nsnp]; n_groups];
                    let mut counts = vec![vec![0u32; nsnp]; n_groups];
                    let mut dosages = vec![0.0f32; nsnp];
                    tgeno.decode_individual(ind_index, &mut dosages);
                    let start = ind_group_offsets[ii];
                    let end = ind_group_offsets[ii + 1];
                    for &gid in &ind_group_ids[start..end] {
                        for i in 0..nsnp {
                            let d = dosages[i];
                            if d.is_finite() {
                                sums[gid][i] += (d as f64) / 2.0;
                                counts[gid][i] += 1;
                            }
                        }
                    }
                    (sums, counts)
                })
                .collect();

            let mut merged_sums = vec![vec![0.0f64; nsnp]; n_groups];
            let mut merged_counts = vec![vec![0u32; nsnp]; n_groups];
            for (sums, counts) in partials {
                for gid in 0..n_groups {
                    for i in 0..nsnp {
                        merged_sums[gid][i] += sums[gid][i];
                        merged_counts[gid][i] += counts[gid][i];
                    }
                }
            }

            let mut out = Vec::with_capacity(n_groups);
            for gid in 0..n_groups {
                let mut freq = vec![f32::NAN; nsnp];
                for i in 0..nsnp {
                    if merged_counts[gid][i] > 0 {
                        freq[i] = (merged_sums[gid][i] / merged_counts[gid][i] as f64) as f32;
                    }
                }
                out.push(freq);
            }
            Ok(out)
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;

    Ok(freqs
        .into_iter()
        .map(|freq| PyArray1::from_vec(py, freq).into())
        .collect())
}
