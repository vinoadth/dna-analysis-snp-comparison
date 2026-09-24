use std::collections::HashMap;

use numpy::PyArray1;
use pyo3::prelude::*;
use pyo3::types::PyList;

fn query_allele1_freq(
    q_ref: u8,
    q_alt: u8,
    dosage: f32,
    allele1: u8,
    allele2: u8,
    haploid: bool,
) -> Option<f32> {
    if dosage < 0.0 {
        return None;
    }
    let copies = if haploid { 1.0 } else { 2.0 };
    let q_alt_freq = dosage / copies;
    if q_ref == allele2 && q_alt == allele1 {
        return Some(q_alt_freq);
    }
    if q_ref == allele1 && q_alt == allele2 {
        return Some(1.0 - q_alt_freq);
    }
    None
}

fn normalize_rsid(raw: &str) -> Option<String> {
    let token = raw.split(';').next()?.split(',').next()?.trim();
    if token.is_empty() || token == "." {
        return None;
    }
    let lower = token.to_ascii_lowercase();
    if lower.starts_with("rs") && lower.len() > 2 {
        return Some(lower);
    }
    None
}

#[pyfunction]
pub fn align_query_to_panel(
    py: Python<'_>,
    panel_chroms: Vec<String>,
    panel_positions: Vec<i64>,
    panel_allele1: Vec<String>,
    panel_allele2: Vec<String>,
    query_chroms: Vec<String>,
    query_positions: Vec<i64>,
    query_ref: Vec<String>,
    query_alt: Vec<String>,
    query_dosage: Vec<f64>,
) -> PyResult<(Py<PyArray1<f32>>, usize)> {
    let n = panel_chroms.len();
    if panel_positions.len() != n || panel_allele1.len() != n || panel_allele2.len() != n {
        return Err(pyo3::exceptions::PyValueError::new_err("panel arrays length mismatch"));
    }
    let qn = query_chroms.len();
    if query_positions.len() != qn
        || query_ref.len() != qn
        || query_alt.len() != qn
        || query_dosage.len() != qn
    {
        return Err(pyo3::exceptions::PyValueError::new_err("query arrays length mismatch"));
    }

    let mut query_map: HashMap<(String, i64), (u8, u8, f32, bool)> = HashMap::with_capacity(qn);
    for i in 0..qn {
        let chrom = query_chroms[i].clone();
        let pos = query_positions[i];
        let refb = query_ref[i].as_bytes().first().copied().unwrap_or(b'N');
        let altb = query_alt[i].as_bytes().first().copied().unwrap_or(b'N');
        if refb.is_ascii_alphabetic() && altb.is_ascii_alphabetic() {
            query_map.insert(
                (chrom, pos),
                (
                    refb.to_ascii_uppercase(),
                    altb.to_ascii_uppercase(),
                    query_dosage[i] as f32,
                    false,
                ),
            );
        }
    }

    let mut out = vec![f32::NAN; n];
    let mut n_overlap = 0usize;
    for i in 0..n {
        let a1 = panel_allele1[i].as_bytes().first().copied().unwrap_or(b'N').to_ascii_uppercase();
        let a2 = panel_allele2[i].as_bytes().first().copied().unwrap_or(b'N').to_ascii_uppercase();
        if let Some((q_ref, q_alt, dosage, haploid)) =
            query_map.get(&(panel_chroms[i].clone(), panel_positions[i]))
        {
            if let Some(freq) = query_allele1_freq(*q_ref, *q_alt, *dosage, a1, a2, *haploid) {
                out[i] = freq;
                n_overlap += 1;
            }
        }
    }
    Ok((PyArray1::from_vec(py, out).into(), n_overlap))
}

#[pyfunction]
pub fn score_pgs_batch<'py>(
    py: Python<'py>,
    rsids: Vec<String>,
    chroms: Vec<String>,
    positions: Vec<i64>,
    effects: Vec<String>,
    others: Vec<String>,
    weights: Vec<f64>,
    offsets: Vec<usize>,
    q_rsids: Vec<String>,
    q_chroms: Vec<String>,
    q_positions: Vec<i64>,
    q_ref: Vec<String>,
    q_alt: Vec<String>,
    q_dosage: Vec<f64>,
) -> PyResult<(Py<PyList>, Py<PyList>, Py<PyList>)> {
    let n = rsids.len();
    if chroms.len() != n
        || positions.len() != n
        || effects.len() != n
        || others.len() != n
        || weights.len() != n
    {
        return Err(pyo3::exceptions::PyValueError::new_err("variant arrays length mismatch"));
    }
    let qn = q_chroms.len();
    if q_positions.len() != qn || q_ref.len() != qn || q_alt.len() != qn || q_dosage.len() != qn {
        return Err(pyo3::exceptions::PyValueError::new_err("query arrays length mismatch"));
    }

    let mut by_pos: HashMap<(String, i64), (u8, u8, f32)> = HashMap::with_capacity(qn);
    let mut by_rsid: HashMap<String, (u8, u8, f32)> = HashMap::new();
    for i in 0..qn {
        let refb = q_ref[i].as_bytes().first().copied().unwrap_or(b'N').to_ascii_uppercase();
        let altb = q_alt[i].as_bytes().first().copied().unwrap_or(b'N').to_ascii_uppercase();
        let dosage = q_dosage[i] as f32;
        by_pos.insert((q_chroms[i].clone(), q_positions[i]), (refb, altb, dosage));
        if let Some(rsid) = normalize_rsid(&q_rsids[i]) {
            by_rsid.insert(rsid, (refb, altb, dosage));
        }
    }

    let complement = |b: u8| -> u8 {
        match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        }
    };

    let match_dosage = |row: (u8, u8, f32), effect: u8, other: u8| -> Option<f32> {
        let (refb, altb, dosage) = row;
        if !effect.is_ascii_alphabetic() || !altb.is_ascii_alphabetic() || !refb.is_ascii_alphabetic() {
            return None;
        }
        if altb == effect && (other == 0 || refb == other) {
            return Some(dosage);
        }
        if refb == effect && (other == 0 || altb == other) {
            return Some(2.0 - dosage);
        }
        if other != 0 && complement(altb) == effect && complement(refb) == other {
            return Some(dosage);
        }
        if other != 0 && complement(refb) == effect && complement(altb) == other {
            return Some(2.0 - dosage);
        }
        None
    };

    let batch_count = offsets.len().saturating_sub(1);
    let mut totals = vec![0.0f64; batch_count];
    let mut used = vec![0usize; batch_count];
    let mut counts = vec![0usize; batch_count];

    for batch in 0..batch_count {
        let start = offsets[batch];
        let end = offsets[batch + 1];
        counts[batch] = end - start;
        for vi in start..end {
            let effect = effects[vi].as_bytes().first().copied().unwrap_or(b'N').to_ascii_uppercase();
            let other = others[vi].as_bytes().first().copied().unwrap_or(b'0');
            let other = if others[vi].is_empty() { 0 } else { other.to_ascii_uppercase() };
            let row = if !rsids[vi].is_empty() {
                normalize_rsid(&rsids[vi]).and_then(|id| by_rsid.get(&id).copied())
            } else {
                None
            }
            .or_else(|| {
                if !chroms[vi].is_empty() && positions[vi] > 0 {
                    by_pos.get(&(chroms[vi].clone(), positions[vi])).copied()
                } else {
                    None
                }
            });
            if let Some(row) = row {
                if let Some(d) = match_dosage(row, effect, other) {
                    totals[batch] += d as f64 * weights[vi];
                    used[batch] += 1;
                }
            }
        }
    }

    let totals_py = PyList::new(py, totals)?.into();
    let used_py = PyList::new(py, used)?.into();
    let counts_py = PyList::new(py, counts)?.into();
    Ok((totals_py, used_py, counts_py))
}
