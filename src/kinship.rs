use std::collections::HashMap;

use pyo3::prelude::*;
use pyo3::types::PyDict;

const MIN_OVERLAP: usize = 200;
const THIN_BP: i64 = 100_000;
const THIN_IF_AT_LEAST: usize = 2500;
const MIN_IGC: f64 = 0.15;
const MIN_GQ: i32 = 7;

fn is_autosome(chrom: &str) -> bool {
    chrom.parse::<u32>().map(|c| (1..=22).contains(&c)).unwrap_or(false)
}

fn normalize_rsid(raw: &str) -> Option<String> {
    let token = raw.split(';').next()?.split(',').next()?.trim();
    if token.is_empty() || token == "." {
        return None;
    }
    let lower = token.to_ascii_lowercase();
    if lower.starts_with("rs") && lower.len() > 2 {
        Some(lower)
    } else {
        None
    }
}

fn diploid_alt_dosage(genotype: &str) -> Option<f32> {
    let token = genotype.split(':').next().unwrap_or("").replace('|', "/");
    if token.contains('/') {
        let parts: Vec<&str> = token.split('/').collect();
        if parts.len() != 2 || parts.iter().any(|p| *p == ".") {
            return None;
        }
        let a: i32 = parts[0].parse().ok()?;
        let b: i32 = parts[1].parse().ok()?;
        return Some((a + b) as f32);
    }
    match token.as_str() {
        "0" => Some(0.0),
        "1" => Some(2.0),
        _ => None,
    }
}

fn qc_ok(filter: &str, igc: f64, gq: i32) -> bool {
    let filt = filter.trim();
    if !filt.is_empty() && filt != "." && !filt.eq_ignore_ascii_case("pass") && !filt.to_ascii_uppercase().starts_with("PASS") {
        return false;
    }
    if igc >= 0.0 && igc < MIN_IGC {
        return false;
    }
    if gq >= 0 && gq < MIN_GQ {
        return false;
    }
    true
}

fn aligned_dosages(
    l_ref: u8,
    l_alt: u8,
    l_dosage: f32,
    r_ref: u8,
    r_alt: u8,
    r_dosage: f32,
) -> Option<(f32, f32)> {
    let complement = |b: u8| -> u8 {
        match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            x => x,
        }
    };
    if l_ref == r_ref && l_alt == r_alt {
        return Some((l_dosage, r_dosage));
    }
    if l_ref == r_alt && l_alt == r_ref {
        return Some((l_dosage, 2.0 - r_dosage));
    }
    let lr = complement(l_ref);
    let la = complement(l_alt);
    if lr == r_ref && la == r_alt {
        return Some((l_dosage, r_dosage));
    }
    if lr == r_alt && la == r_ref {
        return Some((l_dosage, 2.0 - r_dosage));
    }
    None
}

#[pyfunction]
pub fn king_kinship<'py>(
    py: Python<'py>,
    q_chroms: Vec<String>,
    q_positions: Vec<i64>,
    q_ref: Vec<String>,
    q_alt: Vec<String>,
    q_dosage: Vec<f64>,
    q_genotype: Vec<String>,
    q_rsid: Vec<String>,
    q_filter: Vec<String>,
    q_igc: Vec<f64>,
    q_gq: Vec<i32>,
    o_chroms: Vec<String>,
    o_positions: Vec<i64>,
    o_ref: Vec<String>,
    o_alt: Vec<String>,
    o_dosage: Vec<f64>,
    o_genotype: Vec<String>,
    o_rsid: Vec<String>,
    o_filter: Vec<String>,
    o_igc: Vec<f64>,
    o_gq: Vec<i32>,
) -> PyResult<Bound<'py, PyDict>> {
    let qn = q_chroms.len();
    let on = o_chroms.len();
    if q_positions.len() != qn
        || q_ref.len() != qn
        || q_alt.len() != qn
        || q_dosage.len() != qn
        || q_genotype.len() != qn
        || q_rsid.len() != qn
        || q_filter.len() != qn
        || q_igc.len() != qn
        || q_gq.len() != qn
    {
        return Err(pyo3::exceptions::PyValueError::new_err("query arrays length mismatch"));
    }
    if o_positions.len() != on
        || o_ref.len() != on
        || o_alt.len() != on
        || o_dosage.len() != on
        || o_genotype.len() != on
        || o_rsid.len() != on
        || o_filter.len() != on
        || o_igc.len() != on
        || o_gq.len() != on
    {
        return Err(pyo3::exceptions::PyValueError::new_err("other arrays length mismatch"));
    }

    let mut other_by_pos: HashMap<(String, i64), usize> = HashMap::new();
    let mut other_by_rsid: HashMap<String, usize> = HashMap::new();
    for i in 0..on {
        if is_autosome(&o_chroms[i]) {
            other_by_pos.insert((o_chroms[i].clone(), o_positions[i]), i);
            if let Some(rsid) = normalize_rsid(&o_rsid[i]) {
                other_by_rsid.entry(rsid).or_insert(i);
            }
        }
    }

    struct Pair {
        chrom: String,
        pos: i64,
        da: f32,
        db: f32,
        how: &'static str,
    }
    let mut pairs: Vec<Pair> = Vec::new();
    let mut n_qc_dropped = 0usize;
    let mut used_other: HashMap<(String, i64), bool> = HashMap::new();

    for i in 0..qn {
        if !is_autosome(&q_chroms[i]) {
            continue;
        }
        let key = (q_chroms[i].clone(), q_positions[i]);
        let (how, o_idx) = if let Some(&idx) = other_by_pos.get(&key) {
            ("pos", idx)
        } else if let Some(rsid) = normalize_rsid(&q_rsid[i]) {
            match other_by_rsid.get(&rsid) {
                Some(&idx) => ("rsid", idx),
                None => continue,
            }
        } else {
            continue;
        };
        let o_chrom = o_chroms[o_idx].clone();
        let o_pos = o_positions[o_idx];
        let rkey = (o_chrom, o_pos);
        if used_other.contains_key(&rkey) {
            continue;
        }
        if !qc_ok(&q_filter[i], q_igc[i], q_gq[i]) || !qc_ok(&o_filter[o_idx], o_igc[o_idx], o_gq[o_idx]) {
            n_qc_dropped += 1;
            continue;
        }
        let l_ref = q_ref[i].as_bytes().first().copied().unwrap_or(b'N');
        let l_alt = q_alt[i].as_bytes().first().copied().unwrap_or(b'N');
        let mut l_dosage = q_dosage[i] as f32;
        if l_dosage < 0.0 {
            l_dosage = diploid_alt_dosage(&q_genotype[i]).unwrap_or(-1.0);
        }
        let r_ref = o_ref[o_idx].as_bytes().first().copied().unwrap_or(b'N');
        let r_alt = o_alt[o_idx].as_bytes().first().copied().unwrap_or(b'N');
        let mut r_dosage = o_dosage[o_idx] as f32;
        if r_dosage < 0.0 {
            r_dosage = diploid_alt_dosage(&o_genotype[o_idx]).unwrap_or(-1.0);
        }
        if let Some((da, db)) = aligned_dosages(l_ref, l_alt, l_dosage, r_ref, r_alt, r_dosage) {
            used_other.insert(rkey, true);
            pairs.push(Pair {
                chrom: key.0,
                pos: key.1,
                da,
                db,
                how,
            });
        }
    }

    let n_pruned = if pairs.len() >= THIN_IF_AT_LEAST {
        let original_len = pairs.len();
        pairs.sort_by(|a, b| {
            let ac = a.chrom.parse::<u32>().unwrap_or(99);
            let bc = b.chrom.parse::<u32>().unwrap_or(99);
            (ac, a.pos).cmp(&(bc, b.pos))
        });
        let mut kept: Vec<Pair> = Vec::new();
        let mut last_chrom = String::new();
        let mut last_pos = i64::MIN / 2;
        for p in &pairs {
            if p.chrom != last_chrom || p.pos - last_pos >= THIN_BP {
                kept.push(Pair {
                    chrom: p.chrom.clone(),
                    pos: p.pos,
                    da: p.da,
                    db: p.db,
                    how: p.how,
                });
                last_chrom = p.chrom.clone();
                last_pos = p.pos;
            }
        }
        let pruned = if kept.len() >= MIN_OVERLAP {
            original_len - kept.len()
        } else {
            0
        };
        if pruned > 0 {
            pairs = kept;
        }
        pruned
    } else {
        0
    };

    let n_matched_pos = pairs.iter().filter(|p| p.how == "pos").count();
    let n_matched_rsid = pairs.iter().filter(|p| p.how == "rsid").count();
    let mut ibs0 = 0usize;
    let mut ibs1 = 0usize;
    let mut ibs2 = 0usize;
    let mut shared_sum = 0.0f64;
    let mut het_a = 0usize;
    let mut het_b = 0usize;
    let mut both_het = 0usize;
    for p in &pairs {
        let da = p.da;
        let db = p.db;
        let shared = 2.0 - (da - db).abs();
        if (da == 0.0 || da == 2.0) && (db == 0.0 || db == 2.0) && da != db {
            ibs0 += 1;
        } else if da == 1.0 && db == 1.0 {
            ibs2 += 1;
            both_het += 1;
        } else if da != db {
            ibs1 += 1;
        } else {
            ibs2 += 1;
        }
        if da == 1.0 {
            het_a += 1;
        }
        if db == 1.0 {
            het_b += 1;
        }
        shared_sum += shared as f64;
    }
    let n = ibs0 + ibs1 + ibs2;
    let mean_ibs = if n > 0 { shared_sum / (2.0 * n as f64) } else { 0.0 };
    let het_sum = het_a + het_b;
    let kinship = if het_sum == 0 {
        if ibs0 == 0 { 0.5 } else { 0.0 }
    } else {
        (both_het as f64 - 2.0 * ibs0 as f64) / het_sum as f64
    };

    let out = PyDict::new(py);
    out.set_item("n_snps", n)?;
    out.set_item("n_matched_pos", n_matched_pos)?;
    out.set_item("n_matched_rsid", n_matched_rsid)?;
    out.set_item("n_qc_dropped", n_qc_dropped)?;
    out.set_item("n_pruned", n_pruned)?;
    out.set_item("het_rate_query", if n > 0 { het_a as f64 / n as f64 } else { 0.0 })?;
    out.set_item("het_rate_other", if n > 0 { het_b as f64 / n as f64 } else { 0.0 })?;
    out.set_item("mean_ibs", mean_ibs)?;
    out.set_item("kinship", kinship)?;
    out.set_item("ibs0", ibs0)?;
    out.set_item("ibs1", ibs1)?;
    out.set_item("ibs2", ibs2)?;
    Ok(out)
}
