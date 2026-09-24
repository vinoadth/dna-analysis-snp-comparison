mod align;
mod kinship;
mod tgeno;

use pyo3::prelude::*;

#[pymodule]
fn dna_compare_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(tgeno::population_allele1_freq, m)?)?;
    m.add_function(wrap_pyfunction!(tgeno::build_population_freqs_single_pass, m)?)?;
    m.add_function(wrap_pyfunction!(align::align_query_to_panel, m)?)?;
    m.add_function(wrap_pyfunction!(align::score_pgs_batch, m)?)?;
    m.add_function(wrap_pyfunction!(kinship::king_kinship, m)?)?;
    Ok(())
}
