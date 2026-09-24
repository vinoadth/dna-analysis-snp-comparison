import tempfile
import unittest
from pathlib import Path

import numpy as np

from dna_compare._native import (
    align_query_to_panel_native,
    build_population_freqs_single_pass_native,
    native_available,
    population_allele1_freq_native,
)
from dna_compare.align import PanelAlignment, align_query_to_panel
from dna_compare.config import default_settings
from dna_compare.eigenstrat import AadrPanel
from dna_compare.vcf_parser import parse_vcf

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "samples" / "demo.snps.vcf"


class NativeExtensionTests(unittest.TestCase):
    def test_native_module_import(self):
        if native_available():
            import dna_compare_rs  # noqa: F401

    @unittest.skipUnless(
        Path(default_settings().aadr_geno).exists(),
        "AADR panel not present",
    )
    def test_population_freq_matches_python(self):
        settings = default_settings()
        panel = AadrPanel(settings)
        group = ("Mbuti",)
        py_freq = panel._population_allele1_freq_python(group)
        indices = [rec.index for rec in panel.samples_for_population("Mbuti")]
        rs_freq = population_allele1_freq_native(str(settings.aadr_geno), indices, panel.geno().nsnp)
        if rs_freq is None:
            self.skipTest("Rust extension not built")
        np.testing.assert_allclose(rs_freq, py_freq, rtol=1e-5, atol=1e-5, equal_nan=True)

    @unittest.skipUnless(
        Path(default_settings().aadr_geno).exists(),
        "AADR panel not present",
    )
    def test_alignment_matches_python(self):
        settings = default_settings()
        panel = AadrPanel(settings)
        _summary, index = parse_vcf(FIXTURE)
        py_align = align_query_to_panel(panel, index)
        assert py_align is not None
        native = align_query_to_panel_native(panel, index)
        if native is None:
            self.skipTest("Rust extension not built")
        rs_freqs, rs_overlap = native
        self.assertEqual(rs_overlap, py_align.n_overlap)
        np.testing.assert_allclose(rs_freqs, py_align.query_freqs, rtol=1e-5, atol=1e-5, equal_nan=True)

    @unittest.skipUnless(
        Path(default_settings().aadr_geno).exists(),
        "AADR panel not present",
    )
    def test_single_pass_freq_matches_python(self):
        settings = default_settings()
        with tempfile.TemporaryDirectory() as tmp:
            settings.cache_dir = Path(tmp)
            panel = AadrPanel(settings)
            groups = [("Mbuti",), ("Yoruba",)]
            expected = [panel._population_allele1_freq_python(group) for group in groups]
            ind_indices, ind_group_offsets, ind_group_ids, ind_to_groups = panel._single_pass_membership(
                groups
            )
            batch = build_population_freqs_single_pass_native(
                str(settings.aadr_geno),
                panel.geno().nsnp,
                ind_indices,
                ind_group_offsets,
                ind_group_ids,
                len(groups),
            )
            if batch is None:
                batch = panel._build_freq_caches_single_pass_python(groups, ind_indices, ind_to_groups)
            for got, want in zip(batch, expected, strict=True):
                np.testing.assert_allclose(got, want, rtol=1e-5, atol=1e-5, equal_nan=True)

            built, skipped = panel.build_all_population_freq_caches(groups)
            self.assertEqual(built, len(groups))
            self.assertEqual(skipped, 0)
            for group in groups:
                cached = panel.population_allele1_freq(group)
                np.testing.assert_allclose(
                    cached,
                    panel._population_allele1_freq_python(group),
                    rtol=1e-5,
                    atol=1e-5,
                    equal_nan=True,
                )

    def test_panel_alignment_dataclass(self):
        arr = np.array([1.0, np.nan, 0.5], dtype=np.float32)
        align = PanelAlignment(query_freqs=arr, n_overlap=2)
        self.assertEqual(align.n_overlap, 2)


if __name__ == "__main__":
    unittest.main()
