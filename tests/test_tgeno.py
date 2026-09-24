import unittest

from dna_compare.config import default_settings
from dna_compare.eigenstrat import HEADER_SIZE, PackedTGeno


class TGenoTests(unittest.TestCase):
    def test_tgeno_header_math_matches_aadr_file(self):
        path = default_settings().aadr_geno
        if not path.exists():
            self.skipTest("AADR geno not present")
        reader = PackedTGeno(path)
        self.assertEqual(reader.nind, 27594)
        self.assertEqual(reader.nsnp, 584131)
        self.assertEqual(path.stat().st_size, HEADER_SIZE + reader.bytes_per_ind * reader.nind)
        geno = reader.read_individual(0)
        self.assertEqual(geno.shape, (reader.nsnp,))
        self.assertTrue(set(int(x) for x in set(geno[:200].tolist())).issubset({0, 1, 2, 3}))


if __name__ == "__main__":
    unittest.main()
