"""SNP VCF comparison library (HTML dashboard now, JSON API later)."""

from dna_compare.models import AnalysisResult
from dna_compare.service import AnalysisService

__all__ = ["AnalysisService", "AnalysisResult"]
