from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass
class VariantRow:
    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str
    genotype: str
    dosage_alt: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class VcfSummary:
    sample_id: str
    n_records: int
    n_snps: int
    n_non_snp_skipped: int
    n_samples_in_file: int
    chrom_counts: dict[str, int] = field(default_factory=dict)
    preview: list[VariantRow] = field(default_factory=list)
    reference: str | None = None
    source: str | None = None
    contig_lengths: dict[str, int] = field(default_factory=dict)
    assembly: str = "unknown"
    assembly_evidence: str = ""
    lifted_to: str | None = None
    n_lifted: int = 0
    n_unmapped: int = 0
    snp_catalog_id: str | None = None
    n_catalog: int = 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "n_records": self.n_records,
            "n_snps": self.n_snps,
            "n_non_snp_skipped": self.n_non_snp_skipped,
            "n_samples_in_file": self.n_samples_in_file,
            "chrom_counts": self.chrom_counts,
            "preview": [row.to_dict() for row in self.preview],
            "reference": self.reference,
            "source": self.source,
            "contig_lengths": self.contig_lengths,
            "assembly": self.assembly,
            "assembly_evidence": self.assembly_evidence,
            "lifted_to": self.lifted_to,
            "n_lifted": self.n_lifted,
            "n_unmapped": self.n_unmapped,
            "snp_catalog_id": self.snp_catalog_id,
            "n_catalog": self.n_catalog,
        }


@dataclass
class ReferenceStatus:
    key: str
    filename: str
    present: bool
    path: str | None = None
    note: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class HomininEstimate:
    label: str
    method: str
    percent: float | None
    n_snps: int
    detail: dict[str, Any] = field(default_factory=dict)
    available: bool = True
    message: str | None = None

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        return payload


@dataclass
class PopulationEstimate:
    population: str
    percent: float
    n_snps: int
    mean_ibs: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


CasteEstimate = PopulationEstimate


@dataclass
class AdditionalDetail:
    key: str
    topic: str
    finding: str
    evidence: str = ""
    source: str = ""
    available: bool = True
    confidence: str = ""
    n_snps: int = 0
    n_markers: int = 0
    coverage_status: str = ""
    parent: str | None = None
    message: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class DiseaseEstimate:
    trait: str
    pgs_id: str
    score: float | None
    n_snps: int
    n_score: int
    coverage_pct: float | None = None
    coverage_status: str = ""
    trait_category: str = "disease"
    percentile: float | None = None
    relative_level: str = ""
    citation: str = ""
    available: bool = True
    message: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class CommunityRefMatch:
    population: str
    percent: float
    n_snps: int = 0
    panel: str = ""
    ref_aasi_label: str = "Ref AASI"
    ref_steppe_label: str = "Ref Steppe"
    ref_indus_label: str = "Ref Indus"
    ref_east_asian_label: str = ""
    sample_aasi_label: str = "This AASI"
    sample_steppe_label: str = "This Steppe"
    sample_indus_label: str = "This Indus"
    sample_east_asian_label: str = ""
    aasi_in_range: bool = False
    steppe_in_range: bool = False
    indus_in_range: bool = False
    east_asian_in_range: bool = False
    aasi_score: float | None = None
    steppe_score: float | None = None
    indus_score: float | None = None
    east_asian_score: float | None = None
    y_score: float | None = None
    sample_aasi: float | None = None
    sample_steppe: float | None = None
    sample_indus: float | None = None
    sample_east_asian: float | None = None
    ref_aasi: str = ""
    ref_steppe: str = ""
    ref_indus: str = ""
    ref_east_asian: str = ""
    ref_y: str = ""
    sample_y: str | None = None
    y_note: str = ""
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class HaplogroupMarkerCall:
    haplogroup: str
    marker: str
    rsid: str
    chrom: str
    pos: int
    ancestral: str
    derived: str
    observed: str | None
    status: str
    genotype: str | None = None
    backbone: bool = True
    qual: float | None = None
    gq: int | None = None
    dp: int | None = None
    igc: float | None = None
    ad: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class HaplogroupGroupCount:
    n: int
    n_called: int
    percent: float | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class HaplogroupRow:
    haplogroup: str
    marker: str
    sample_status: str
    groups: dict[str, HaplogroupGroupCount] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "haplogroup": self.haplogroup,
            "marker": self.marker,
            "sample_status": self.sample_status,
            "groups": {name: item.to_dict() for name, item in self.groups.items()},
        }


@dataclass
class HaplogroupResult:
    available: bool
    sample_best: str | None = None
    markers: list[HaplogroupMarkerCall] = field(default_factory=list)
    rows: list[HaplogroupRow] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    mt_available: bool = False
    mt_sample_best: str | None = None
    mt_markers: list[HaplogroupMarkerCall] = field(default_factory=list)
    mt_rows: list[HaplogroupRow] = field(default_factory=list)
    mt_notes: list[str] = field(default_factory=list)
    status_notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "available": self.available,
            "sample_best": self.sample_best,
            "markers": [item.to_dict() for item in self.markers],
            "rows": [item.to_dict() for item in self.rows],
            "notes": self.notes,
            "mt_available": self.mt_available,
            "mt_sample_best": self.mt_sample_best,
            "mt_markers": [item.to_dict() for item in self.mt_markers],
            "mt_rows": [item.to_dict() for item in self.mt_rows],
            "mt_notes": self.mt_notes,
            "status_notes": self.status_notes,
        }


@dataclass
class RelatednessResult:
    available: bool
    other_filename: str | None = None
    other_sample_id: str | None = None
    query_sample_id: str | None = None
    n_snps: int = 0
    n_matched_pos: int = 0
    n_matched_rsid: int = 0
    n_qc_dropped: int = 0
    n_pruned: int = 0
    het_rate_query: float | None = None
    het_rate_other: float | None = None
    reliability: str | None = None
    mean_ibs: float | None = None
    kinship: float | None = None
    ibs0: int = 0
    ibs1: int = 0
    ibs2: int = 0
    ibs0_pct: float | None = None
    ibs1_pct: float | None = None
    ibs2_pct: float | None = None
    shared_pct: float | None = None
    relationship: str | None = None
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ComparisonBlock:
    kind: str
    available: bool
    estimates: list[Any] = field(default_factory=list)
    missing_files: list[ReferenceStatus] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    hidden: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "available": self.available,
            "estimates": [item.to_dict() for item in self.estimates],
            "missing_files": [item.to_dict() for item in self.missing_files],
            "notes": self.notes,
            "hidden": self.hidden,
        }


@dataclass
class AnalysisResult:
    """JSON-serializable payload for a future POST /analyze API and HTML UI."""

    ok: bool
    source_filename: str
    vcf: VcfSummary | None
    hominin: ComparisonBlock
    caste: ComparisonBlock
    populations: ComparisonBlock
    ancestry: ComparisonBlock
    haplogroups: HaplogroupResult
    ancestry_5: ComparisonBlock | None = None
    community_ref: ComparisonBlock | None = None
    disease: ComparisonBlock | None = None
    additional: ComparisonBlock | None = None
    relatedness: RelatednessResult | None = None
    errors: list[str] = field(default_factory=list)
    timings: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        relatedness = self.relatedness or RelatednessResult(
            available=False,
            notes=["No second VCF uploaded."],
        )
        return {
            "ok": self.ok,
            "source_filename": self.source_filename,
            "vcf": None if self.vcf is None else self.vcf.to_dict(),
            "hominin": self.hominin.to_dict(),
            "caste": self.caste.to_dict(),
            "populations": self.populations.to_dict(),
            "ancestry": self.ancestry.to_dict(),
            "ancestry_5": (
                self.ancestry_5 or ComparisonBlock(kind="ancestry_5", available=False)
            ).to_dict(),
            "haplogroups": self.haplogroups.to_dict(),
            "community_ref": (
                self.community_ref or ComparisonBlock(kind="community_ref", available=False)
            ).to_dict(),
            "disease": (
                self.disease or ComparisonBlock(kind="disease", available=False)
            ).to_dict(),
            "additional": (
                self.additional or ComparisonBlock(kind="additional", available=False)
            ).to_dict(),
            "relatedness": relatedness.to_dict(),
            "errors": self.errors,
            "timings": self.timings,
        }
