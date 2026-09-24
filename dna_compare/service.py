from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import BinaryIO, TextIO

from dna_compare.align import PanelAlignment, align_query_to_panel
from dna_compare.comparisons import (
    compare_ancestry,
    compare_ancestry_5source,
    compare_caste,
    compare_haplogroups,
    compare_hominin,
    compare_populations,
    compare_relatedness,
    compare_additional,
    compare_disease,
    score_community_reference,
)
from dna_compare.assembly import assembly_note, detect_assembly
from dna_compare.config import Settings, default_settings
from dna_compare.eigenstrat import AadrPanel
from dna_compare.genotype_index import build_rsid_index
from dna_compare.genotype_parser import parse_genotype_file
from dna_compare.liftover import ensure_hg38_to_hg19_chain, lift_query_index, load_chain
from dna_compare.models import AnalysisResult, ComparisonBlock, HaplogroupResult, RelatednessResult
from dna_compare.raw23_parser import Raw23ParseError
from dna_compare.snp_catalog import publish_snp_list
from dna_compare.population_groups import collect_population_groups
from dna_compare.timing import AnalyzeTimings, analyze_timings_enabled


class AnalysisService:
    """Stable entry point for CLI today and FastAPI/Flask later."""

    def __init__(self, settings: Settings | None = None):
        self.settings = settings or default_settings()
        self._panel: AadrPanel | None = None
        self._chain = None

    @property
    def panel(self) -> AadrPanel | None:
        panel = AadrPanel(self.settings)
        if not panel.available:
            return None
        if self._panel is None:
            self._panel = panel
        return self._panel

    def analyze(
        self,
        source: Path | str | TextIO | BinaryIO,
        *,
        filename: str | None = None,
        compare_hominin_flag: bool | None = None,
        compare_caste_flag: bool | None = None,
        compare_populations_flag: bool | None = None,
        compare_ancestry_flag: bool | None = None,
        compare_haplogroups_flag: bool | None = None,
        compare_disease_flag: bool | None = None,
        other_source: Path | str | TextIO | BinaryIO | None = None,
        other_filename: str | None = None,
    ) -> AnalysisResult:
        timings = AnalyzeTimings(enabled=analyze_timings_enabled(self.settings.analyze_timings))
        do_hominin = self.settings.compare_hominin if compare_hominin_flag is None else compare_hominin_flag
        do_caste = self.settings.compare_caste if compare_caste_flag is None else compare_caste_flag
        do_pops = (
            self.settings.compare_populations if compare_populations_flag is None else compare_populations_flag
        )
        do_ancestry = self.settings.compare_ancestry if compare_ancestry_flag is None else compare_ancestry_flag
        do_haplo = (
            self.settings.compare_haplogroups if compare_haplogroups_flag is None else compare_haplogroups_flag
        )
        do_disease = self.settings.compare_disease if compare_disease_flag is None else compare_disease_flag
        source_name = filename or (Path(source).name if isinstance(source, (str, Path)) else "upload.vcf")
        empty_pops = ComparisonBlock(kind="populations", available=False)
        empty_ancestry = ComparisonBlock(kind="ancestry", available=False)
        empty_ancestry_5 = ComparisonBlock(kind="ancestry_5", available=False)
        try:
            with timings.stage("parse"):
                summary, index = parse_genotype_file(
                    source,
                    filename=source_name,
                    preview_limit=self.settings.variant_preview_limit,
                )
        except Raw23ParseError as exc:
            return AnalysisResult(
                ok=False,
                source_filename=source_name,
                vcf=None,
                hominin=ComparisonBlock(kind="hominin", available=False, notes=[str(exc)]),
                caste=ComparisonBlock(kind="caste", available=False),
                populations=empty_pops,
                ancestry=empty_ancestry,
                ancestry_5=empty_ancestry_5,
                haplogroups=HaplogroupResult(available=False, notes=[str(exc)]),
                community_ref=ComparisonBlock(kind="community_ref", available=False, notes=[str(exc)]),
                disease=ComparisonBlock(kind="disease", available=False, notes=[str(exc)]),
                additional=ComparisonBlock(kind="additional", available=False, notes=[str(exc)]),
                relatedness=RelatednessResult(available=False, notes=[str(exc)]),
                errors=[str(exc)],
            )
        except Exception as exc:  # noqa: BLE001
            return AnalysisResult(
                ok=False,
                source_filename=source_name,
                vcf=None,
                hominin=ComparisonBlock(kind="hominin", available=False, notes=[str(exc)]),
                caste=ComparisonBlock(kind="caste", available=False),
                populations=empty_pops,
                ancestry=empty_ancestry,
                ancestry_5=empty_ancestry_5,
                haplogroups=HaplogroupResult(available=False, notes=[str(exc)]),
                community_ref=ComparisonBlock(kind="community_ref", available=False, notes=[str(exc)]),
                disease=ComparisonBlock(kind="disease", available=False, notes=[str(exc)]),
                additional=ComparisonBlock(kind="additional", available=False, notes=[str(exc)]),
                relatedness=RelatednessResult(available=False, notes=[str(exc)]),
                errors=[f"Failed to parse genotype file: {exc}"],
            )

        if filename:
            source_name = filename
        detect_assembly(summary, assume=self.settings.assume_assembly)
        with timings.stage("liftover"):
            index = self._maybe_lift(summary, index)
        with timings.stage("rsid_index"):
            by_rsid = build_rsid_index(index)

        panel = self.panel
        alignment: PanelAlignment | None = None
        needs_panel = do_hominin or do_caste or do_pops or do_ancestry
        if panel is not None and panel.available and needs_panel:
            with timings.stage("prefetch_freqs"):
                panel.prefetch_population_freqs(collect_population_groups(self.settings))
            with timings.stage("align_panel"):
                alignment = align_query_to_panel(panel, index)

        hominin = ComparisonBlock(kind="hominin", available=False, notes=["Hominin comparison disabled."])
        caste = ComparisonBlock(kind="caste", available=False, notes=["Caste comparison disabled."])
        populations = ComparisonBlock(kind="populations", available=False, notes=["Population comparison disabled."])
        ancestry = ComparisonBlock(kind="ancestry", available=False, notes=["Ancestry comparison disabled."])
        ancestry_5 = ComparisonBlock(kind="ancestry_5", available=False, notes=["Ancestry comparison disabled."])
        haplogroups = HaplogroupResult(available=False, notes=["Haplogroup comparison disabled."])
        disease = ComparisonBlock(kind="disease", available=False, notes=["Disease comparison disabled."])
        additional = ComparisonBlock(kind="additional", available=False)

        def _run_hominin():
            if not do_hominin:
                return hominin
            return compare_hominin(index, panel=panel, settings=self.settings, alignment=alignment)

        def _run_caste():
            if not do_caste:
                return caste
            return compare_caste(index, panel=panel, settings=self.settings, alignment=alignment)

        def _run_populations():
            if not do_pops:
                return populations
            return compare_populations(index, panel=panel, settings=self.settings, alignment=alignment)

        def _run_ancestry():
            if not do_ancestry:
                return ancestry
            return compare_ancestry(index, panel=panel, settings=self.settings, alignment=alignment)

        def _run_ancestry_5():
            if not do_ancestry:
                return ancestry_5
            return compare_ancestry_5source(index, panel=panel, settings=self.settings, alignment=alignment)

        def _run_haplogroups():
            if not do_haplo:
                return haplogroups
            return compare_haplogroups(index, settings=self.settings)

        def _run_disease():
            if not do_disease:
                return disease
            return compare_disease(index, settings=self.settings, by_rsid=by_rsid)

        def _run_additional():
            return compare_additional(
                index,
                settings=self.settings,
                other_index=other_index,
                other_filename=other_name,
            )

        with timings.stage("comparisons"):
            other_index: dict[tuple[str, int], dict] | None = None
            other_name: str | None = None
            relatedness = RelatednessResult(
                available=False,
                query_sample_id=summary.sample_id,
                notes=["No second VCF uploaded. Add a parent, relative, or any other SNP VCF to estimate closeness."],
            )

            if other_source is not None:
                relatedness, other_index, other_name = self._compare_other(
                    summary,
                    index,
                    other_source,
                    other_filename,
                    source_name,
                    query_by_rsid=by_rsid,
                )

            jobs = {
                "hominin": _run_hominin,
                "caste": _run_caste,
                "populations": _run_populations,
                "ancestry": _run_ancestry,
                "ancestry_5": _run_ancestry_5,
                "haplogroups": _run_haplogroups,
                "disease": _run_disease,
                "additional": _run_additional,
            }
            results: dict[str, object] = {}
            with ThreadPoolExecutor(max_workers=min(8, len(jobs))) as pool:
                futures = {name: pool.submit(fn) for name, fn in jobs.items()}
                for name, fut in futures.items():
                    results[name] = fut.result()
            hominin = results["hominin"]  # type: ignore[assignment]
            caste = results["caste"]  # type: ignore[assignment]
            populations = results["populations"]  # type: ignore[assignment]
            ancestry = results["ancestry"]  # type: ignore[assignment]
            ancestry_5 = results["ancestry_5"]  # type: ignore[assignment]
            haplogroups = results["haplogroups"]  # type: ignore[assignment]
            disease = results["disease"]  # type: ignore[assignment]
            additional = results["additional"]  # type: ignore[assignment]

        with timings.stage("community_ref"):
            community_ref = score_community_reference(
                ancestry if do_ancestry else None,
                haplogroups if do_haplo else None,
                ancestry_5=ancestry_5 if do_ancestry else None,
                caste=caste if do_caste else None,
                filename=source_name,
                settings=self.settings,
            )
        with timings.stage("snp_catalog"):
            publish_snp_list(summary, index, filename=source_name, settings=self.settings)

        note = assembly_note(summary)
        if note:
            hominin.notes = [note] + list(hominin.notes)
            caste.notes = [note] + list(caste.notes)
            populations.notes = [note] + list(populations.notes)
            ancestry.notes = [note] + list(ancestry.notes)
            ancestry_5.notes = [note] + list(ancestry_5.notes)
            haplogroups.notes = [note] + list(haplogroups.notes)
            relatedness.notes = [note] + list(relatedness.notes)
            community_ref.notes = [note] + list(community_ref.notes)
            disease.notes = [note] + list(disease.notes)
            additional.notes = [note] + list(additional.notes)

        timing_payload = timings.to_dict() if timings.enabled else None
        return AnalysisResult(
            ok=True,
            source_filename=source_name,
            vcf=summary,
            hominin=hominin,
            caste=caste,
            populations=populations,
            ancestry=ancestry,
            ancestry_5=ancestry_5,
            haplogroups=haplogroups,
            community_ref=community_ref,
            disease=disease,
            additional=additional,
            relatedness=relatedness,
            timings=timing_payload,
        )

    def _maybe_lift(self, summary, index: dict[tuple[str, int], dict]) -> dict[tuple[str, int], dict]:
        if summary.assembly != "GRCh38":
            return index
        chain_path = ensure_hg38_to_hg19_chain(
            self.settings.liftover_chain,
            download=self.settings.auto_download_chain,
        )
        if chain_path is None:
            return index
        try:
            if self._chain is None:
                self._chain = load_chain(chain_path)
            lifted, stats = lift_query_index(index, self._chain)
        except (OSError, ValueError):
            return index
        summary.lifted_to = "GRCh37"
        summary.n_lifted = int(stats["n_lifted"])
        summary.n_unmapped = int(stats["n_unmapped"])
        return lifted

    def _compare_other(
        self,
        summary,
        index: dict[tuple[str, int], dict],
        other_source: Path | str | TextIO | BinaryIO | None,
        other_filename: str | None,
        query_filename: str,
        query_by_rsid: dict[str, dict] | None = None,
    ) -> tuple[RelatednessResult, dict[tuple[str, int], dict] | None, str | None]:
        if other_source is None:
            return (
                RelatednessResult(
                    available=False,
                    query_sample_id=summary.sample_id,
                    notes=["No second VCF uploaded. Add a parent, relative, or any other SNP VCF to estimate closeness."],
                ),
                None,
                None,
            )
        other_name = other_filename or (
            Path(other_source).name if isinstance(other_source, (str, Path)) else "other.vcf"
        )
        try:
            other_summary, other_index = parse_genotype_file(
                other_source,
                filename=other_name,
                preview_limit=self.settings.variant_preview_limit,
            )
        except Raw23ParseError as exc:
            return (
                RelatednessResult(
                    available=False,
                    other_filename=other_name,
                    query_sample_id=summary.sample_id,
                    notes=[str(exc)],
                ),
                None,
                other_name,
            )
        except Exception as exc:  # noqa: BLE001
            return (
                RelatednessResult(
                    available=False,
                    other_filename=other_name,
                    query_sample_id=summary.sample_id,
                    notes=[f"Failed to parse second genotype file: {exc}"],
                ),
                None,
                other_name,
            )
        detect_assembly(other_summary, assume=self.settings.assume_assembly)
        other_index = self._maybe_lift(other_summary, other_index)
        other_by_rsid = build_rsid_index(other_index)
        return (
            compare_relatedness(
                index,
                other_index,
                other_filename=other_name,
                other_sample_id=other_summary.sample_id,
                query_sample_id=summary.sample_id,
                query_filename=query_filename,
                query_assembly=summary.assembly,
                other_assembly=other_summary.assembly,
                query_lifted_to=summary.lifted_to,
                other_lifted_to=other_summary.lifted_to,
                query_by_rsid=query_by_rsid,
                other_by_rsid=other_by_rsid,
            ),
            other_index,
            other_name,
        )
