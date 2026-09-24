from __future__ import annotations

import mmap
import threading
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from dna_compare._native import build_population_freqs_single_pass_native
from dna_compare.config import CODE_TO_CHROM, Settings, default_settings

HEADER_SIZE = 48
PACKED_MISSING = 3


@dataclass
class IndRecord:
    index: int
    sample_id: str
    sex: str
    population: str


@dataclass
class SnpRecord:
    index: int
    rsid: str
    chrom: str
    gpos: float
    pos: int
    allele1: str
    allele2: str


class PackedTGeno:
    """Random-access reader for AADR TGENO (individual-major packed ancestrymap)."""

    def __init__(self, geno_path: Path):
        self.path = Path(geno_path)
        self._file = self.path.open("rb")
        header = self._file.read(HEADER_SIZE)
        if not header.startswith(b"TGENO"):
            self._file.close()
            raise ValueError(f"{self.path} is not a TGENO packed file (got {header[:16]!r})")
        parts = header.split()
        self.nind = int(parts[1])
        self.nsnp = int(parts[2])
        self.bytes_per_ind = (self.nsnp + 3) // 4
        expected = HEADER_SIZE + self.bytes_per_ind * self.nind
        size = self.path.stat().st_size
        if size != expected:
            self._file.close()
            raise ValueError(
                f"TGENO size mismatch for {self.path}: expected {expected} bytes, got {size}"
            )
        self._mmap = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
        self._shifts = np.array([6, 4, 2, 0], dtype=np.uint8)

    def close(self) -> None:
        if getattr(self, "_mmap", None) is not None:
            self._mmap.close()
            self._mmap = None
        if getattr(self, "_file", None) is not None:
            self._file.close()
            self._file = None

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:  # noqa: BLE001
            pass

    def read_individual(self, ind_index: int) -> np.ndarray:
        if ind_index < 0 or ind_index >= self.nind:
            raise IndexError(ind_index)
        offset = HEADER_SIZE + ind_index * self.bytes_per_ind
        packed = np.frombuffer(
            self._mmap,
            dtype=np.uint8,
            count=self.bytes_per_ind,
            offset=offset,
        ).copy()
        unpacked = ((packed[:, None] >> self._shifts) & 3).ravel()
        return unpacked[: self.nsnp].astype(np.uint8)


class AadrPanel:
    def __init__(self, settings: Settings | None = None):
        self.settings = settings or default_settings()
        self._inds: list[IndRecord] | None = None
        self._snps: list[SnpRecord] | None = None
        self._id_to_ind: dict[str, IndRecord] | None = None
        self._pop_to_inds: dict[str, list[IndRecord]] | None = None
        self._geno: PackedTGeno | None = None
        self._cache_lock = threading.Lock()

    @property
    def available(self) -> bool:
        return (
            self.settings.aadr_geno.exists()
            and self.settings.aadr_ind.exists()
            and self.settings.aadr_snp.exists()
        )

    def inds(self) -> list[IndRecord]:
        if self._inds is None:
            records = []
            with self.settings.aadr_ind.open() as handle:
                for i, line in enumerate(handle):
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    records.append(IndRecord(i, parts[0], parts[1], parts[2]))
            self._inds = records
            self._id_to_ind = {rec.sample_id: rec for rec in records}
            pops: dict[str, list[IndRecord]] = {}
            for rec in records:
                pops.setdefault(rec.population, []).append(rec)
            self._pop_to_inds = pops
        return self._inds

    def snps(self) -> list[SnpRecord]:
        if self._snps is None:
            records = []
            with self.settings.aadr_snp.open() as handle:
                for i, line in enumerate(handle):
                    parts = line.split()
                    if len(parts) < 6:
                        continue
                    chrom_raw = parts[1]
                    chrom = CODE_TO_CHROM.get(int(chrom_raw), chrom_raw) if chrom_raw.isdigit() else chrom_raw
                    records.append(
                        SnpRecord(
                            index=i,
                            rsid=parts[0],
                            chrom=str(chrom),
                            gpos=float(parts[2]),
                            pos=int(parts[3]),
                            allele1=parts[4],
                            allele2=parts[5],
                        )
                    )
            self._snps = records
        return self._snps

    def sample_by_id(self, sample_id: str) -> IndRecord | None:
        self.inds()
        assert self._id_to_ind is not None
        return self._id_to_ind.get(sample_id)

    def samples_for_population(self, population: str, limit: int | None = None) -> list[IndRecord]:
        self.inds()
        assert self._pop_to_inds is not None
        recs = self._pop_to_inds.get(population, [])
        cap = self.settings.max_samples_per_pop if limit is None else limit
        return recs[:cap]

    def geno(self) -> PackedTGeno:
        if self._geno is None:
            self._geno = PackedTGeno(self.settings.aadr_geno)
        return self._geno

    def dosage_allele1(self, sample_id: str) -> np.ndarray:
        rec = self.sample_by_id(sample_id)
        if rec is None:
            raise KeyError(sample_id)
        g = self.geno().read_individual(rec.index).astype(np.float32)
        g[g == PACKED_MISSING] = np.nan
        return g

    def _population_allele1_freq_python(self, populations: tuple[str, ...]) -> np.ndarray:
        dosages = []
        for pop in populations:
            for rec in self.samples_for_population(pop):
                dosages.append(self.dosage_allele1(rec.sample_id))
        if not dosages:
            raise KeyError(f"No AADR samples for {populations}")
        stacked = np.vstack(dosages)
        with np.errstate(all="ignore"):
            return np.nanmean(stacked / 2.0, axis=0).astype(np.float32)

    @staticmethod
    def freq_cache_path(settings: Settings, populations: tuple[str, ...]) -> Path:
        cache_name = "aadr_freq_" + "_".join(sorted(populations)) + ".npy"
        return settings.cache_dir / cache_name

    def _freq_cache_path(self, populations: tuple[str, ...]) -> Path:
        return self.freq_cache_path(self.settings, populations)

    def _single_pass_membership(
        self, groups: list[tuple[str, ...]]
    ) -> tuple[list[int], list[int], list[int], dict[int, list[int]]]:
        ind_to_groups: dict[int, list[int]] = {}
        for gi, group in enumerate(groups):
            for pop in group:
                for rec in self.samples_for_population(pop):
                    ind_to_groups.setdefault(rec.index, []).append(gi)
        if not ind_to_groups:
            missing = ", ".join("_".join(sorted(g)) for g in groups[:3])
            raise KeyError(f"No AADR samples for population groups (e.g. {missing})")

        ind_indices = sorted(ind_to_groups)
        ind_group_offsets = [0]
        ind_group_ids: list[int] = []
        for ind_index in ind_indices:
            for gi in ind_to_groups[ind_index]:
                ind_group_ids.append(gi)
            ind_group_offsets.append(len(ind_group_ids))
        return ind_indices, ind_group_offsets, ind_group_ids, ind_to_groups

    def _build_freq_caches_single_pass_python(
        self,
        groups: list[tuple[str, ...]],
        ind_indices: list[int],
        ind_to_groups: dict[int, list[int]],
    ) -> list[np.ndarray]:
        nsnp = self.geno().nsnp
        n_groups = len(groups)
        sums = [np.zeros(nsnp, dtype=np.float64) for _ in range(n_groups)]
        counts = [np.zeros(nsnp, dtype=np.uint32) for _ in range(n_groups)]
        geno = self.geno()
        for ind_index in ind_indices:
            packed = geno.read_individual(ind_index)
            valid = packed != PACKED_MISSING
            dose = packed.astype(np.float64) / 2.0
            for gi in ind_to_groups[ind_index]:
                sums[gi][valid] += dose[valid]
                counts[gi][valid] += 1
        freqs: list[np.ndarray] = []
        for gi in range(n_groups):
            freq = np.full(nsnp, np.nan, dtype=np.float32)
            ok = counts[gi] > 0
            freq[ok] = (sums[gi][ok] / counts[gi][ok]).astype(np.float32)
            freqs.append(freq)
        return freqs

    def _build_freq_caches_single_pass(self, groups: list[tuple[str, ...]]) -> int:
        """Decode each needed AADR individual once and write all missing .npy caches."""
        if not groups:
            return 0
        pending = [tuple(sorted(group)) for group in groups if group]
        pending = [group for group in pending if not self._freq_cache_path(group).exists()]
        if not pending:
            return 0

        ind_indices, ind_group_offsets, ind_group_ids, ind_to_groups = self._single_pass_membership(
            pending
        )
        nsnp = self.geno().nsnp
        freqs = build_population_freqs_single_pass_native(
            str(self.settings.aadr_geno),
            nsnp,
            ind_indices,
            ind_group_offsets,
            ind_group_ids,
            len(pending),
        )
        if freqs is None:
            freqs = self._build_freq_caches_single_pass_python(pending, ind_indices, ind_to_groups)

        self.settings.cache_dir.mkdir(parents=True, exist_ok=True)
        for group, freq in zip(pending, freqs, strict=True):
            np.save(self._freq_cache_path(group), freq)
        return len(pending)

    def build_all_population_freq_caches(
        self,
        population_groups: list[tuple[str, ...]] | None = None,
        *,
        force: bool = False,
    ) -> tuple[int, int]:
        """Build AADR frequency caches in one TGENO pass. Returns (built, skipped)."""
        from dna_compare.population_groups import collect_population_groups

        groups = population_groups if population_groups is not None else collect_population_groups(
            self.settings
        )
        unique = {tuple(sorted(group)) for group in groups if group}
        if force:
            for group in unique:
                path = self._freq_cache_path(group)
                if path.exists():
                    path.unlink()
        with self._cache_lock:
            missing = [group for group in unique if not self._freq_cache_path(group).exists()]
            skipped = len(unique) - len(missing)
            built = self._build_freq_caches_single_pass(missing)
        return built, skipped

    def population_allele1_freq(self, populations: tuple[str, ...] | str) -> np.ndarray:
        if isinstance(populations, str):
            populations = (populations,)
        populations = tuple(sorted(populations))
        cache_path = self._freq_cache_path(populations)
        if cache_path.exists():
            return np.load(cache_path)

        with self._cache_lock:
            if cache_path.exists():
                return np.load(cache_path)
            self._build_freq_caches_single_pass([populations])
            if not cache_path.exists():
                raise KeyError(f"No AADR samples for {populations}")
            return np.load(cache_path)

    def prefetch_population_freqs(self, population_groups: list[tuple[str, ...]]) -> None:
        unique = {tuple(sorted(group)) for group in population_groups if group}
        missing = [group for group in unique if not self._freq_cache_path(group).exists()]
        if not missing:
            return
        with self._cache_lock:
            missing = [group for group in unique if not self._freq_cache_path(group).exists()]
            if missing:
                self._build_freq_caches_single_pass(missing)
