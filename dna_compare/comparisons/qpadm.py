from __future__ import annotations

import numpy as np


def f4(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    """Patterson f4(A,B;C,D) = mean((pA-pB)*(pC-pD)) on aligned allele frequencies."""
    return float(np.mean((a - b) * (c - d)))


def qpadm_weights(
    target: np.ndarray,
    sources: list[np.ndarray],
    rights: list[np.ndarray],
) -> tuple[np.ndarray, np.ndarray, float]:
    """qpAdm-style weights: target ≈ sum w_i source_i, sum(w)=1, using f4 vs outgroups.

    Returns (raw weights, clipped-renormalized weights, residual RMS of the f4 fit).
    """
    if len(sources) < 2:
        raise ValueError("qpAdm needs at least two sources")
    if len(rights) < 3:
        raise ValueError("qpAdm needs at least three right/outgroup populations")

    k = len(sources)
    s0 = sources[0]
    y_rows: list[float] = []
    x_rows: list[list[float]] = []
    for i, ra in enumerate(rights):
        for rb in rights[i + 1 :]:
            y_rows.append(f4(target, s0, ra, rb))
            x_rows.append([f4(sources[j], s0, ra, rb) for j in range(1, k)])
    y = np.asarray(y_rows, dtype=float)
    x = np.asarray(x_rows, dtype=float)
    rest, residuals, *_ = np.linalg.lstsq(x, y, rcond=None)
    raw = np.concatenate([[1.0 - float(rest.sum())], rest.astype(float)])
    fitted = x @ rest
    rss = float(np.mean((y - fitted) ** 2))
    clip = np.clip(raw, 0.0, None)
    total = float(clip.sum())
    normed = clip / total if total > 0 else np.zeros_like(clip)
    return raw, normed, rss
