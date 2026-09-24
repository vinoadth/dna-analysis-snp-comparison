from __future__ import annotations

import re

_RSID_RE = re.compile(r"(rs\d+)", re.I)


def normalize_rsid(value: object) -> str | None:
    raw = str(value or "").strip()
    if not raw or raw == ".":
        return None
    token = raw.split(";")[0].split(",")[0].strip()
    token = re.sub(r"^(gsa-|exm-|ilm-)", "", token, flags=re.I)
    match = _RSID_RE.search(token)
    if match:
        return match.group(1).lower()
    if token.lower().startswith("rs") and len(token) > 2:
        return token.lower()
    return None


def build_rsid_index(index: dict[tuple[str, int], dict]) -> dict[str, dict]:
    by_rsid: dict[str, dict] = {}
    for row in index.values():
        rsid = normalize_rsid(row.get("rsid"))
        if rsid:
            by_rsid[rsid] = row
    return by_rsid
