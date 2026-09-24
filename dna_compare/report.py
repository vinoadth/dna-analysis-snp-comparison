from __future__ import annotations

import json
from pathlib import Path

STATIC_DIR = Path(__file__).resolve().parent / "static"


def render_report_html(payload: dict) -> str:
    """Standalone HTML report with the same charts as the live dashboard."""
    bootstrap = (STATIC_DIR / "vendor" / "cerulean.min.css").read_text(encoding="utf-8")
    css = (STATIC_DIR / "app.css").read_text(encoding="utf-8")
    js = (STATIC_DIR / "app.js").read_text(encoding="utf-8")
    data = json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")
    title = payload.get("source_filename") or "SNP comparison"
    return f"""<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>{_escape(title)}</title>
    <style>{bootstrap}</style>
    <style>{css}</style>
  </head>
  <body>
    <nav class="navbar navbar-dark bg-primary">
      <div class="container-fluid px-3 px-xl-4">
        <span class="navbar-brand mb-0 h1">SNP comparison</span>
      </div>
    </nav>
    <main class="container-fluid px-3 px-xl-4 py-4">
      <p class="lead text-secondary page-intro">Embedded analysis of {_escape(str(title))}. Mixture weights use overlapping AADR Human Origins SNPs.</p>
      <div id="results" class="results-wide"></div>
    </main>
    <script>window.ANALYSIS_PAYLOAD = {data};</script>
    <script>{js}</script>
  </body>
</html>
"""


def write_report(payload: dict, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_report_html(payload), encoding="utf-8")
    return path


def _escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
