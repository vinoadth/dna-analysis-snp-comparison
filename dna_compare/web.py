from __future__ import annotations

import io
import json
from email.parser import BytesParser
from email.policy import HTTP
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse
from webbrowser import open as open_browser

from dna_compare.config import SAMPLE_DIR, UPLOAD_DIR, default_settings
from dna_compare.report import STATIC_DIR
from dna_compare.service import AnalysisService
from dna_compare.snp_catalog import catalog_path, filter_catalog_tsv, query_snp_catalog

MIME = {
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".txt": "text/plain; charset=utf-8",
}


class FormPart:
    def __init__(self, value: bytes, filename: str | None = None) -> None:
        self.value = value
        self.filename = filename
        self.file = io.BytesIO(value)


class MultipartForm:
    """cgi.FieldStorage stand-in; cgi was removed from the stdlib in Python 3.13."""

    def __init__(self, parts: dict[str, FormPart]) -> None:
        self._parts = parts

    def __contains__(self, name: str) -> bool:
        return name in self._parts

    def __getitem__(self, name: str) -> FormPart:
        return self._parts[name]

    def getfirst(self, name: str, default: str | None = None) -> str | None:
        part = self._parts.get(name)
        if part is None:
            return default
        return part.value.decode("utf-8", errors="replace")


def parse_multipart(headers, body: bytes) -> MultipartForm:
    content_type = headers.get("Content-Type") or ""
    envelope = b"Content-Type: " + content_type.encode("ascii", "surrogateescape") + b"\r\n\r\n" + body
    msg = BytesParser(policy=HTTP).parsebytes(envelope)
    parts: dict[str, FormPart] = {}
    if not msg.is_multipart():
        return MultipartForm(parts)
    for part in msg.iter_parts():
        if part.get_content_disposition() != "form-data":
            continue
        name = part.get_param("name", header="content-disposition")
        if isinstance(name, tuple):
            name = name[-1]
        if not name:
            continue
        payload = part.get_payload(decode=True)
        if payload is None:
            payload = b""
        elif isinstance(payload, str):
            payload = payload.encode("utf-8")
        filename = part.get_filename()
        parts[str(name)] = FormPart(payload, filename)
    return MultipartForm(parts)


def _truthy(value: str | None, default: bool = True) -> bool:
    if value is None:
        return default
    return value.strip().lower() not in {"0", "false", "no", "off"}


def _sample_path(name: str | None) -> Path | None:
    if not name:
        return None
    safe = Path(name).name
    path = SAMPLE_DIR / safe
    if not safe or not path.is_file() or path.resolve().parent != SAMPLE_DIR.resolve():
        return None
    return path


def _save_upload(filename: str, data: bytes) -> Path:
    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    dest = UPLOAD_DIR / Path(filename).name
    dest.write_bytes(data)
    return dest


def _static_file(url_path: str) -> Path | None:
    rel = url_path[len("/static/") :].lstrip("/")
    if not rel or ".." in Path(rel).parts:
        return None
    path = (STATIC_DIR / rel).resolve()
    root = STATIC_DIR.resolve()
    if path.is_file() and (path == root or root in path.parents):
        return path
    return None


class AnalyzeHandler(BaseHTTPRequestHandler):
    server_version = "DnaCompare/1.0"

    def log_message(self, fmt: str, *args) -> None:
        print(f"{self.address_string()} - {fmt % args}")

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path in {"/", "/index.html"}:
            self._send_file(STATIC_DIR / "index.html")
            return
        if parsed.path.startswith("/static/"):
            path = _static_file(parsed.path)
            if path is not None:
                self._send_file(path)
                return
            self._json(404, {"error": "not found"})
            return
        if parsed.path == "/api/samples":
            names = sorted(
                p.name
                for p in SAMPLE_DIR.iterdir()
                if p.is_file() and (p.suffix == ".vcf" or p.name.endswith(".vcf.gz"))
            )
            self._json(200, {"samples": names})
            return
        if parsed.path in {"/api/snps", "/api/snps.tsv"}:
            self._send_snp_list(parsed)
            return
        self._json(404, {"error": "not found"})

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        length = int(self.headers.get("Content-Length") or 0)
        body = self.rfile.read(length) if length else b""
        if parsed.path == "/api/analyze":
            query = parse_qs(parsed.query)
            ctype = self.headers.get("Content-Type") or ""
            if ctype.lower().startswith("multipart/form-data"):
                form = parse_multipart(self.headers, body)
                flags = {
                    "compare_hominin_flag": _truthy(form.getfirst("hominin", "true")),
                    "compare_caste_flag": _truthy(form.getfirst("caste", "true")),
                    "compare_populations_flag": _truthy(form.getfirst("populations", "true")),
                    "compare_ancestry_flag": _truthy(form.getfirst("ancestry", "true")),
                    "compare_haplogroups_flag": _truthy(form.getfirst("haplogroups", "true")),
                    "compare_disease_flag": _truthy(form.getfirst("disease", "true")),
                }
                primary_item = form["vcf"] if "vcf" in form else None
                other_item = form["other"] if "other" in form else None
                sample = _sample_path(form.getfirst("sample"))
                other_sample = _sample_path(form.getfirst("other_sample"))
                if primary_item is not None and getattr(primary_item, "filename", None) and primary_item.file:
                    dest = _save_upload(primary_item.filename, primary_item.file.read())
                    filename = Path(primary_item.filename).name
                    source = dest
                elif sample is not None:
                    source = sample
                    filename = sample.name
                else:
                    self._json(400, {"error": "Choose a bundled sample or upload a SNP VCF / 23andMe raw .txt file."})
                    return
                other_source = None
                other_name = None
                if other_item is not None and getattr(other_item, "filename", None) and other_item.file:
                    other_dest = _save_upload(other_item.filename, other_item.file.read())
                    other_source = other_dest
                    other_name = Path(other_item.filename).name
                elif other_sample is not None:
                    other_source = other_sample
                    other_name = other_sample.name
                settings = default_settings()
                settings.auto_download_pgs = True
                settings.analyze_timings = True
                payload = AnalysisService(settings).analyze(
                    source,
                    filename=filename,
                    other_source=other_source,
                    other_filename=other_name,
                    **flags,
                ).to_dict()
                self._json(200 if payload["ok"] else 400, payload)
                return
            filename = Path(self.headers.get("X-Filename") or "upload.vcf").name
            flags = {
                "compare_hominin_flag": _truthy((query.get("hominin") or ["true"])[0]),
                "compare_caste_flag": _truthy((query.get("caste") or ["true"])[0]),
                "compare_populations_flag": _truthy((query.get("populations") or ["true"])[0]),
                "compare_ancestry_flag": _truthy((query.get("ancestry") or ["true"])[0]),
                "compare_haplogroups_flag": _truthy((query.get("haplogroups") or ["true"])[0]),
                "compare_disease_flag": _truthy((query.get("disease") or ["true"])[0]),
            }
            dest = _save_upload(filename, body)
            settings = default_settings()
            settings.auto_download_pgs = True
            settings.analyze_timings = True
            payload = AnalysisService(settings).analyze(dest, filename=filename, **flags).to_dict()
            self._json(200 if payload["ok"] else 400, payload)
            return
        if parsed.path == "/api/analyze-sample":
            try:
                data = json.loads(body.decode("utf-8") or "{}")
            except json.JSONDecodeError:
                self._json(400, {"error": "invalid json"})
                return
            path = _sample_path(str(data.get("filename") or ""))
            if path is None:
                self._json(400, {"error": "unknown sample"})
                return
            other_path = _sample_path(str(data.get("other_filename") or ""))
            flags = {
                "compare_hominin_flag": _truthy(str(data.get("hominin", True))),
                "compare_caste_flag": _truthy(str(data.get("caste", True))),
                "compare_populations_flag": _truthy(str(data.get("populations", True))),
                "compare_ancestry_flag": _truthy(str(data.get("ancestry", True))),
                "compare_haplogroups_flag": _truthy(str(data.get("haplogroups", True))),
                "compare_disease_flag": _truthy(str(data.get("disease", True))),
            }
            settings = default_settings()
            settings.auto_download_pgs = True
            settings.analyze_timings = True
            payload = AnalysisService(settings).analyze(
                path,
                filename=path.name,
                other_source=other_path,
                other_filename=None if other_path is None else other_path.name,
                **flags,
            ).to_dict()
            self._json(200 if payload["ok"] else 400, payload)
            return
        self._json(404, {"error": "not found"})

    def _send_snp_list(self, parsed) -> None:
        query = parse_qs(parsed.query)
        filename = (query.get("file") or query.get("catalog") or [""])[0]
        if not filename:
            self._json(400, {"error": "file query parameter is required"})
            return
        path = catalog_path(filename, default_settings())
        called_raw = (query.get("called_only") or ["false"])[0]
        called_only = str(called_raw).strip().lower() in {"1", "true", "yes", "on"}
        if parsed.path.endswith(".tsv"):
            if not path.is_file():
                self._json(404, {"error": "SNP list not found. Analyze a VCF first."})
                return
            data = filter_catalog_tsv(
                path,
                called_only=called_only,
                chrom=(query.get("chrom") or [None])[0],
            )
            name = path.name
            if called_only:
                name = path.stem + ".called.tsv"
            self.send_response(200)
            self.send_header("Content-Type", "text/tab-separated-values; charset=utf-8")
            self.send_header(
                "Content-Disposition",
                f'attachment; filename="{name}"',
            )
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        payload = query_snp_catalog(
            path,
            chrom=(query.get("chrom") or [None])[0],
            q=(query.get("q") or [None])[0],
            called_only=called_only,
            offset=int((query.get("offset") or ["0"])[0] or 0),
            limit=int((query.get("limit") or ["200"])[0] or 200),
        )
        self._json(200 if payload.get("ok") else 404, payload)

    def _send_file(self, path: Path) -> None:
        data = path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", MIME.get(path.suffix, "application/octet-stream"))
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _json(self, status: int, payload: dict) -> None:
        data = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)


def serve(host: str = "127.0.0.1", port: int = 8765, open_browser_tab: bool = True) -> None:
    httpd = ThreadingHTTPServer((host, port), AnalyzeHandler)
    url = f"http://{host}:{port}/"
    print(f"Dashboard at {url}", flush=True)
    if open_browser_tab:
        open_browser(url)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping.")
        httpd.shutdown()
