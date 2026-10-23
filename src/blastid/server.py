"""로컬 전용 HTTPS 서버. Excel 작업창은 여기로만 말하고, NCBI 호출은 이 프로세스가 한다."""
from __future__ import annotations

import base64
import re
import threading
import uuid
from datetime import datetime
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from blastid import __version__
from blastid.fields import INPUT_FIELDS
from blastid.pipeline import (
    assess,
    blast_queries,
    clean_seq,
    combine_review,
    default_outdir,
    fasta_label,
    unique_labels,
)
from blastid.policy import DISCLAIMER_EN, DISCLAIMER_KO, describe_window
from blastid.policy import ET as ET_ZONE
from blastid.template import load_template

RID_RE = re.compile(r"RID\s+(\S+)")
URL_RE = re.compile(r"https://blast\.ncbi\.nlm\.nih\.gov/\S+")


def row_fill(template, level: str | None, backend: str | None) -> str | None:
    """템플릿 색. 로컬 행은 깃발이 없으면 local_fill 을 돌려준다."""
    style = template.style or {}
    if level == "strange":
        return str(style.get("strange_fill") or "F8CBAD")
    if level == "ambiguous":
        return str(style.get("ambiguous_fill") or "FFF2CC")
    if backend == "local":
        return str(style.get("local_fill") or "D6E3F0")
    return None


def addin_dir() -> Path:
    env = __import__("os").environ.get("BLASTID_ADDIN", "").strip()
    packaged = Path(__file__).resolve().parent / "addin"
    candidates = []
    if env:
        candidates.append(Path(env))
    candidates.append(packaged)
    for cand in candidates:
        if (cand / "taskpane.html").is_file():
            return cand
    raise RuntimeError(
        "작업창 파일 taskpane.html 이 설치본에 없습니다. "
        "blastid를 소스에서 다시 설치하거나 BLASTID_ADDIN 으로 addin 폴더를 지정하세요."
    )


class JobRequest(BaseModel):
    rows: list[dict] = Field(default_factory=list)
    bold_only: bool = False
    only_empty: bool = False
    min_ident: float = 90.0
    min_cover: float = 80.0
    alt_margin: float = 2.0
    batch_size: int = 10
    email: str = ""
    off_hours: bool = False
    template_yaml: str | None = None
    min_seq_len: int = 50
    dry_run: bool = False
    backend: str = "auto"
    remote_limit: int = 100
    local_db: str = "ITS_RefSeq_Fungi"
    max_wait: float = 60


class _Opt:
    def __init__(self, body: JobRequest):
        self.bold_only = body.bold_only
        self.only_empty = body.only_empty
        self.min_ident = body.min_ident
        self.min_cover = body.min_cover
        self.alt_margin = body.alt_margin
        self.batch_size = max(1, int(body.batch_size or 10))
        self.email = body.email or ""
        self.off_hours = body.off_hours
        self.template_yaml = body.template_yaml
        self.min_seq_len = body.min_seq_len
        self.dry_run = body.dry_run
        self.backend = body.backend or "auto"
        self.remote_limit = body.remote_limit
        self.local_db = body.local_db or "ITS_RefSeq_Fungi"
        self.max_wait = body.max_wait
        self.blastn = None
        self.reuse_rids = None


def _columns(template) -> list[dict]:
    return [
        {
            "field": col.field,
            "header": col.header,
            "aliases": list(col.aliases),
            "number_format": col.number_format,
            "role": col.role,
        }
        for col in template.columns
    ]


def identify_rows(rows: list[dict], options: _Opt, log, outdir: Path | None = None) -> dict:
    """시트에서 읽은 행을 BLAST하고, 머리글 기준으로 쓸 값을 돌려준다."""
    outdir = Path(outdir) if outdir else default_outdir(None)
    outdir.mkdir(parents=True, exist_ok=True)
    if options.template_yaml:
        yml = outdir / "upload-template.yaml"
        yml.write_text(options.template_yaml, encoding="utf-8")
        template = load_template(yml)
    else:
        template = load_template(None)
    queries = []
    for raw in rows:
        excel_row = int(raw.get("row") or 0)
        if options.bold_only and not raw.get("bold"):
            continue
        if options.only_empty and str(raw.get("name") or "").strip():
            continue
        sid = str(raw.get("id") or f"row{excel_row}").strip() or f"row{excel_row}"
        seq = clean_seq(raw.get("seq"), options.min_seq_len)
        if not seq:
            log(f"건너뜀 {excel_row}행 {sid}: 서열이 아니거나 {options.min_seq_len}bp 미만")
            continue
        queries.append(
            dict(
                key=len(queries),
                row=excel_row,
                label=fasta_label(sid, f"row{excel_row}"),
                sample_id=sid,
                seq=seq,
                merge=None,
                no=None,
            )
        )
    unique_labels(queries)
    payload = {"columns": _columns(template), "updates": [], "dry_run": bool(options.dry_run), "stopped": None}
    _ncbi, results, stopped = blast_queries(queries, options, outdir, log)
    payload["stopped"] = stopped
    if results is None:
        return payload
    shots = outdir / "shots"
    shots.mkdir(exist_ok=True)
    formats = {c.field: c.number_format for c in template.columns if c.field and c.number_format}
    for query in queries:
        if query["key"] not in results:
            continue
        result = results[query["key"]]
        level, comment = combine_review(
            *assess(result["hits"], options.min_ident, options.min_cover, options.alt_margin),
            None,
        )
        if result.get("backend") == "local":
            from blastid.localblast import local_note

            note = local_note(result.get("db") or "")
            comment = f"{comment} {note}".strip() if comment else note
        elif comment:
            comment = f"{comment} (RID {result['rid']})"
        from blastid.fields import build_field_values

        values = build_field_values(
            sample_id=query["sample_id"],
            sequence=query["seq"],
            blast=result,
            level=level,
            comment=comment,
            merge=None,
        )
        kept = {}
        for key, value in values.items():
            if key in INPUT_FIELDS or key == "screenshot" or value is None:
                continue
            kept[key] = value
        shot = None
        if level:
            from blastid.render import render_table

            png = shots / f"row{query['row']}.png"
            try:
                render_table(
                    png,
                    f"{query['row']}행 {query['sample_id']}  |  RID {result['rid']}  |  query {result['qlen']} bp",
                    result["hits"],
                )
                shot = base64.b64encode(png.read_bytes()).decode("ascii")
            except Exception as exc:
                log(f"  스크린샷 실패 {query['row']}행: {exc}")
        payload["updates"].append(
            {
                "row": query["row"],
                "backend": result.get("backend") or "remote",
                "level": level,
                "fill": row_fill(template, level, result.get("backend")),
                "values": kept,
                "formats": {k: formats[k] for k in kept if k in formats},
                "comment": comment,
                "screenshot_png": shot,
                "rid": result["rid"],
                "rid_url": result["rid_url"],
            }
        )
    return payload


class Job:
    def __init__(self):
        self.id = uuid.uuid4().hex[:12]
        self.status = "queued"
        self.events: list[dict] = []
        self.result: dict | None = None
        self.error: str | None = None
        self.lock = threading.Lock()

    def log(self, message: str) -> None:
        match = RID_RE.search(message or "")
        url = URL_RE.search(message or "")
        event = {"message": message, "rid": match.group(1) if match else None, "url": url.group(0) if url else None}
        with self.lock:
            self.events.append(event)
            if "대기" in message and self.status == "running":
                self.status = "waiting"
            elif self.status == "waiting" and "들어왔습니다" in message:
                self.status = "running"

    def snapshot(self) -> dict:
        with self.lock:
            return {
                "id": self.id,
                "status": self.status,
                "events": list(self.events),
                "result": self.result,
                "error": self.error,
            }


class JobStore:
    def __init__(self):
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()

    def add(self, job: Job) -> None:
        with self._lock:
            self._jobs[job.id] = job

    def get(self, job_id: str) -> Job | None:
        with self._lock:
            return self._jobs.get(job_id)


JOBS = JobStore()


def create_app() -> FastAPI:
    app = FastAPI(title="blastid", version=__version__)
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["https://localhost:8787", "https://127.0.0.1:8787"],
        allow_credentials=False,
        allow_methods=["GET", "POST"],
        allow_headers=["*"],
    )

    @app.get("/api/usage")
    def usage():
        from blastid.policy import UsageLog

        log = UsageLog()
        return {"count": log.count(), "limit": 100, "window": describe_window(datetime.now(ET_ZONE))}

    @app.get("/api/health")
    def health():
        return {
            "ok": True,
            "tool": "blastid",
            "version": __version__,
            "window": describe_window(datetime.now(ET_ZONE)),
        }

    @app.get("/api/disclaimer")
    def disclaimer():
        return {"ko": DISCLAIMER_KO, "en": DISCLAIMER_EN, "window": describe_window(datetime.now(ET_ZONE))}

    @app.get("/api/template")
    def default_template():
        path = load_template(None).source
        return PlainTextResponse(Path(path).read_text(encoding="utf-8"), media_type="text/yaml; charset=utf-8")

    @app.post("/api/jobs")
    def start_job(body: JobRequest):
        job = Job()
        JOBS.add(job)

        def run():
            job.status = "running"
            try:
                outdir = default_outdir(None)
                job.result = identify_rows(body.rows, _Opt(body), job.log, outdir)
                stopped = (job.result or {}).get("stopped")
                if stopped and not (job.result or {}).get("updates"):
                    job.status = "error"
                    job.error = stopped
                else:
                    job.status = "done"
                    if stopped:
                        job.error = stopped
                    else:
                        job.log("완료")
            except (Exception, SystemExit) as exc:
                job.status = "error"
                job.error = str(exc) or exc.__class__.__name__
                job.log(f"오류: {job.error}")

        threading.Thread(target=run, name=f"blastid-{job.id}", daemon=True).start()
        return {"id": job.id}

    @app.get("/api/jobs/{job_id}")
    def job_status(job_id: str):
        job = JOBS.get(job_id)
        if job is None:
            raise HTTPException(404, "작업을 찾지 못했습니다.")
        return job.snapshot()

    static = addin_dir()

    @app.get("/manifest.xml")
    def manifest():
        return FileResponse(static / "manifest.xml", media_type="application/xml")

    app.mount("/", StaticFiles(directory=str(static), html=False), name="addin")
    return app


def ensure_dev_cert(directory: Path | None = None) -> tuple[Path, Path]:
    """Office가 요구하는 HTTPS용 인증서. 공식 dev-certs가 있으면 그것을 쓴다."""
    office = Path.home() / ".office-addin-dev-certs"
    for folder in (office,):
        cert, key = folder / "localhost.crt", folder / "localhost.key"
        if cert.exists() and key.exists():
            return cert, key
    directory = directory or (Path.home() / ".blastid" / "dev-certs")
    cert, key = directory / "localhost.crt", directory / "localhost.key"
    if cert.exists() and key.exists():
        return cert, key
    directory.mkdir(parents=True, exist_ok=True)
    import subprocess

    subprocess.run(
        [
            "openssl", "req", "-x509", "-newkey", "rsa:2048",
            "-keyout", str(key), "-out", str(cert),
            "-days", "825", "-nodes",
            "-subj", "/CN=localhost",
            "-addext", "subjectAltName=DNS:localhost,IP:127.0.0.1",
        ],
        check=True,
        capture_output=True,
    )
    return cert, key


def serve(host: str = "127.0.0.1", port: int = 8787) -> None:
    import uvicorn

    if host not in {"127.0.0.1", "localhost"}:
        raise SystemExit("blastid serve 는 이 컴퓨터(127.0.0.1)에서만 엽니다.")
    try:
        addin_dir()
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    cert, key = ensure_dev_cert()
    print(f"blastid 작업창: https://localhost:{port}/taskpane.html")
    print(f"매니페스트: https://localhost:{port}/manifest.xml")
    print(f"인증서: {cert}")
    print("Excel에서 이 인증서를 신뢰한 뒤 src/blastid/addin/manifest.xml 을 사이드로드하세요.")
    print("pip로 설치한 경우에도 같은 파일이 패키지에 들어 있고, 없으면 서버가 여기서 멈춥니다.")
    print(DISCLAIMER_KO)
    uvicorn.run(
        create_app(),
        host="127.0.0.1",
        port=port,
        ssl_certfile=str(cert),
        ssl_keyfile=str(key),
        log_level="info",
    )
