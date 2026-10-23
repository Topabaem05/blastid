"""로컬 NCBI BLAST+ (megablast). 원격 웹 검색과 같은 표시값 규칙을 쓴다.

IP를 바꾸거나 원격 한도를 우회하지 않는다. 한도를 넘기면 이 컴퓨터의 BLAST+만 쓴다.
"""
from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import tarfile
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from blastid.scoring import web_values

DEFAULT_DB = "ITS_RefSeq_Fungi"
BLAST_VERSION = "2.17.0+"
BLAST_FTP = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST"
SETUP_HINT = "blastid setup-local"


class LocalBlastUnavailable(RuntimeError):
    """BLAST+ 또는 로컬 DB가 없어 남은 검색을 진행할 수 없다."""

    def __init__(self, message: str, partial: dict | None = None):
        super().__init__(message)
        self.partial = partial or {}


def blast_home() -> Path:
    return Path.home() / ".blastid"


def config_path() -> Path:
    return blast_home() / "local.json"


def db_dir() -> Path:
    return blast_home() / "db"


def load_config() -> dict:
    path = config_path()
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def save_config(data: dict) -> None:
    path = config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")


def _which(name: str) -> str | None:
    found = shutil.which(name)
    return found


def find_blastn(explicit: str | None = None) -> str | None:
    if explicit:
        path = Path(explicit)
        if path.is_file():
            return str(path)
        found = _which(explicit)
        if found:
            return found
    cfg = load_config().get("blastn")
    if cfg and Path(cfg).is_file():
        return str(cfg)
    found = _which("blastn")
    if found:
        return found
    for cand in blast_home().glob("ncbi-blast*/bin/blastn"):
        if cand.is_file():
            return str(cand)
    nested = list(blast_home().glob("ncbi-blast*/*/bin/blastn"))
    if nested:
        return str(nested[0])
    return None


def find_makeblastdb(blastn: str | None = None) -> str | None:
    if blastn:
        sibling = Path(blastn).with_name("makeblastdb")
        if sibling.is_file():
            return str(sibling)
    return _which("makeblastdb")


def _db_exists(prefix: Path) -> bool:
    parent = prefix.parent
    if not parent.is_dir():
        return False
    stem = prefix.name
    for path in parent.glob(stem + "*"):
        if path.suffix.lower() in {".nsq", ".nin", ".ndb", ".nal", ".nhr"}:
            return True
    return False


def resolve_db(name: str | None) -> Path | None:
    """NCBI DB 이름 또는 makeblastdb 접두 경로."""
    name = (name or DEFAULT_DB).strip()
    candidates = [Path(name)]
    env = os.environ.get("BLASTDB", "").strip()
    if env:
        candidates.append(Path(env) / name)
    candidates.append(db_dir() / name)
    cfg_dir = load_config().get("db_dir")
    if cfg_dir:
        candidates.append(Path(cfg_dir) / name)
    for cand in candidates:
        if _db_exists(cand):
            return cand
    return None


def source_label(db: str, *, local: bool, rid: str | None = None) -> str:
    if local:
        return f"로컬 BLAST+ ({Path(db).name or db})"
    shown = db or "core_nt"
    return f"NCBI 웹 {shown} RID {rid or ''}".strip()


def local_note(db: str) -> str:
    label = Path(db).name or db
    return (
        f"로컬 BLAST+ ({label}) 결과입니다. "
        "웹 core_nt와 데이터베이스가 달라 1위 hit이 다를 수 있습니다."
    )


def missing_message(db: str, blastn: str | None) -> str:
    parts = []
    if not blastn:
        parts.append("blastn 실행 파일이 없습니다")
    if resolve_db(db) is None:
        parts.append(f"로컬 데이터베이스 {db} 가 없습니다")
    why = ", ".join(parts) or "로컬 BLAST+를 사용할 수 없습니다"
    return (
        f"{why}. 원격 24시간 한도를 넘기지 않도록 여기서 중단합니다. "
        f"설치: {SETUP_HINT} --db {db}  "
        "(공식 NCBI BLAST+와 update_blastdb.pl). "
        "이미 기록된 행은 결과 파일에 남겨 두었습니다. "
        "이어서 하려면 그 파일을 --only-empty 로 다시 실행하세요."
    )


def _sciname(defn: str) -> str:
    parts = (defn or "").split()
    if len(parts) >= 2 and parts[0][:1].isupper() and parts[1][:1].islower():
        if len(parts) >= 3 and parts[2].lower() in {"sp.", "aff.", "cf."}:
            return " ".join(parts[:3])
        return " ".join(parts[:2])
    return (defn or "").strip()


def _versioned(hit_id: str, accession: str) -> str:
    for text in (hit_id or "", accession or ""):
        token = text.split("|")[-2] if text.count("|") >= 2 else text
        token = token.split()[0]
        if "." in token and token.split(".")[-1].isdigit():
            return token
    return (accession or hit_id or "").split()[0]


def xml_to_results(xml_text: str, queries: list[dict], db: str) -> dict:
    """BLAST+ XML을 웹과 같은 web_values 결과 dict로 바꾼다."""
    root = ET.fromstring(xml_text)
    iterations = list(root.iter("Iteration"))
    if len(iterations) != len(queries):
        raise RuntimeError(f"로컬 BLAST 쿼리 {len(queries)}개인데 결과 {len(iterations)}개")
    label = Path(db).name or db
    out = {}
    for query, it in zip(queries, iterations):
        title = (it.findtext("Iteration_query-def") or "").split()
        if title and title[0] != query["label"].split()[0]:
            raise RuntimeError(f"로컬 BLAST 쿼리 순서 불일치 {title} != {query['label']}")
        qlen = int(it.findtext("Iteration_query-len") or len(query["seq"]))
        hits = []
        container = it.find("Iteration_hits")
        for hit in ([] if container is None else list(container.findall("Hit"))):
            hsps = []
            for hsp in hit.findall("Hit_hsps/Hsp"):
                hsps.append(
                    {
                        "bit_score": float(hsp.findtext("Hsp_bit-score") or 0),
                        "evalue": float(hsp.findtext("Hsp_evalue") or 0),
                        "identity": int(hsp.findtext("Hsp_identity") or 0),
                        "align_len": int(hsp.findtext("Hsp_align-len") or 0),
                        "query_from": int(hsp.findtext("Hsp_query-from") or 1),
                        "query_to": int(hsp.findtext("Hsp_query-to") or 1),
                    }
                )
            if not hsps:
                continue
            acc = _versioned(hit.findtext("Hit_id") or "", hit.findtext("Hit_accession") or "")
            defn = hit.findtext("Hit_def") or ""
            taxid = hit.findtext("Hit_taxid")
            desc = {
                "id": f"lcl|{acc}|",
                "accession": acc,
                "sciname": _sciname(defn),
                "taxid": int(taxid) if taxid and str(taxid).isdigit() else taxid or None,
                "title": defn,
            }
            hits.append({"len": int(hit.findtext("Hit_len") or 0), "description": [desc], "hsps": hsps})
        values = [web_values(h, qlen) for h in hits]
        out[query["key"]] = dict(
            rid=f"local:{label}",
            rid_url="",
            query_index=queries.index(query) + 1,
            qlen=qlen,
            seq_len=len(query["seq"]),
            db=label,
            local_db=str(db),
            backend="local",
            source=source_label(label, local=True),
            params={"expect": 0.05, "filter": "L;m;", "sc_match": 1, "sc_mismatch": -2, "gap": "0 0"},
            blast_version="BLAST+ local",
            hits=values,
            elapsed=0,
        )
    return out


def blastn_command(blastn: str, query: Path, db: Path, out: Path) -> list[str]:
    return [
        blastn,
        "-task", "megablast",
        "-query", str(query),
        "-db", str(db),
        "-word_size", "28",
        "-reward", "1",
        "-penalty", "-2",
        "-gapopen", "0",
        "-gapextend", "0",
        "-dust", "yes",
        "-soft_masking", "true",
        "-evalue", "0.05",
        "-max_target_seqs", "100",
        "-max_hsps", "100",
        "-outfmt", "5",
        "-out", str(out),
    ]


def run_local(queries: list[dict], db: str, outdir: Path, log=print, blastn: str | None = None) -> dict:
    program = find_blastn(blastn)
    resolved = resolve_db(db)
    if not program or resolved is None:
        raise LocalBlastUnavailable(missing_message(db, program))
    outdir = Path(outdir)
    raw = outdir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    fasta = raw / "local_query.fasta"
    fasta.write_text("".join(f">{q['label']}\n{q['seq']}\n" for q in queries), encoding="utf-8")
    xml_path = raw / "local.xml"
    env = os.environ.copy()
    env["BLASTDB"] = str(resolved.parent)
    cmd = blastn_command(program, fasta, resolved, xml_path)
    log(f"로컬 BLAST+ ({resolved.name}) {len(queries)}개 서열")
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=3600)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        detail = tail[-1] if tail else f"exit {proc.returncode}"
        raise RuntimeError(f"로컬 blastn 실패: {detail}")
    return xml_to_results(xml_path.read_text(encoding="utf-8", errors="replace"), queries, str(resolved))


def make_local_db(fasta: Path, out_prefix: Path, blastn: str | None = None) -> Path:
    """사용자 FASTA로 makeblastdb. 헤더는 ACC.1 학명 형태를 권장."""
    program = find_makeblastdb(find_blastn(blastn))
    if not program:
        raise LocalBlastUnavailable(missing_message(str(out_prefix), find_blastn(blastn)))
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [program, "-in", str(fasta), "-dbtype", "nucl", "-parse_seqids", "-out", str(out_prefix)],
        capture_output=True,
        text=True,
        timeout=3600,
    )
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        raise RuntimeError(f"makeblastdb 실패: {tail[-1] if tail else proc.returncode}")
    return out_prefix


def _blast_archive_name() -> str:
    system = platform.system()
    machine = platform.machine().lower()
    if system == "Linux" and machine in {"x86_64", "amd64"}:
        return f"ncbi-blast-{BLAST_VERSION}-x64-linux.tar.gz"
    if system == "Linux" and machine in {"aarch64", "arm64"}:
        return f"ncbi-blast-{BLAST_VERSION}-aarch64-linux.tar.gz"
    if system == "Darwin" and machine in {"arm64", "aarch64"}:
        return f"ncbi-blast-{BLAST_VERSION}-aarch64-macosx.tar.gz"
    if system == "Darwin":
        return f"ncbi-blast-{BLAST_VERSION}-x64-macosx.tar.gz"
    if system == "Windows":
        return f"ncbi-blast-{BLAST_VERSION}-x64-win64.tar.gz"
    raise SystemExit(f"이 플랫폼용 공식 BLAST+ 이름을 모릅니다: {system} {machine}")


def ensure_blast_binaries(log=print) -> str:
    found = find_blastn(None)
    if found and find_makeblastdb(found):
        log(f"blastn: {found}")
        return found
    name = _blast_archive_name()
    url = f"{BLAST_FTP}/{name}"
    dest = blast_home() / "ncbi-blast"
    dest.mkdir(parents=True, exist_ok=True)
    archive = dest / name
    log(f"공식 NCBI BLAST+ 받는 중: {url}")
    urllib.request.urlretrieve(url, archive)
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(dest)
    found = find_blastn(None)
    if not found:
        raise SystemExit("압축을 풀었지만 blastn을 찾지 못했습니다.")
    cfg = load_config()
    cfg["blastn"] = found
    cfg["db_dir"] = str(db_dir())
    save_config(cfg)
    log(f"blastn 설치: {found}")
    return found


def setup_local(db: str = DEFAULT_DB, fasta: str | None = None, name: str | None = None, taxdb: bool = True, log=print) -> None:
    """공식 BLAST+와 NCBI update_blastdb.pl, 또는 사용자 FASTA로 로컬 DB를 준비한다."""
    blastn = ensure_blast_binaries(log)
    folder = db_dir()
    folder.mkdir(parents=True, exist_ok=True)
    if fasta:
        label = name or Path(fasta).stem
        prefix = make_local_db(Path(fasta), folder / label, blastn)
        log(f"사용자 데이터베이스: {prefix}")
        cfg = load_config()
        cfg["db"] = str(prefix)
        cfg["db_dir"] = str(folder)
        save_config(cfg)
        return
    script = Path(blastn).with_name("update_blastdb.pl")
    if not script.is_file():
        raise SystemExit(f"update_blastdb.pl 이 blastn 옆에 없습니다: {script}. {SETUP_HINT} 로 BLAST+를 설치하세요.")
    perl = shutil.which("perl")
    if not perl:
        raise SystemExit("update_blastdb.pl 을 실행하려면 perl이 필요합니다.")
    targets = [db]
    if taxdb and db != "taxdb":
        targets.append("taxdb")
    for target in targets:
        log(f"NCBI 데이터베이스 받는 중: {target} → {folder}")
        proc = subprocess.run(
            [perl, str(script), "--decompress", target],
            cwd=folder,
            capture_output=True,
            text=True,
            timeout=7200,
        )
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()
            raise SystemExit(f"{target} 다운로드 실패: {tail[-1] if tail else proc.returncode}")
        log(f"{target} 준비됨")
    cfg = load_config()
    cfg["db"] = db
    cfg["db_dir"] = str(folder)
    cfg["blastn"] = blastn
    save_config(cfg)
