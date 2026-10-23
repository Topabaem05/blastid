"""NCBI BLAST URL API와 Taxonomy E-utilities.

이메일과 API 키는 환경변수 ``NCBI_EMAIL``, ``NCBI_API_KEY`` 가
비어 있지 않을 때만 요청에 넣는다. 없어도 동작한다.
"""
from __future__ import annotations

import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

from blastid import __version__

TOOL_NAME = "blastid"
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
WEB_LINK = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}"

# 웹 BLAST blastn 기본값(core_nt / megablast). RID report.params 로 확인한 값.
WEB_PARAMS = [
    ("PROGRAM", "blastn"),
    ("MEGABLAST", "on"),
    ("DATABASE", "core_nt"),
    ("EXPECT", "0.05"),
    ("WORD_SIZE", "28"),
    ("NUCL_REWARD", "1"),
    ("NUCL_PENALTY", "-2"),
    ("GAPCOSTS", "0 0"),
    ("FILTER", "L"),
    ("FILTER", "m"),
    ("HITLIST_SIZE", "100"),
    ("SHORT_QUERY_ADJUST", "on"),
]
SUBMIT_GAP = 10.0
FIRST_POLL = 20.0
POLL_GAP = 60.0
RETRY_HTTP = {429, 500, 502, 503, 504}


def ncbi_email() -> str:
    return os.environ.get("NCBI_EMAIL", "").strip()


def ncbi_api_key() -> str:
    return os.environ.get("NCBI_API_KEY", "").strip()


def blast_identity_params() -> list[tuple[str, str]]:
    """BLAST URL에 붙이는 tool/email. api_key 는 E-utilities 전용."""
    params = [("tool", TOOL_NAME)]
    email = ncbi_email()
    if email:
        params.append(("email", email))
    return params


class NCBI:
    def __init__(self, log=print):
        self.log = log
        self._last_put = 0.0
        self._last_poll: dict[str, float] = {}
        self.email = ncbi_email()
        self.api_key = ncbi_api_key()
        try:
            from Bio import Entrez
        except ImportError as e:
            raise SystemExit("Biopython이 필요합니다: pip install biopython") from e
        Entrez.tool = TOOL_NAME
        if self.email:
            Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
        self.Entrez = Entrez
        if not self.email:
            import warnings

            warnings.filterwarnings("ignore", message=".*[Ee]mail.*", module="Bio.Entrez")

    def _req(self, params, post=False, timeout=180, tries=4) -> str:
        params = list(params) + blast_identity_params()
        data = urllib.parse.urlencode(params)
        headers = {"User-Agent": f"{TOOL_NAME}/{__version__}"}
        for k in range(tries):
            try:
                if post:
                    req = urllib.request.Request(BLAST_URL, data=data.encode(), headers=headers)
                else:
                    req = urllib.request.Request(BLAST_URL + "?" + data, headers=headers)
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    return resp.read().decode("utf-8", "replace")
            except urllib.error.HTTPError as e:
                if e.code not in RETRY_HTTP or k == tries - 1:
                    body = e.read().decode("utf-8", "replace")[:500]
                    raise RuntimeError(f"NCBI HTTP {e.code}: {body}") from e
                self.log(f"  NCBI HTTP {e.code}; {30 * (k + 1)}초 후 재시도")
                time.sleep(30 * (k + 1))
            except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
                if k == tries - 1:
                    raise
                self.log(f"  네트워크 오류({e}); {30 * (k + 1)}초 후 재시도")
                time.sleep(30 * (k + 1))
        raise RuntimeError("NCBI 요청 실패")

    def submit(self, fasta_text: str) -> tuple[str, int]:
        wait = self._last_put + SUBMIT_GAP - time.time()
        if wait > 0:
            time.sleep(wait)
        t = self._req([("CMD", "Put"), ("QUERY", fasta_text)] + WEB_PARAMS, post=True)
        self._last_put = time.time()
        rid = re.search(r"^\s*RID = (\S+)", t, re.M)
        rtoe = re.search(r"^\s*RTOE = (\d+)", t, re.M)
        if not rid:
            raise RuntimeError("BLAST 제출 실패: " + re.sub(r"<[^>]+>", " ", t)[:500])
        return rid.group(1), int(rtoe.group(1)) if rtoe else 30

    def status(self, rid: str) -> str:
        t = self._req([("CMD", "Get"), ("FORMAT_OBJECT", "SearchInfo"), ("RID", rid)])
        self._last_poll[rid] = time.time()
        m = re.search(r"Status=(\w+)", t)
        return m.group(1) if m else "UNKNOWN"

    def get_json(self, rid: str) -> dict:
        t = self._req(
            [
                ("CMD", "Get"),
                ("RID", rid),
                ("FORMAT_TYPE", "JSON2_S"),
                ("FORMAT_OBJECT", "Alignment"),
                ("DESCRIPTIONS", "100"),
                ("ALIGNMENTS", "100"),
            ]
        )
        try:
            return json.loads(t)
        except json.JSONDecodeError as e:
            raise RuntimeError(f"RID {rid} JSON 해석 실패: {t[:300]}") from e

    def get_xml(self, rid: str) -> str:
        return self._req(
            [
                ("CMD", "Get"),
                ("RID", rid),
                ("FORMAT_TYPE", "XML"),
                ("DESCRIPTIONS", "100"),
                ("ALIGNMENTS", "100"),
            ]
        )

    def wait_all(self, jobs: dict[str, tuple[float, int]]) -> dict[str, float]:
        """jobs: {rid: (submit_time, rtoe)}. 모두 READY 될 때까지 확인."""
        pending = dict(jobs)
        done: dict[str, float] = {}
        while pending:
            now = time.time()
            nxt = []
            for rid, (t0, rtoe) in list(pending.items()):
                if rid not in self._last_poll:
                    due = t0 + (max(FIRST_POLL, rtoe) if rtoe else 0)
                else:
                    due = self._last_poll[rid] + POLL_GAP
                if due <= now:
                    st = self.status(rid)
                    self.log(f"  RID {rid}: {st} ({time.time() - t0:.0f}초 경과)")
                    if st == "READY":
                        done[rid] = time.time()
                        del pending[rid]
                    elif st in ("FAILED", "UNKNOWN"):
                        raise RuntimeError(f"RID {rid} 상태 {st} (NCBI 측 실패/만료)")
                    else:
                        nxt.append(self._last_poll[rid] + POLL_GAP)
                else:
                    nxt.append(due)
            if pending:
                time.sleep(max(1.0, min(nxt) - time.time()))
        return done

    def taxonomy(self, taxids) -> dict:
        taxids = sorted({str(t) for t in taxids if t})
        out: dict = {}
        if not taxids:
            return out
        for i in range(0, len(taxids), 100):
            chunk = taxids[i : i + 100]
            recs = None
            for k in range(4):
                try:
                    h = self.Entrez.efetch(db="taxonomy", id=",".join(chunk), retmode="xml")
                    recs = self.Entrez.read(h)
                    h.close()
                    break
                except Exception as e:
                    if k == 3:
                        raise
                    self.log(f"  taxonomy 재시도({e})")
                    time.sleep(5 * (k + 1))
            if isinstance(recs, dict):
                recs = [recs]
            for r in recs or []:
                ids = [str(r["TaxId"])] + [str(x) for x in r.get("AkaTaxIds", [])]
                info = {
                    "taxid": str(r["TaxId"]),
                    "name": str(r["ScientificName"]),
                    "lineage_full": str(r.get("Lineage", "")),
                }
                for x in ids:
                    out[x] = info
        return out


def rid_url(rid: str) -> str:
    return WEB_LINK.format(rid=rid)


def check_params(results: dict, log=print) -> None:
    for r in results.values():
        p = r.get("params") or {}
        if (
            r.get("db") != "core_nt"
            or p.get("filter") != "L;m;"
            or p.get("expect") != 0.05
            or p.get("sc_match") != 1
            or p.get("sc_mismatch") != -2
        ):
            log(f"  경고: RID {r.get('rid')} 설정이 웹 기본과 다름: db={r.get('db')} params={p}")


def blast_reports(payload: dict) -> list:
    reps = payload["BlastOutput2"]
    if isinstance(reps, dict):
        return [reps]
    return list(reps)
