"""Sanger .ab1 읽기, 품질 트리밍, 정·역방향 짝 맞추기, 겹침 합의 서열.

기본 백엔드는 Biopython PairwiseAligner 로 끝 갭을 공짜로 둔 overlap 정렬이다.
염기가 다르면 Phred 품질 차이가 ``conflict_q`` 이상일 때 높은 쪽을 고르고,
그렇지 않으면 IUPAC 모호 부호를 남긴다. 겹침이 너무 짧거나 일치도가 낮으면
평균 품질이 더 좋은 단일 리드로 돌아간다.

``--merger tracy`` 는 gear-genomics tracy ``consensus`` 가 PATH에 있을 때
그 합의 서열을 쓴다. tracy는 트레이스에서 염기를 다시 호출한다.
"""
from __future__ import annotations

import csv
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

GAP = set("-. ")
ACGT = set("ACGT")

# 1bp indel 은 갭으로, SNP 는 불일치로 남기도록 맞춘 점수.
MATCH = 2
MISMATCH = -5
GAP_OPEN = -3
GAP_EXTEND = -1

FORWARD_DIRS = {
    "f", "fwd", "fw", "for", "forward",
    "its1f", "its5", "its1", "its5f", "lr1",
}
REVERSE_DIRS = {
    "r", "rev", "rv", "reverse",
    "its4", "its4r", "its4b", "nl4", "lr3",
}
# 긴 토큰을 앞에 둔다 (its4r 가 its4 보다 먼저).
_DIR_ALT = (
    "forward|reverse|fwd|rev|for|fw|rv|"
    "its1f|its5f|its5|its1|its4r|its4b|its4|lr1|lr3|nl4|f|r"
)
DEFAULT_PAIR_REGEXES = [
    rf"^(?P<id>.+?)[_\-.](?P<dir>{_DIR_ALT})$",
]

IUPAC_SETS = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "S": set("CG"), "W": set("AT"),
    "K": set("GT"), "M": set("AC"),
    "B": set("CGT"), "D": set("AGT"), "H": set("ACT"), "V": set("ACG"),
    "N": set("ACGT"),
}
IUPAC_CODE = {frozenset(v): k for k, v in IUPAC_SETS.items()}


class TracyUnavailable(RuntimeError):
    pass


class TracyError(RuntimeError):
    pass


@dataclass
class SangerRead:
    path: Path
    sample_id: str
    seq: str
    quals: list[int]
    raw_len: int
    trim_start: int
    trim_end: int
    trim_mean_q: float
    trim_q20: float
    n_count: int
    no_quality: bool
    warning: str = ""
    direction: str = "U"
    stem: str = ""

    @property
    def trim_len(self) -> int:
        return len(self.seq)


@dataclass
class MergeResult:
    seq: str
    quals: list[int]
    mode: str  # consensus, fallback, single
    backend: str  # biopython, tracy, none
    overlap_len: int | None = None
    overlap_ident: float | None = None  # 0–1
    n_conflicts: int | None = None
    n_resolved: int | None = None
    n_ambiguities: int | None = None
    n_indels: int | None = None
    warning: str = ""
    direction: str = "U"
    trimmed_len: int = 0
    raw_len: int | None = None
    trim_mean_q: float | None = None
    sample_id: str = ""
    reads: list[SangerRead] = field(default_factory=list)

    def field_values(self) -> dict:
        return {
            "overlap_len": self.overlap_len,
            "overlap_ident": self.overlap_ident,
            "n_conflicts": self.n_conflicts,
            "n_resolved": self.n_resolved,
            "n_ambiguities": self.n_ambiguities,
            "n_indels": self.n_indels,
            "merge_mode": self.mode,
            "merge_backend": self.backend,
            "trimmed_len": self.trimmed_len,
            "raw_len": self.raw_len,
            "trim_mean_q": self.trim_mean_q,
            "merge_warning": self.warning or None,
            "direction": self.direction,
        }


def mott_trim(quals: list[int], cutoff: float = 0.05) -> tuple[int, int]:
    """Mott 알고리즘. (cutoff - 오류확률) 누적합이 최대인 구간을 반환. 끝은 exclusive."""
    best, cur, s, bs, be = 0.0, 0.0, 0, 0, 0
    for i, q in enumerate(quals):
        cur += cutoff - 10 ** (-q / 10.0)
        if cur <= 0:
            cur, s = 0.0, i + 1
        elif cur > best:
            best, bs, be = cur, s, i + 1
    return bs, be


def read_ab1(path: Path, cutoff: float = 0.05, trimmer: str = "mott") -> SangerRead:
    from Bio import SeqIO

    rec = SeqIO.read(str(path), "abi")
    q = list(rec.letter_annotations.get("phred_quality") or [])
    noq = (not q) or max(q) == 0
    raw = str(rec.seq).upper().replace("U", "T")
    if noq:
        s, e, seq = 0, len(raw), raw
    elif trimmer == "biopython":
        trimmed = SeqIO.read(str(path), "abi-trim")
        seq = str(trimmed.seq).upper().replace("U", "T")
        s = raw.find(seq)
        s = max(s, 0)
        e = s + len(seq)
    else:
        s, e = mott_trim(q, cutoff)
        seq = raw[s:e]
    tq = q[s:e] if q and not noq else [0] * len(seq)
    if len(tq) < len(seq):
        tq = tq + [0] * (len(seq) - len(tq))
    tq = tq[: len(seq)]
    return SangerRead(
        path=Path(path),
        sample_id=Path(path).stem,
        stem=Path(path).stem,
        seq=seq,
        quals=[int(x) for x in tq],
        raw_len=len(raw),
        trim_start=s + 1,
        trim_end=e,
        trim_mean_q=round(sum(tq) / len(tq), 1) if tq else 0.0,
        trim_q20=round(100 * sum(1 for x in tq if x >= 20) / len(tq), 1) if tq else 0.0,
        n_count=seq.count("N"),
        no_quality=noq,
        direction="U",
    )


def collect_ab1(paths: list[str | Path]) -> list[Path]:
    out: list[Path] = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            for f in sorted(p.rglob("*")):
                if f.is_file() and f.suffix.lower() in {".ab1", ".abi"}:
                    out.append(f)
        elif p.is_file():
            out.append(p)
    # 중복 제거, 순서 유지
    seen = set()
    uniq = []
    for f in out:
        key = str(f.resolve())
        if key not in seen:
            seen.add(key)
            uniq.append(f)
    return uniq


def _compile_patterns(patterns: list[str]) -> list[re.Pattern]:
    compiled = []
    for p in patterns:
        try:
            cre = re.compile(p, re.I)
        except re.error as e:
            raise SystemExit(f"짝 정규식이 잘못되었습니다: {p} ({e})") from e
        if "id" not in cre.groupindex or "dir" not in cre.groupindex:
            raise SystemExit(f"짝 정규식에 (?P<id>...) 와 (?P<dir>...) 가 필요합니다: {p}")
        compiled.append(cre)
    return compiled


def classify_direction(token: str) -> str:
    d = token.lower()
    if d in FORWARD_DIRS:
        return "F"
    if d in REVERSE_DIRS:
        return "R"
    return "U"


def match_pair_name(stem: str, patterns: list[re.Pattern]):
    for cre in patterns:
        m = cre.fullmatch(stem)
        if not m:
            continue
        return m.group("id"), classify_direction(m.group("dir"))
    return None


def _resolve_pair_path(raw: str, base: Path) -> Path:
    p = Path(raw.strip())
    if not p.is_absolute():
        cand = (base / p).resolve()
        if cand.exists():
            return cand
        cand2 = (Path.cwd() / p).resolve()
        if cand2.exists():
            return cand2
        return cand
    return p


def read_pairs_csv(path: Path) -> list[tuple[str, Path, Path]]:
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        if not reader.fieldnames:
            raise SystemExit(f"짝 CSV에 헤더가 없습니다: {path}")
        fields = {(n or "").strip().lower(): n for n in reader.fieldnames}
        def col(*names):
            for n in names:
                if n in fields:
                    return fields[n]
            return None
        c_f = col("forward", "fwd", "f", "r1", "정방향")
        c_r = col("reverse", "rev", "r", "r2", "역방향")
        c_id = col("id", "sample", "name", "시료")
        if not c_f or not c_r:
            raise SystemExit("짝 CSV에는 forward, reverse 열이 필요합니다.")
        base = path.parent
        for rec in reader:
            fp = (rec.get(c_f) or "").strip()
            rp = (rec.get(c_r) or "").strip()
            if not fp and not rp:
                continue
            if not fp or not rp:
                raise SystemExit(f"짝 CSV에서 한쪽 파일이 비어 있습니다: {rec}")
            sid = (rec.get(c_id) or "").strip() if c_id else ""
            rows.append((sid, _resolve_pair_path(fp, base), _resolve_pair_path(rp, base)))
    return rows


def pair_reads(
    files: list[Path],
    patterns: list[str] | None = None,
    pairs_csv: str | None = None,
    log=print,
    explicit: list[tuple[str, Path, Path]] | None = None,
) -> tuple[list[tuple[str, Path, Path]], list[Path]]:
    """반환: ([(id, forward, reverse), ...], [unpaired paths])."""
    patterns_c = _compile_patterns(patterns or DEFAULT_PAIR_REGEXES)
    used: set[str] = set()
    pairs: list[tuple[str, Path, Path]] = []
    listed = explicit
    if listed is None and pairs_csv:
        listed = read_pairs_csv(Path(pairs_csv))
    if listed:
        for sid, fp, rp in listed:
            if not fp.exists() or not rp.exists():
                raise SystemExit(f"짝 CSV의 파일을 찾을 수 없습니다: {fp} / {rp}")
            if fp.resolve() == rp.resolve():
                raise SystemExit(f"정방향과 역방향이 같은 파일입니다: {fp}")
            if not sid:
                sid = fp.stem
            pairs.append((sid, fp, rp))
            used.add(str(fp.resolve()))
            used.add(str(rp.resolve()))
            log(f"짝 CSV: {sid}  ←  {fp.name} + {rp.name}")
    groups: dict[str, dict[str, list[Path]]] = {}
    singles: list[Path] = []
    for f in files:
        if str(f.resolve()) in used:
            continue
        hit = match_pair_name(f.stem, patterns_c)
        if not hit or hit[1] == "U":
            singles.append(f)
            continue
        sid, direc = hit
        groups.setdefault(sid, {"F": [], "R": []})
        groups[sid][direc].append(f)
    for sid, g in groups.items():
        fs, rs = g["F"], g["R"]
        if len(fs) == 1 and len(rs) == 1:
            pairs.append((sid, fs[0], rs[0]))
            log(f"짝: {sid}  ←  {fs[0].name} + {rs[0].name}")
        else:
            if fs and rs:
                log(f"짝을 하나로 정할 수 없어 각각 처리: {sid} (F {len(fs)}개, R {len(rs)}개)")
            for f in fs + rs:
                singles.append(f)
    return pairs, singles


def _aligner() -> PairwiseAligner:
    al = PairwiseAligner(
        mode="global",
        match_score=MATCH,
        mismatch_score=MISMATCH,
        open_gap_score=GAP_OPEN,
        extend_gap_score=GAP_EXTEND,
    )
    al.end_gap_score = 0
    return al


def _gapped(seq_a: str, seq_b: str, aligned) -> tuple[str, str]:
    blocks_a = [tuple(int(x) for x in row) for row in aligned[0]]
    blocks_b = [tuple(int(x) for x in row) for row in aligned[1]]
    out_a: list[str] = []
    out_b: list[str] = []
    pa = pb = 0
    for (sa, ea), (sb, eb) in zip(blocks_a, blocks_b):
        while pa < sa:
            out_a.append(seq_a[pa])
            out_b.append("-")
            pa += 1
        while pb < sb:
            out_a.append("-")
            out_b.append(seq_b[pb])
            pb += 1
        if ea - sa != eb - sb:
            raise RuntimeError("정렬 블록 길이가 서로 다릅니다.")
        out_a.append(seq_a[sa:ea])
        out_b.append(seq_b[sb:eb])
        pa, pb = ea, eb
    while pa < len(seq_a):
        out_a.append(seq_a[pa])
        out_b.append("-")
        pa += 1
    while pb < len(seq_b):
        out_a.append("-")
        out_b.append(seq_b[pb])
        pb += 1
    return "".join(out_a), "".join(out_b)


def _expand(base: str) -> set[str]:
    return set(IUPAC_SETS.get(base, set("ACGT")))


def _code(bases: set[str]) -> str:
    return IUPAC_CODE.get(frozenset(bases), "N")


def overlap_consensus(
    seq_a: str,
    qual_a: list[int],
    seq_b: str,
    qual_b: list[int],
    *,
    conflict_q: float = 5,
    indel_q: float = 15,
) -> dict:
    """두 서열은 이미 같은 방향(역방향은 reverse-complement 된 상태)."""
    seq_a = seq_a.upper().replace("U", "T")
    seq_b = seq_b.upper().replace("U", "T")
    qa = list(qual_a) + [0] * max(0, len(seq_a) - len(list(qual_a)))
    qb = list(qual_b) + [0] * max(0, len(seq_b) - len(list(qual_b)))
    qa, qb = qa[: len(seq_a)], qb[: len(seq_b)]
    if not seq_a or not seq_b:
        return _empty_consensus()
    aln = _aligner().align(seq_a, seq_b)[0]
    ga, gb = _gapped(seq_a, seq_b, aln.aligned)
    cols = []
    ia = ib = 0
    for ca, cb in zip(ga, gb):
        if ca == "-":
            q1 = None
        else:
            q1 = qa[ia]
            ia += 1
        if cb == "-":
            q2 = None
        else:
            q2 = qb[ib]
            ib += 1
        cols.append((ca, cb, q1, q2))
    dual = [i for i, (a, b, _, _) in enumerate(cols) if a != "-" and b != "-"]
    if not dual:
        return _empty_consensus()
    lo, hi = dual[0], dual[-1]
    out = []
    oq = []
    overlap_len = 0
    compared = matches = conflicts = resolved = indels = 0
    for i, (a, b, q1, q2) in enumerate(cols):
        inside = lo <= i <= hi
        if a != "-" and b != "-":
            overlap_len += 1
            base, _was_conflict, was_resolved, _amb = _call_pair(a, q1, b, q2, conflict_q)
            if a in ACGT and b in ACGT:
                compared += 1
                if a == b:
                    matches += 1
                else:
                    conflicts += 1
                    if was_resolved:
                        resolved += 1
            out.append(base)
            oq.append(_qual_of(base, a, q1, b, q2, was_resolved))
            continue
        if not inside:
            if a != "-":
                out.append(a)
                oq.append(q1 or 0)
            elif b != "-":
                out.append(b)
                oq.append(q2 or 0)
            continue
        # 겹침 안의 indel
        base_c, qual = (a, q1) if a != "-" else (b, q2)
        indels += 1
        if base_c != "-" and (qual or 0) >= indel_q:
            out.append(base_c)
            oq.append(qual or 0)
    seq = "".join(out)
    amb = sum(1 for c in seq if c not in ACGT)
    ident = (matches / compared) if compared else 0.0
    return dict(
        seq=seq,
        quals=oq,
        overlap_len=overlap_len,
        overlap_ident=ident,
        n_conflicts=conflicts,
        n_resolved=resolved,
        n_ambiguities=amb,
        n_indels=indels,
        compared=compared,
    )


def _empty_consensus() -> dict:
    return dict(
        seq="", quals=[], overlap_len=0, overlap_ident=0.0,
        n_conflicts=0, n_resolved=0, n_ambiguities=0, n_indels=0, compared=0,
    )


def _call_pair(a, q1, b, q2, conflict_q: float):
    sa, sb = _expand(a), _expand(b)
    inter = sa & sb
    q1 = q1 or 0
    q2 = q2 or 0
    if inter:
        if sa == sb:
            return _code(sa), False, False, len(sa) != 1
        if sa < sb and q1 >= q2:
            return _code(sa), False, False, len(sa) != 1
        if sb < sa and q2 >= q1:
            return _code(sb), False, False, len(sb) != 1
        union = sa | sb
        return _code(union), False, False, len(union) != 1
    # 서로 배타적인 염기
    if q1 >= q2 + conflict_q:
        return a, True, True, a not in ACGT
    if q2 >= q1 + conflict_q:
        return b, True, True, b not in ACGT
    union = sa | sb
    return _code(union), True, False, True


def _qual_of(base, a, q1, b, q2, resolved: bool) -> int:
    q1 = int(q1 or 0)
    q2 = int(q2 or 0)
    if base in ACGT and a == b:
        return min(93, max(q1, q2))
    if resolved:
        return max(1, abs(q1 - q2))
    return min(q1, q2)


def _better(a: SangerRead, b: SangerRead) -> SangerRead:
    ka = (a.trim_mean_q, a.trim_len, a.trim_q20)
    kb = (b.trim_mean_q, b.trim_len, b.trim_q20)
    return a if ka >= kb else b


def _fmt_ident(frac: float | None) -> str:
    if frac is None:
        return "?"
    return f"{frac * 100:.2f}%"


def merge_pair(
    forward: SangerRead,
    reverse: SangerRead,
    *,
    min_overlap: int = 50,
    min_overlap_ident: float = 90.0,
    conflict_q: float = 5,
    indel_q: float = 15,
    backend: str = "biopython",
    tracy_trim: int = 4,
    sample_id: str | None = None,
    log=print,
) -> MergeResult:
    sid = sample_id or forward.sample_id
    if backend == "tracy":
        try:
            return _merge_tracy(
                forward, reverse,
                min_overlap=min_overlap,
                min_overlap_ident=min_overlap_ident,
                tracy_trim=tracy_trim,
                sample_id=sid,
                log=log,
            )
        except TracyUnavailable:
            raise
        except TracyError as e:
            log(f"  tracy 실패({e}). Biopython 합침으로 넘어갑니다.")
    rc = str(Seq(reverse.seq).reverse_complement())
    rq = list(reversed(reverse.quals))
    stats = overlap_consensus(
        forward.seq, forward.quals, rc, rq, conflict_q=conflict_q, indel_q=indel_q,
    )
    ident_pct = stats["overlap_ident"] * 100
    ok = (
        stats["overlap_len"] >= min_overlap
        and stats["compared"] > 0
        and ident_pct + 1e-9 >= min_overlap_ident
        and stats["seq"]
    )
    if ok:
        amb = stats["n_ambiguities"]
        warning = ""
        if amb:
            warning = (
                f"합침: 겹침 {stats['overlap_len']}bp, 일치도 {_fmt_ident(stats['overlap_ident'])}, "
                f"충돌 {stats['n_conflicts']} (품질로 해소 {stats['n_resolved']}), "
                f"미해결 모호염기 {amb}."
            )
        return MergeResult(
            seq=stats["seq"],
            quals=stats["quals"],
            mode="consensus",
            backend="biopython",
            overlap_len=stats["overlap_len"],
            overlap_ident=round(stats["overlap_ident"], 6),
            n_conflicts=stats["n_conflicts"],
            n_resolved=stats["n_resolved"],
            n_ambiguities=amb,
            n_indels=stats["n_indels"],
            warning=warning,
            direction="FR",
            trimmed_len=len(stats["seq"]),
            raw_len=forward.raw_len + reverse.raw_len,
            trim_mean_q=round(sum(stats["quals"]) / len(stats["quals"]), 1) if stats["quals"] else 0,
            sample_id=sid,
            reads=[forward, reverse],
        )
    if stats["overlap_len"] < min_overlap:
        why = f"겹침 {stats['overlap_len']}bp < {min_overlap}bp"
    elif stats["compared"] == 0:
        why = "겹침 구간에 A/C/G/T 비교 자리가 없음"
    else:
        why = f"겹침 일치도 {_fmt_ident(stats['overlap_ident'])} < {min_overlap_ident:g}%"
    chosen = _better(forward, reverse)
    which = "정방향" if chosen is forward else "역방향"
    warning = (
        f"[애매] 정·역방향 합침 실패: {why}. "
        f"더 좋은 단일 리드({which}, {chosen.trim_len}bp, 평균Q {chosen.trim_mean_q})를 사용."
    )
    return MergeResult(
        seq=chosen.seq,
        quals=list(chosen.quals),
        mode="fallback",
        backend="biopython",
        overlap_len=stats["overlap_len"],
        overlap_ident=round(stats["overlap_ident"], 6) if stats["compared"] else 0.0,
        n_conflicts=stats["n_conflicts"],
        n_resolved=stats["n_resolved"],
        n_ambiguities=stats["n_ambiguities"],
        n_indels=stats["n_indels"],
        warning=warning,
        direction="F" if chosen is forward else "R",
        trimmed_len=chosen.trim_len,
        raw_len=chosen.raw_len,
        trim_mean_q=chosen.trim_mean_q,
        sample_id=sid,
        reads=[forward, reverse],
    )


def single_result(read: SangerRead, warning: str = "") -> MergeResult:
    return MergeResult(
        seq=read.seq,
        quals=list(read.quals),
        mode="single",
        backend="none",
        warning=warning or read.warning,
        direction=read.direction,
        trimmed_len=read.trim_len,
        raw_len=read.raw_len,
        trim_mean_q=read.trim_mean_q,
        sample_id=read.sample_id,
        reads=[read],
    )


def _merge_tracy(forward, reverse, *, min_overlap, min_overlap_ident, tracy_trim, sample_id, log) -> MergeResult:
    exe = shutil.which("tracy")
    if not exe:
        raise TracyUnavailable(
            "tracy 가 PATH에 없습니다. "
            "https://github.com/gear-genomics/tracy/releases 의 정적 바이너리를 설치하거나 "
            "--merger biopython 을 쓰세요."
        )
    import tempfile

    with tempfile.TemporaryDirectory(prefix="blastid_tracy_") as tmp:
        prefix = str(Path(tmp) / "cons")
        cmd = [
            exe, "consensus",
            "-a",
            "-o", prefix,
            "-b", sample_id[:40] or "consensus",
            "-t", str(int(tracy_trim)),
            "-c", str(int(min_overlap)),
            "-f", str(float(min_overlap_ident) / 100.0),
            str(forward.path),
            str(reverse.path),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if proc.returncode != 0:
            msg = (proc.stderr or proc.stdout or "").strip().splitlines()
            raise TracyError(msg[-1] if msg else f"exit {proc.returncode}")
        fa = Path(prefix + ".fa")
        if not fa.exists():
            raise TracyError("cons.fa 가 없습니다.")
        seq = _first_fasta_seq(fa.read_text(encoding="utf-8", errors="replace"))
        seq = re.sub(r"[^A-Za-z]", "", seq).upper().replace("U", "T")
        align_fa = Path(prefix + ".align.fa")
        stats = _empty_consensus()
        if align_fa.exists():
            seqs = _fasta_seqs(align_fa.read_text(encoding="utf-8", errors="replace"))
            if len(seqs) >= 2:
                s1 = re.sub(r"[^A-Za-z\-]", "", seqs[0]).upper()
                s2 = re.sub(r"[^A-Za-z\-]", "", seqs[1]).upper()
                s1u, s2u = s1.replace("-", ""), s2.replace("-", "")
                if s1u and s2u:
                    if len(s1) == len(s2) and ("-" in s1 or "-" in s2):
                        stats = _column_stats(s1, s2)
                    else:
                        stats = overlap_consensus(s1u, [40] * len(s1u), s2u, [40] * len(s2u))
        ident_pct = stats["overlap_ident"] * 100
        ok = bool(seq) and stats["overlap_len"] >= min_overlap and ident_pct + 1e-9 >= min_overlap_ident
        if not ok:
            why = (
                f"tracy 합의 사용 안 함 (겹침 {stats['overlap_len']}bp, "
                f"일치도 {_fmt_ident(stats['overlap_ident'])}, 합의 {len(seq)}bp)"
            )
            log(f"  {why}")
            raise TracyError(why)
        amb = sum(1 for c in seq if c not in ACGT)
        warning = ""
        if amb:
            warning = (
                f"합침(tracy): 겹침 {stats['overlap_len']}bp, "
                f"일치도 {_fmt_ident(stats['overlap_ident'])}, 모호염기 {amb}."
            )
        return MergeResult(
            seq=seq,
            quals=[40] * len(seq),
            mode="consensus",
            backend="tracy",
            overlap_len=stats["overlap_len"],
            overlap_ident=round(stats["overlap_ident"], 6),
            n_conflicts=stats["n_conflicts"],
            n_resolved=stats["n_resolved"],
            n_ambiguities=amb,
            n_indels=stats["n_indels"],
            warning=warning,
            direction="FR",
            trimmed_len=len(seq),
            raw_len=forward.raw_len + reverse.raw_len,
            trim_mean_q=None,
            sample_id=sample_id,
            reads=[forward, reverse],
        )


def _column_stats(ga: str, gb: str) -> dict:
    """이미 맞춰 둔 같은 길이의 정렬(갭 포함)."""
    assert len(ga) == len(gb)
    dual = [i for i, (a, b) in enumerate(zip(ga, gb)) if a not in GAP and b not in GAP]
    if not dual:
        return _empty_consensus()
    compared = matches = conflicts = 0
    for i in dual:
        a, b = ga[i], gb[i]
        if a in ACGT and b in ACGT:
            compared += 1
            if a == b:
                matches += 1
            else:
                conflicts += 1
    return dict(
        seq="",
        quals=[],
        overlap_len=len(dual),
        overlap_ident=(matches / compared) if compared else 0.0,
        n_conflicts=conflicts,
        n_resolved=0,
        n_ambiguities=0,
        n_indels=sum(1 for a, b in zip(ga, gb) if (a in GAP) ^ (b in GAP)),
        compared=compared,
    )


def _fasta_seqs(text: str) -> list[str]:
    seqs = []
    cur: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def _first_fasta_seq(text: str) -> str:
    seqs = _fasta_seqs(text)
    return seqs[0] if seqs else ""


def assemble_files(
    paths: list[str | Path],
    *,
    pairs_csv: str | None = None,
    pair_regexes: list[str] | None = None,
    no_merge: bool = False,
    trimmer: str = "mott",
    trim_cutoff: float = 0.05,
    min_len: int = 150,
    min_q: float = 25,
    min_overlap: int = 50,
    min_overlap_ident: float = 90.0,
    conflict_q: float = 5,
    indel_q: float = 15,
    merger: str = "biopython",
    tracy_trim: int = 4,
    log=print,
) -> list[MergeResult]:
    files = collect_ab1(paths)
    explicit = read_pairs_csv(Path(pairs_csv)) if pairs_csv else None
    if explicit:
        for _sid, fp, rp in explicit:
            files.extend([fp, rp])
        seen: set[str] = set()
        deduped: list[Path] = []
        for f in files:
            key = str(Path(f).resolve())
            if key not in seen:
                seen.add(key)
                deduped.append(Path(f))
        files = deduped
    if not files:
        raise SystemExit("ab1/abi 파일이 없습니다.")
    reads: dict[str, SangerRead] = {}
    for f in files:
        try:
            r = read_ab1(f, trim_cutoff, trimmer)
        except Exception as e:
            log(f"읽기 실패 {f}: {e}")
            continue
        warns = []
        if r.no_quality:
            warns.append("품질값 없음(트리밍 안 함)")
        if r.trim_len < min_len:
            warns.append(f"트리밍 후 {r.trim_len}bp < {min_len}")
        if not r.no_quality and r.trim_mean_q < min_q:
            warns.append(f"평균 Q {r.trim_mean_q} < {min_q}")
        r.warning = "; ".join(warns)
        reads[str(f.resolve())] = r
        log(
            f"{f.name}: {r.raw_len}bp → {r.trim_len}bp ({r.trim_start}-{r.trim_end}), "
            f"평균Q {r.trim_mean_q}, Q20 {r.trim_q20}%, N {r.n_count}"
            + (f"  ⚠ {r.warning}" if r.warning else "")
        )
    if no_merge:
        out = []
        for r in reads.values():
            if r.trim_len < min_len or not r.seq:
                log(f"건너뜀 {r.stem}: 트리밍 후 길이가 {min_len}bp 미만")
                continue
            out.append(single_result(r))
        return out
    pairs, _singles = pair_reads(
        [r.path for r in reads.values()], pair_regexes, log=log, explicit=explicit,
    )
    # pair_reads에 넘긴 것은 Path 키. reads 키는 resolve().
    results: list[MergeResult] = []
    used = set()
    for sid, fp, rp in pairs:
        fk, rk = str(fp.resolve()), str(rp.resolve())
        fr, rr = reads.get(fk), reads.get(rk)
        if fr is None or rr is None:
            log(f"짝 {sid}: 읽기에 실패해 건너뜁니다.")
            continue
        used.add(fk)
        used.add(rk)
        fr.direction, rr.direction = "F", "R"
        fr.sample_id = rr.sample_id = sid
        if fr.trim_len < min_len and rr.trim_len < min_len:
            log(f"{sid}: 양쪽 다 {min_len}bp 미만이라 BLAST에서 제외")
            continue
        if fr.trim_len < min_len or rr.trim_len < min_len:
            keep = rr if fr.trim_len < min_len else fr
            short = "정방향" if fr.trim_len < min_len else "역방향"
            warn = (
                f"[애매] {short} 리드가 {min_len}bp 미만이라 합치지 않고 "
                f"{'역방향' if keep is rr else '정방향'}만 사용 ({keep.trim_len}bp)."
            )
            res = single_result(keep, warn)
            res.sample_id = sid
            res.mode = "fallback"
            res.direction = keep.direction
            results.append(res)
            log(f"{sid}: {warn}")
            continue
        try:
            res = merge_pair(
                fr, rr,
                min_overlap=min_overlap,
                min_overlap_ident=min_overlap_ident,
                conflict_q=conflict_q,
                indel_q=indel_q,
                backend=merger,
                tracy_trim=tracy_trim,
                sample_id=sid,
                log=log,
            )
        except TracyUnavailable as e:
            raise SystemExit(str(e)) from e
        if res.trimmed_len < min_len:
            log(f"{sid}: 합침 결과 {res.trimmed_len}bp < {min_len}, BLAST에서 제외")
            continue
        extra = f"  ⚠ {res.warning}" if res.warning else ""
        if res.mode == "consensus":
            log(
                f"{sid}: 합침({res.backend}) {res.trimmed_len}bp, "
                f"겹침 {res.overlap_len}bp, 일치도 {_fmt_ident(res.overlap_ident)}, "
                f"충돌 {res.n_conflicts}, 모호 {res.n_ambiguities}{extra}"
            )
        else:
            log(f"{sid}: {res.warning}")
        results.append(res)
    for f, r in reads.items():
        if f in used:
            continue
        # singles 목록은 짝이 안 된 파일. used에 없는 것 전부 단일 처리.
        if r.trim_len < min_len or not r.seq:
            log(f"건너뜀 {r.stem}: 트리밍 후 길이가 {min_len}bp 미만")
            continue
        hit = None
        if not no_merge:
            cre = _compile_patterns(pair_regexes or DEFAULT_PAIR_REGEXES)
            hit = match_pair_name(r.stem, cre)
        if hit and hit[1] in {"F", "R"}:
            r.direction = hit[1]
            r.sample_id = r.stem
        results.append(single_result(r))
        log(f"단일 리드: {r.stem} ({r.direction}, {r.trim_len}bp)")
    # 같은 BLAST 라벨이 생기지 않게 sample_id 중복을 나중에 pipeline에서 처리
    return results


def write_qc(outdir: Path, results: list[MergeResult]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fasta = outdir / "merged.fasta"
    fastq = outdir / "merged.fastq"
    with open(fasta, "w", encoding="utf-8") as fa, open(fastq, "w", encoding="utf-8") as fq:
        for res in results:
            if not res.seq:
                continue
            fa.write(
                f">{res.sample_id} mode={res.mode} backend={res.backend} "
                f"len={res.trimmed_len} overlap={res.overlap_len} "
                f"ident={_fmt_ident(res.overlap_ident)}\n{res.seq}\n"
            )
            quals = res.quals or [0] * len(res.seq)
            if len(quals) < len(res.seq):
                quals = quals + [0] * (len(res.seq) - len(quals))
            qstr = "".join(chr(33 + max(0, min(93, int(q)))) for q in quals[: len(res.seq)])
            fq.write(f"@{res.sample_id}\n{res.seq}\n+\n{qstr}\n")
    reads_rows = []
    for res in results:
        for r in res.reads:
            reads_rows.append(r)
    # 단일 리드가 여러 결과에 안 들어가면 위에서 포함됨
    keys = [
        "file", "id", "direction", "raw_len", "trim_start", "trim_end", "trim_len",
        "trim_mean_q", "trim_q20", "n_count", "warning",
    ]
    with open(outdir / "ab1_qc.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        seen = set()
        for r in reads_rows:
            key = str(r.path)
            if key in seen:
                continue
            seen.add(key)
            w.writerow(
                {
                    "file": str(r.path),
                    "id": r.sample_id,
                    "direction": r.direction,
                    "raw_len": r.raw_len,
                    "trim_start": r.trim_start,
                    "trim_end": r.trim_end,
                    "trim_len": r.trim_len,
                    "trim_mean_q": r.trim_mean_q,
                    "trim_q20": r.trim_q20,
                    "n_count": r.n_count,
                    "warning": r.warning,
                }
            )
    mkeys = [
        "id", "mode", "backend", "direction", "trimmed_len", "overlap_len",
        "overlap_ident", "n_conflicts", "n_resolved", "n_ambiguities", "n_indels", "warning",
    ]
    with open(outdir / "merge_qc.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=mkeys)
        w.writeheader()
        for res in results:
            w.writerow(
                {
                    "id": res.sample_id,
                    "mode": res.mode,
                    "backend": res.backend,
                    "direction": res.direction,
                    "trimmed_len": res.trimmed_len,
                    "overlap_len": res.overlap_len if res.overlap_len is not None else "",
                    "overlap_ident": res.overlap_ident if res.overlap_ident is not None else "",
                    "n_conflicts": res.n_conflicts if res.n_conflicts is not None else "",
                    "n_resolved": res.n_resolved if res.n_resolved is not None else "",
                    "n_ambiguities": res.n_ambiguities if res.n_ambiguities is not None else "",
                    "n_indels": res.n_indels if res.n_indels is not None else "",
                    "warning": res.warning,
                }
            )
