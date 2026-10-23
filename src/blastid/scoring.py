"""웹 BLAST Descriptions 표와 같은 값 계산, 문제 행 판정.

공식은 NCBI C++ Toolkit ``align_format_util.cpp``
(GetScoreString, GetSeqAlignSetCalcParams, GetPercentIdentity)와
기존에 손으로 맞춘 프로토타입을 그대로 따른다.
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET

UNRES_RE = re.compile(
    r"\b(sp|aff|cf|nr)\.|uncultured|unidentified|environmental|fungal sp|lichen sp",
    re.I,
)


def fmt_evalue(e: float) -> str:
    if e < 1.0e-180:
        return "0.0"
    if e < 1.0e-99:
        return "%2.0e" % e
    if e < 0.0009:
        return "%3.0e" % e
    if e < 0.1:
        return "%4.3f" % e
    if e < 1.0:
        return "%3.2f" % e
    if e < 10.0:
        return "%2.1f" % e
    return "%2.0f" % e


def fmt_bits(b: float, total: bool = False) -> str:
    if b > 99999:
        return "%5.3e" % b
    if b > 99.9:
        return "%d" % int(b)  # (long)bit_score → 0 쪽으로 자름
    return ("%2.1f" if total else "%4.1f") % b


def union_len(intervals: list[tuple[int, int]]) -> int:
    """1-based inclusive 구간의 합집합 길이. 맞닿은 구간(끝+1)은 이어 붙인다."""
    iv = sorted(intervals)
    tot = 0
    cs, ce = iv[0]
    for a, b in iv[1:]:
        if a > ce + 1:
            tot += ce - cs + 1
            cs, ce = a, b
        else:
            ce = max(ce, b)
    return tot + ce - cs + 1


def _num(s: str):
    f = float(s)
    return int(f) if f.is_integer() else f


def acc_version(desc: dict) -> str:
    m = re.search(r"\|([A-Za-z0-9_]+\.\d+)\|?", desc.get("id", "") or "")
    return m.group(1) if m else (desc.get("accession") or "")


def web_values(hit: dict, qlen: int) -> dict:
    """JSON2_S hit 하나 → 웹 표에 보이는 문자열과 엑셀에 넣는 숫자."""
    hsps = hit["hsps"]
    best, hb = hsps[0], 0.0
    for h in hsps:  # 동점이면 앞 HSP (WB-1175: 최고 bit score HSP의 identity)
        if h["bit_score"] > hb:
            hb, best = h["bit_score"], h
    total = sum(h["bit_score"] for h in hsps)
    cov = union_len(
        [
            (min(h["query_from"], h["query_to"]), max(h["query_from"], h["query_to"]))
            for h in hsps
        ]
    )
    qc = int(0.5 + 100 * cov / float(qlen))
    if qc > 100:
        qc = min(99, qc)
    if best["identity"] == best["align_len"]:
        pid = 100.0
    else:
        pid = 100.0 * best["identity"] / best["align_len"]
    d = hit["description"][0]
    v = dict(
        name=d.get("sciname", "") or "",
        taxid=d.get("taxid"),
        title=d.get("title", "") or "",
        max_s=fmt_bits(hb),
        total_s=fmt_bits(total, True),
        qc_s=f"{qc}%",
        e_s=fmt_evalue(best["evalue"]),
        pi_s="%.2f%%" % pid,
        acc_len=int(hit["len"]),
        accession=acc_version(d),
    )
    v.update(_stored(v))
    return v


def _stored(v: dict) -> dict:
    return dict(
        max=_num(v["max_s"]),
        total=_num(v["total_s"]),
        qc=round(int(v["qc_s"][:-1]) / 100, 4),
        evalue=0 if v["e_s"] == "0.0" else float(v["e_s"]),
        pi=round(float(v["pi_s"][:-1]) / 100, 6),
    )


def genus(name: str | None) -> str:
    return (name or "").split(" ")[0]


def assess(hits: list[dict], thr_pi: float = 90.0, thr_qc: float = 80.0, alt_margin: float = 2.0):
    """hits: web_values 목록(웹 순서). 반환 (level, comment).

    level 은 None, 'ambiguous'([애매]), 'strange'([이상]) 중 하나.
    """
    if not hits:
        return "strange", "[이상] BLAST hit 없음(No significant similarity). 서열 확인 필요."
    top = hits[0]
    tpi = top["pi"] * 100
    tqc = top["qc"] * 100
    good = [h for h in hits if h["pi"] * 100 >= thr_pi and h["qc"] * 100 >= thr_qc]

    def s(h, i=None):
        r = f"{h['name']} {h['accession']} {h['pi_s']}/{h['qc_s']}"
        return (f"{i}위 " + r) if i else r

    rank = {id(h): i + 1 for i, h in enumerate(hits)}
    if tpi < thr_pi or tqc < thr_qc:
        why = []
        if tpi < thr_pi:
            why.append(f"Per. Ident {thr_pi:.0f}% 미만")
        if tqc < thr_qc:
            why.append(f"Query Cover {thr_qc:.0f}% 미만")
        alt, nm = [], set()
        for h in good:
            if h["name"] not in nm:
                nm.add(h["name"])
                alt.append(h)
        alt = alt[:4]
        txt = (
            f"[이상] 1위 {top['name']} {top['pi_s']}/{top['qc_s']} ({', '.join(why)}). "
            + (
                "90%↑·Cover 80%↑ 후보: " + "; ".join(s(h, rank[id(h)]) for h in alt) + ". "
                if alt
                else "90%↑·Cover 80%↑ hit 없음. "
            )
            + "서열 품질·혼입 확인 또는 재시퀀싱 권장."
        )
        hi = max(hits, key=lambda h: (h["pi"], h["max"]))
        if hi is not top and hi["pi"] > top["pi"]:
            txt += f" 최고 Ident hit: {s(hi, rank[id(hi)])}."
        return "strange", txt
    tg = genus(top["name"])
    other_genus, seen = [], set()
    for h in hits[1:]:
        if (
            h in good
            and genus(h["name"]) != tg
            and h["name"] not in seen
            and not UNRES_RE.search(h["name"])
            and h["pi"] * 100 >= tpi - alt_margin
        ):
            seen.add(h["name"])
            other_genus.append(h)
    ties = []
    for h in hits[1:]:
        if (
            h["max"] == top["max"]
            and h["pi"] == top["pi"]
            and h["name"] != top["name"]
            and h["name"] not in [t["name"] for t in ties]
        ):
            ties.append(h)
    unresolved = bool(UNRES_RE.search(top["name"]))
    named = []
    if unresolved:
        for h in good[1:]:
            if not UNRES_RE.search(h["name"]) and h["name"] not in [x["name"] for x in named]:
                named.append(h)
    if not other_genus and not ties and not unresolved:
        return None, ""
    parts = [f"[애매] 1위 {top['name']} {top['pi_s']}/{top['qc_s']}."]
    if unresolved:
        parts.append(
            "1위가 종 미확정(sp./aff./cf./uncultured) 서열 → 속까지만 신뢰."
            + (
                " 이름 있는 종 후보(90%↑·Cover 80%↑): "
                + "; ".join(s(h, rank[id(h)]) for h in named[:4])
                + "."
                if named
                else " 90%↑·Cover 80%↑ 이름 있는 종 후보 없음."
            )
        )
    if ties:
        parts.append(
            "1위와 점수·Ident 동점인 다른 종: " + "; ".join(s(h, rank[id(h)]) for h in ties[:5]) + "."
        )
    if other_genus:
        parts.append(
            f"다른 속 후보(90%↑·Cover 80%↑, 1위와 Ident 차 {alt_margin:g}%p 이내): "
            + "; ".join(s(h, rank[id(h)]) for h in other_genus[:5])
            + "."
        )
    hi = max(hits, key=lambda h: (h["pi"], h["max"]))
    if hi["pi"] * 100 >= tpi + 1:
        parts.append(
            f"참고: 더 높은 Ident hit {s(hi, rank[id(hi)])}"
            + (" (Cover 낮음)." if hi["qc"] * 100 < thr_qc else ".")
        )
    parts.append("종(또는 속) 판정 주의.")
    return "ambiguous", " ".join(parts)


def alternatives_text(hits: list[dict], n: int = 5) -> str:
    if not hits:
        return ""
    parts = []
    seen = {hits[0]["name"]}
    for i, h in enumerate(hits[1:], start=2):
        if h["name"] in seen:
            continue
        seen.add(h["name"])
        parts.append(f"{i}위 {h['name']} {h['accession']} {h['pi_s']}/{h['qc_s']}")
        if len(parts) >= n:
            break
    return "; ".join(parts)


def fungi_lineage(full: str | None) -> str:
    full = full or ""
    i = full.find("Fungi")
    return full[i:] if i >= 0 else full


def parse_xml_report(xml_text: str) -> list[dict | None]:
    """BLAST XML(다른 포맷, 다른 파서)로 쿼리별 1위 hit 값을 다시 계산."""
    root = ET.fromstring(xml_text)
    out: list[dict | None] = []
    for it in root.iter("Iteration"):
        qlen = int(it.findtext("Iteration_query-len"))
        hits = it.find("Iteration_hits")
        h = hits.find("Hit") if hits is not None else None
        if h is None:
            out.append(None)
            continue
        hs = h.find("Hit_hsps").findall("Hsp")
        bits = [float(x.findtext("Hsp_bit-score")) for x in hs]
        bi = max(range(len(hs)), key=lambda i: (bits[i], -i))
        b = hs[bi]
        qf = [
            (
                min(int(x.findtext("Hsp_query-from")), int(x.findtext("Hsp_query-to"))),
                max(int(x.findtext("Hsp_query-from")), int(x.findtext("Hsp_query-to"))),
            )
            for x in hs
        ]
        ident, alen = int(b.findtext("Hsp_identity")), int(b.findtext("Hsp_align-len"))
        hid = h.findtext("Hit_id") or ""
        src = hid.split("|")[3] if hid.count("|") >= 3 else hid
        m = re.search(r"([A-Za-z0-9_]+\.\d+)", src)
        raw_qc = int(0.5 + 100 * union_len(qf) / qlen)
        if raw_qc > 100:
            raw_qc = min(99, raw_qc)
        ev = fmt_evalue(float(b.findtext("Hsp_evalue")))
        out.append(
            dict(
                qlen=qlen,
                accession=m.group(1) if m else h.findtext("Hit_accession"),
                max=_num(fmt_bits(max(bits))),
                total=_num(fmt_bits(sum(bits), True)),
                qc=round(raw_qc / 100, 4),
                evalue=0 if ev == "0.0" else float(ev),
                pi=round(float("%.2f" % (100.0 * ident / alen)) / 100, 6),
                acc_len=int(h.findtext("Hit_len")),
            )
        )
    return out


def flag_label(level: str | None) -> str:
    if level == "strange":
        return "[이상]"
    if level == "ambiguous":
        return "[애매]"
    return ""
