"""웹 표 값 계산. 기대값은 NCBI 표시 규칙을 손으로 계산한 것이다."""
import math

from blastid.scoring import (
    assess,
    fmt_bits,
    fmt_evalue,
    fungi_lineage,
    parse_xml_report,
    union_len,
    web_values,
)
from tests.blastlib import blast_xml, json_hit


def test_fmt_evalue_golden():
    assert fmt_evalue(0) == "0.0"
    assert fmt_evalue(1e-181) == "0.0"
    assert fmt_evalue(1.2e-100) == "1e-100"
    assert fmt_evalue(1.23e-50) == "1e-50"
    assert fmt_evalue(5.5e-20) == "5e-20"
    assert fmt_evalue(0.00089) == "9e-04"
    assert fmt_evalue(0.0009) == "0.001"
    assert fmt_evalue(0.01234) == "0.012"
    assert fmt_evalue(0.0994) == "0.099"
    assert fmt_evalue(0.5) == "0.50"
    assert fmt_evalue(1.26) == "1.3"
    assert fmt_evalue(9.94) == "9.9"
    assert fmt_evalue(10) == "10"
    assert fmt_evalue(42) == "42"


def test_fmt_bits_truncates_above_99_9():
    assert fmt_bits(0) == " 0.0"
    assert fmt_bits(0, total=True) == "0.0"
    assert fmt_bits(9.96) == "10.0"
    assert fmt_bits(99.9) == "99.9"
    assert fmt_bits(99.94) == "99"  # 99.94 > 99.9 이면 내림
    assert fmt_bits(100.9) == "100"
    assert fmt_bits(900.5) == "900"
    assert fmt_bits(500.4) == "500"
    assert fmt_bits(100000) == "1.000e+05"


def test_union_len_merges_touching_and_overlap():
    assert union_len([(1, 10), (11, 20)]) == 20
    assert union_len([(1, 10), (12, 20)]) == 19
    assert union_len([(50, 120), (1, 100)]) == 120
    assert union_len([(1, 60), (100, 300)]) == 261


def test_multi_hsp_matches_hand_calculation_and_xml():
    # 더 낮은 HSP를 앞에 두어 "첫 HSP"가 아니라 최고 bit score HSP를 쓰는지 확인.
    spec = {
        "name": "Lecidella biformis",
        "taxid": 999,
        "acc": "PZ322662.1",
        "title": "Lecidella biformis ITS",
        "length": 540,
        "hsps": [
            {"bit": 100.2, "ident": 50, "alen": 60, "qfrom": 1, "qto": 60, "evalue": 1e-5},
            {"bit": 400.2, "ident": 200, "alen": 210, "qfrom": 100, "qto": 300, "evalue": 1.23e-50},
        ],
    }
    qlen = 320
    v = web_values(json_hit(spec), qlen)
    assert v["name"] == "Lecidella biformis"
    assert v["accession"] == "PZ322662.1"
    assert v["max"] == 400
    assert v["total"] == 500
    assert v["qc"] == 0.82
    assert v["qc_s"] == "82%"
    assert v["pi_s"] == "95.24%"
    assert v["pi"] == 0.9524
    assert v["e_s"] == "1e-50"
    assert math.isclose(v["evalue"], 1e-50)
    assert v["acc_len"] == 540
    # 0.0 은 정수 0 으로 저장
    zero = web_values(
        json_hit(
            {
                **spec,
                "hsps": [{"bit": 10, "ident": 10, "alen": 10, "qfrom": 1, "qto": 10, "evalue": 0.0}],
            }
        ),
        10,
    )
    assert zero["evalue"] == 0 and isinstance(zero["evalue"], int)

    seq = "A" * qlen
    xml = blast_xml([("q", seq, [spec])])
    parsed = parse_xml_report(xml)[0]
    for key in ("max", "total", "qc", "pi", "evalue", "acc_len", "accession"):
        assert parsed[key] == v[key], key


def test_query_cover_above_100_is_clamped():
    spec = {
        "name": "X",
        "taxid": 1,
        "acc": "A.1",
        "length": 10,
        "hsps": [{"bit": 20, "ident": 10, "alen": 12, "qfrom": 1, "qto": 12, "evalue": 0.1}],
    }
    v = web_values(json_hit(spec), 10)
    assert v["qc_s"] == "99%"
    assert v["qc"] == 0.99


def test_fungi_lineage_trims_prefix():
    full = "cellular organisms; Eukaryota; Opisthokonta; Fungi; Dikarya; Ascomycota"
    assert fungi_lineage(full) == "Fungi; Dikarya; Ascomycota"
    assert fungi_lineage("Bacteria; something") == "Bacteria; something"
    assert fungi_lineage("") == ""


def _hit(name, pi=0.99, qc=1.0, maxs=200, acc="A.1"):
    return dict(
        name=name,
        pi=pi,
        qc=qc,
        max=maxs,
        accession=acc,
        pi_s=f"{pi * 100:.2f}%",
        qc_s=f"{round(qc * 100):.0f}%",
    )


def test_assess_flags():
    clean = [_hit("Pyxine copelandii")]
    assert assess(clean) == (None, "")

    low = [_hit("Lecidella biformis", pi=0.8721, qc=0.96)]
    level, text = assess(low)
    assert level == "strange" and text.startswith("[이상]")
    assert "90%" in text and "87.21%" in text

    cover = [_hit("Collema furfuraceum", pi=0.95, qc=0.70)]
    level, text = assess(cover)
    assert level == "strange" and "Query Cover" in text

    sp = [_hit("Trimmatothelopsis sp.", pi=0.9486)]
    level, text = assess(sp)
    assert level == "ambiguous" and "미확정" in text

    tie = [_hit("Alpha beta", maxs=100, pi=0.99), _hit("Gamma delta", maxs=100, pi=0.99)]
    level, text = assess(tie)
    assert level == "ambiguous" and "동점" in text

    near = [_hit("Alpha beta", pi=0.99), _hit("Other genus", pi=0.98)]
    level, text = assess(near, alt_margin=2)
    assert level == "ambiguous" and "다른 속" in text

    far = [_hit("Alpha beta", pi=0.99), _hit("Other genus", pi=0.90)]
    assert assess(far, alt_margin=2)[0] is None

    same_genus = [_hit("Alpha beta", pi=0.99), _hit("Alpha gamma", pi=0.98)]
    assert assess(same_genus)[0] is None

    assert assess([])[0] == "strange"
