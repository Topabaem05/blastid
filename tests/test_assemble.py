"""정·역방향 합침, Mott 트리밍, 파일명 짝."""
import shutil
from pathlib import Path

import pytest
from Bio.Seq import Seq

from blastid.assemble import (
    DEFAULT_PAIR_REGEXES,
    MergeResult,
    SangerRead,
    assemble_files,
    classify_direction,
    match_pair_name,
    merge_pair,
    mott_trim,
    overlap_consensus,
    pair_reads,
    read_ab1,
    _compile_patterns,
)
from tests.abifutil import write_ab1


def _read(seq, quals=None, qmean=40.0, direction="F", name="S") -> SangerRead:
    quals = list(quals) if quals is not None else [40] * len(seq)
    return SangerRead(
        path=Path(f"{name}_{direction}.ab1"),
        sample_id=name,
        stem=f"{name}_{direction}",
        seq=seq,
        quals=quals,
        raw_len=len(seq),
        trim_start=1,
        trim_end=len(seq),
        trim_mean_q=qmean,
        trim_q20=100.0,
        n_count=seq.count("N"),
        no_quality=False,
        direction=direction,
    )


def test_mott_trim_drops_low_quality_ends():
    quals = [2] * 10 + [40] * 120 + [2] * 10
    bs, be = mott_trim(quals, 0.05)
    assert (bs, be) == (10, 130)


def test_consensus_includes_flanks_and_reverse_complements():
    core = "ACGT" * 30
    fwd = "AAAA" + core
    oriented = core + "TTTT"
    rev = str(Seq(oriented).reverse_complement())
    res = merge_pair(_read(fwd, direction="F"), _read(rev, quals=[40] * len(rev), direction="R"))
    assert res.mode == "consensus"
    assert res.backend == "biopython"
    assert res.seq == fwd + "TTTT"
    assert res.overlap_len == len(core)
    assert res.overlap_ident == 1
    assert res.n_conflicts == 0
    assert res.direction == "FR"


def test_quality_resolves_conflict_or_writes_iupac():
    fwd = "ACGTACGTAC"
    oriented = "ACGTGCGTAC"  # index 4 A/G
    low = [40, 40, 40, 40, 10, 40, 40, 40, 40, 40]
    rev_low = str(Seq(oriented).reverse_complement())
    res = merge_pair(
        _read(fwd), _read(rev_low, quals=list(reversed(low)), direction="R"),
        min_overlap=5, min_overlap_ident=80,
    )
    assert res.mode == "consensus"
    assert res.seq == fwd
    assert res.n_conflicts == 1
    assert res.n_resolved == 1
    assert res.n_ambiguities == 0

    both = [40] * 10
    res2 = merge_pair(
        _read(fwd), _read(rev_low, quals=list(reversed(both)), direction="R"),
        min_overlap=5, min_overlap_ident=80,
    )
    assert res2.mode == "consensus"
    assert res2.seq == "ACGTRCGTAC"
    assert res2.n_conflicts == 1
    assert res2.n_resolved == 0
    assert res2.n_ambiguities == 1


def test_internal_indel_kept_only_when_quality_is_high():
    a = "AAAACGTAAA"
    b = "AAAAGTAAA"
    kept = overlap_consensus(a, [40] * len(a), b, [40] * len(b), indel_q=15)
    assert kept["seq"] == "AAAACGTAAA"
    assert kept["n_indels"] >= 1
    assert kept["overlap_ident"] == 1
    dropped = overlap_consensus(a, [5] * len(a), b, [40] * len(b), indel_q=15)
    assert dropped["seq"] == "AAAAGTAAA"


def test_short_overlap_falls_back_to_better_read():
    fwd = "AAA" + "ACGT" * 5 + "CCC"  # 겹침 20bp
    oriented = "ACGT" * 5 + "GGGG"
    rev = str(Seq(oriented).reverse_complement())
    res = merge_pair(
        _read(fwd, qmean=35),
        _read(rev, quals=[40] * len(rev), qmean=20, direction="R"),
        min_overlap=50,
    )
    assert res.mode == "fallback"
    assert res.seq == fwd
    assert res.direction == "F"
    assert res.warning.startswith("[애매]")
    assert "겹침" in res.warning


def test_low_identity_falls_back():
    res = merge_pair(_read("A" * 80, qmean=30), _read("C" * 80, qmean=10, direction="R"), min_overlap=20)
    assert res.mode == "fallback"
    assert res.seq == "A" * 80
    assert res.overlap_len == 0


@pytest.mark.parametrize(
    "stem,ident,direction",
    [
        ("OhF_F", "OhF", "F"),
        ("OhF_R", "OhF", "R"),
        ("OhF-F", "OhF", "F"),
        ("OhF-fwd", "OhF", "F"),
        ("OhF_rev", "OhF", "R"),
        ("OhF_forward", "OhF", "F"),
        ("OhF_reverse", "OhF", "R"),
        ("OhF_ITS1F", "OhF", "F"),
        ("OhF_ITS4", "OhF", "R"),
        ("OhF_ITS4R", "OhF", "R"),
        ("OhF.F", "OhF", "F"),
        ("Sample_01_R", "Sample_01", "R"),
    ],
)
def test_builtin_pair_names(stem, ident, direction):
    hit = match_pair_name(stem, _compile_patterns(DEFAULT_PAIR_REGEXES))
    assert hit == (ident, direction)
    assert classify_direction(hit[1] and {"F": "f", "R": "r"}[direction]) == direction


def test_pair_reads_and_csv(tmp_path: Path):
    files = [tmp_path / "S1_F.ab1", tmp_path / "S1_R.ab1", tmp_path / "lonely.ab1"]
    for f in files:
        f.write_bytes(b"x")
    pairs, singles = pair_reads(files, log=lambda *_: None)
    assert len(pairs) == 1 and pairs[0][0] == "S1"
    assert [p.name for p in singles] == ["lonely.ab1"]

    odd_f = tmp_path / "left.ab1"
    odd_r = tmp_path / "right.ab1"
    odd_f.write_bytes(b"x")
    odd_r.write_bytes(b"x")
    csvp = tmp_path / "pairs.csv"
    csvp.write_text("id,forward,reverse\nCustom,left.ab1,right.ab1\n", encoding="utf-8")
    pairs, singles = pair_reads([odd_f, odd_r], pairs_csv=str(csvp), log=lambda *_: None)
    assert pairs[0][0] == "Custom"
    assert singles == []


def test_read_ab1_and_assemble_roundtrip(tmp_path: Path):
    core = "ACGT" * 40
    fwd = "AAAA" + core
    oriented = core + "GGGG"
    rev = str(Seq(oriented).reverse_complement())
    write_ab1(tmp_path / "Spec_F.ab1", fwd, [40] * len(fwd), "Spec_F")
    write_ab1(tmp_path / "Spec_R.ab1", rev, [40] * len(rev), "Spec_R")
    got = read_ab1(tmp_path / "Spec_F.ab1")
    assert got.seq == fwd
    assert got.trim_len == len(fwd)
    results = assemble_files([tmp_path], min_len=50, log=lambda *_: None)
    assert len(results) == 1
    res = results[0]
    assert isinstance(res, MergeResult)
    assert res.sample_id == "Spec"
    assert res.mode == "consensus"
    assert res.seq == fwd + "GGGG"
    assert res.overlap_ident == 1


def test_unpaired_read_stays_single(tmp_path: Path):
    seq = "ACGT" * 40
    write_ab1(tmp_path / "only.ab1", seq, [40] * len(seq))
    results = assemble_files([tmp_path], min_len=50, log=lambda *_: None)
    assert len(results) == 1
    assert results[0].mode == "single"
    assert results[0].seq == seq


def test_custom_pair_regex(tmp_path: Path):
    seq = "ACGT" * 40
    write_ab1(tmp_path / "AaF.ab1", seq, [40] * len(seq))
    write_ab1(tmp_path / "AaR.ab1", str(Seq(seq).reverse_complement()), [40] * len(seq))
    results = assemble_files(
        [tmp_path],
        pair_regexes=[r"^(?P<id>.+?)(?P<dir>[FR])$"],
        min_len=50,
        log=lambda *_: None,
    )
    assert len(results) == 1 and results[0].mode == "consensus"


@pytest.mark.skipif(shutil.which("tracy") is None, reason="tracy가 PATH에 없음")
def test_tracy_backend_optional(tmp_path: Path):
    core = "ACGT" * 40
    fwd = "AAAA" + core
    rev = str(Seq(core + "GGGG").reverse_complement())
    write_ab1(tmp_path / "T_F.ab1", fwd, [40] * len(fwd))
    write_ab1(tmp_path / "T_R.ab1", rev, [40] * len(rev))
    results = assemble_files(
        [tmp_path], min_len=50, merger="tracy", tracy_trim=1, min_overlap=20, log=lambda *_: None
    )
    assert len(results) == 1
    assert results[0].seq
    assert results[0].backend in {"tracy", "biopython"}
    assert set(results[0].seq) <= set("ACGTNRYSWKMBDHV")
