"""CSV 입력은 엑셀과 같은 머리글로 읽고, UTF-8 BOM과 CP949를 받는다."""
from pathlib import Path

from blastid.pipeline import detect_kind, execute, read_text_auto, workbook_from_csv


def _run(path: Path, tmp_path: Path):
    from tests.test_pipeline import _args

    execute(
        _args(
            kind="excel",
            inputs=[str(path)],
            out_dir=str(tmp_path / "run"),
            dry_run=True,
        )
    )
    logs = list((tmp_path / "run").rglob("run.log"))
    assert logs
    return logs[0].read_text(encoding="utf-8")


def test_utf8_sig_csv_dry_run(tmp_path: Path):
    seq = "ACGT" * 20
    text = f"ID(Voucher),Seq\nS1,{seq}\n"
    path = tmp_path / "samples.csv"
    path.write_bytes(b"\xef\xbb\xbf" + text.encode("utf-8"))
    assert read_text_auto(path).startswith("ID(Voucher)")
    log = _run(path, tmp_path)
    assert "S1" in log
    assert "dry-run" in log
    assert detect_kind([str(path)]) == "excel"


def test_cp949_csv_headers(tmp_path: Path):
    seq = "ACGT" * 20
    text = f"ID,서열\nS2,{seq}\n"
    path = tmp_path / "kr.csv"
    path.write_bytes(text.encode("cp949"))
    wb = workbook_from_csv(path)
    ws = wb.active
    assert ws.cell(1, 1).value == "ID"
    assert ws.cell(1, 2).value == "서열"
    assert ws.cell(2, 1).value == "S2"
    log = _run(path, tmp_path)
    assert "S2" in log
