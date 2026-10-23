from pathlib import Path

import pytest
import yaml

from blastid.template import build_workbook, init_template, load_template

EXPECTED_HEADERS = [
    "No",
    "Collection No.",
    "ID(Voucher)",
    "Lichen name(1차 가동정)",
    "Substratum",
    "과제",
    "ELF",
    "납품",
    "DNA",
    "Seq",
    "Scientific Name",
    "Max Score",
    "Total Score",
    "Query Cover",
    "E value",
    "Per. Ident",
    "Acc. Len",
    "Accession",
    "Lineage",
    "검색 출처",
    "검토 코멘트 / 스크린샷",
    "스크린샷 (링크)",
]


def test_default_template_matches_workbook_layout():
    tmpl = load_template()
    assert [c.header for c in tmpl.columns] == EXPECTED_HEADERS
    assert tmpl.freeze == "D2"
    assert tmpl.row_height == 30
    by_field = {c.field: c for c in tmpl.columns if c.field}
    assert by_field["query_cover"].number_format == "0%"
    assert by_field["per_ident"].number_format == "0.00%"
    assert "검토 코멘트" in by_field["comment"].aliases
    assert "스크린샷" in by_field["screenshot"].aliases
    wb = build_workbook(tmpl)
    ws = wb[tmpl.sheet]
    headers = [ws.cell(1, i).value for i in range(1, len(EXPECTED_HEADERS) + 1)]
    assert headers == EXPECTED_HEADERS
    assert ws.freeze_panes == "D2"
    assert ws.cell(1, 11).fill.fgColor.rgb.endswith("1F4E78")


def test_init_template_writes_yaml_and_xlsx(tmp_path: Path):
    yml, xlsx = init_template(tmp_path / "tmpl")
    assert yml.exists() and xlsx.exists()
    data = yaml.safe_load(yml.read_text(encoding="utf-8"))
    assert data["workbook"] == "template.xlsx"
    tmpl = load_template(yml)
    assert tmpl.workbook_path == xlsx.resolve()
    assert [c.header for c in tmpl.columns] == EXPECTED_HEADERS


def test_unknown_field_is_rejected(tmp_path: Path):
    raw = {
        "columns": [{"header": "X", "field": "not_a_field"}],
    }
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(SystemExit):
        load_template(path)
