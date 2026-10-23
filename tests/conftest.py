import pytest


@pytest.fixture(autouse=True)
def _isolate_ncbi_usage(tmp_path, monkeypatch):
    """테스트가 홈 디렉터리의 24시간 검색 카운터를 건드리지 않게 한다."""
    monkeypatch.setenv("BLASTID_USAGE_FILE", str(tmp_path / "ncbi_usage.json"))
