"""로컬 서버는 작업창이 넘긴 행만 받고, NCBI 고지를 보여 준다."""
import time

from fastapi.testclient import TestClient

from blastid.ncbi import TOOL_NAME, blast_identity_params
from blastid.server import create_app, ensure_dev_cert


def test_dev_cert_is_created(tmp_path):
    cert, key = ensure_dev_cert(tmp_path / "certs")
    assert cert.exists() and key.exists()
    assert b"BEGIN CERTIFICATE" in cert.read_bytes()


def test_health_and_disclaimer():
    client = TestClient(create_app())
    health = client.get("/api/health")
    assert health.status_code == 200
    body = health.json()
    assert body["ok"] is True
    assert body["tool"] == "blastid"
    assert "KST" in body["window"]
    notice = client.get("/api/disclaimer").json()
    assert "U.S. Government" in notice["en"]
    assert "보증하지 않" in notice["ko"]
    assert "tool=blastid" in blast_identity_params("") or TOOL_NAME == "blastid"
    params = blast_identity_params("")
    assert params[0] == ("tool", "blastid")
    assert all(key != "email" for key, _ in params)


def test_job_dry_run_lists_rows_without_ncbi(monkeypatch):
    def boom(*_a, **_k):
        raise AssertionError("dry-run 인데 NCBI 객체가 만들어지면 안 됩니다")

    monkeypatch.setattr("blastid.pipeline.NCBI", boom)
    client = TestClient(create_app())
    seq = "ACGT" * 20
    started = client.post(
        "/api/jobs",
        json={"rows": [{"row": 2, "id": "Synth1", "seq": seq, "bold": False, "name": ""}], "dry_run": True},
    )
    assert started.status_code == 200
    job_id = started.json()["id"]
    body = None
    for _ in range(40):
        body = client.get(f"/api/jobs/{job_id}").json()
        if body["status"] in {"done", "error"}:
            break
        time.sleep(0.05)
    assert body["status"] == "done", body
    text = "\n".join(event["message"] for event in body["events"])
    assert "Synth1" in text
    assert body["result"]["dry_run"] is True
    assert body["result"]["updates"] == []
    manifest = client.get("/manifest.xml")
    assert manifest.status_code == 200
    assert b"TaskPaneApp" in manifest.content
    assert b"taskpane.html" in manifest.content
