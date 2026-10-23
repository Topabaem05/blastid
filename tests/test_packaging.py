"""편집 설치가 아닌 wheel에도 작업창 파일이 들어 있는지 확인한다."""
import os
import subprocess
import sys
from pathlib import Path


def test_wheel_serves_taskpane_outside_the_repo(tmp_path: Path):
    repo = Path(__file__).resolve().parents[1]
    wheel_dir = tmp_path / "wheels"
    wheel_dir.mkdir()
    subprocess.check_call(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            "--no-deps",
            "--no-build-isolation",
            "-w",
            str(wheel_dir),
            str(repo),
        ],
    )
    wheels = list(wheel_dir.glob("blastid-*.whl"))
    assert wheels, "wheel이 만들어지지 않았습니다"
    target = tmp_path / "install"
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "--no-deps", "--target", str(target), str(wheels[0])],
    )
    away = tmp_path / "away"
    away.mkdir()
    script = tmp_path / "check_taskpane.py"
    script.write_text(
        """
import os, sys
sys.path.insert(0, os.environ["BLASTID_TARGET"])
os.chdir(os.environ["BLASTID_AWAY"])
from fastapi.testclient import TestClient
from blastid.server import create_app
client = TestClient(create_app())
page = client.get("/taskpane.html")
if page.status_code != 200:
    raise SystemExit(f"status {page.status_code}")
if "RID 최대 대기" not in page.text:
    raise SystemExit("taskpane.html 내용이 없습니다")
icon = client.get("/assets/icon-32.png")
if icon.status_code != 200 or not icon.content.startswith(b"\\x89PNG"):
    raise SystemExit(f"icon {icon.status_code}")
print("ok")
""",
        encoding="utf-8",
    )
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["BLASTID_TARGET"] = str(target)
    env["BLASTID_AWAY"] = str(away)
    subprocess.check_call([sys.executable, str(script)], cwd=away, env=env)
