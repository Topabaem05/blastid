"""NCBI BLAST / E-utilities 사용 지침.

제출 간격과 RID 폴링은 ``ncbi.py`` 가 지킨다.
여기서는 50건 초과 시 한가한 시간, 24시간 검색 수, 고지문을 다룬다.
검색 1건은 RID 한 번(묶음 제출 1회)이다. NCBI도 묶음을 검색 1회로 센다.
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

ET = ZoneInfo("America/New_York")
KST = ZoneInfo("Asia/Seoul")

# 평일 미국 동부 21:00 이상, 다음날 05:00 미만. 주말은 하루 종일.
NIGHT_START_HOUR = 21
NIGHT_END_HOUR = 5
MANY_SEARCHES = 50
DAILY_SOFT_LIMIT = 100
APPROACH_AT = 80
WINDOW_SECONDS = 24 * 60 * 60

DISCLAIMER_EN = (
    "BLAST is a public-domain service of the National Center for Biotechnology Information (NCBI). "
    "Liability: For documents and software available from NCBI servers, the U.S. Government does not "
    "warrant or assume any legal liability or responsibility for the accuracy, completeness, or usefulness "
    "of any information, apparatus, product, or process disclosed. "
    "Endorsement: NCBI does not endorse or recommend any commercial products, processes, or services. "
    "The views and opinions of authors expressed on NCBI websites do not necessarily state or reflect "
    "those of the U.S. Government, and they may not be used for advertising or product endorsement purposes. "
    "https://www.ncbi.nlm.nih.gov/home/about/policies/ "
    "https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html"
)

DISCLAIMER_KO = (
    "이 기능은 NCBI BLAST와 E-utilities로 동정 결과를 가져옵니다. "
    "미국 정부는 그 정보의 정확성·완전성·유용성을 보증하지 않으며 법적 책임을 지지 않습니다. "
    "NCBI는 특정 상업 제품·방법·서비스를 보증하거나 권하지 않습니다. "
    "저자의 견해가 미국 정부의 입장을 뜻하지는 않으며, 광고나 제품 보증에 쓸 수 없습니다. "
    "BLAST 웹 서비스는 공동 자원입니다. 제출 간격은 10초 이상, 같은 RID 확인은 1분에 한 번 이하, "
    "E-utilities는 초당 3회 이하로 호출합니다. "
    "한 번에 50건(RID 기준)을 넘기면 주말이나 평일 미국 동부 21:00–05:00에 실행하는 것이 NCBI 지침입니다. "
    "24시간 원격 검색이 100건을 넘기게 되면 나머지는 이 컴퓨터의 BLAST+로 실행합니다. "
    "IP를 바꿔 한도를 피하지 않습니다. "
    "그 이상으로 NCBI에 보내면 느린 대기열로 옮기거나 요청을 막을 수 있습니다."
)


def _as_et(moment: datetime) -> datetime:
    if moment.tzinfo is None:
        raise ValueError("시각에는 시간대가 있어야 합니다.")
    return moment.astimezone(ET)


def is_off_hours(moment: datetime) -> bool:
    """주말 종일, 또는 평일 21:00–05:00 미국 동부(서머타임 포함)."""
    et = _as_et(moment)
    if et.weekday() >= 5:
        return True
    return et.hour >= NIGHT_START_HOUR or et.hour < NIGHT_END_HOUR


def next_off_hours_start(moment: datetime) -> datetime:
    """이미 한가한 시간이면 그 시각, 아니면 다음 21:00 ET."""
    et = _as_et(moment)
    if is_off_hours(et):
        return et
    start = et.replace(hour=NIGHT_START_HOUR, minute=0, second=0, microsecond=0)
    if et >= start:
        start = start + timedelta(days=1)
    return start


def _binary_edge(left: datetime, right: datetime, *, left_is_off: bool) -> datetime:
    """left와 right 사이에서 한가한 시간이 바뀌는 시각. 결과는 분 단위."""
    lo, hi = left, right
    while (hi - lo) > timedelta(minutes=1):
        mid = lo + (hi - lo) / 2
        if is_off_hours(mid) == left_is_off:
            lo = mid
        else:
            hi = mid
    edge = hi.astimezone(ET).replace(second=0, microsecond=0)
    # 05:00·21:00 경계로 붙인다.
    candidates = []
    base = edge.replace(minute=0)
    for day in (-1, 0, 1):
        day0 = (base + timedelta(days=day)).replace(hour=0)
        candidates.append(day0.replace(hour=NIGHT_END_HOUR))
        candidates.append(day0.replace(hour=NIGHT_START_HOUR))
    return min(candidates, key=lambda item: abs((item - edge).total_seconds()))


def off_hours_span(moment: datetime) -> tuple[datetime, datetime]:
    """moment가 들어 있는, 또는 다음에 오는 연속 한가한 구간 (ET).

    금요일 밤과 주말, 월요일 새벽은 하나의 구간이다.
    예를 들어 서머타임에는 토요일 10:00 KST부터 월요일 18:00 KST까지 이어진다.
    """
    et = _as_et(moment)
    step = timedelta(minutes=30)
    if is_off_hours(et):
        cursor = et
        while is_off_hours(cursor - step):
            cursor -= step
        start = _binary_edge(cursor - step, cursor, left_is_off=False)
    else:
        start = next_off_hours_start(et)
    cursor = start
    while is_off_hours(cursor + step):
        cursor += step
    end = _binary_edge(cursor, cursor + step, left_is_off=True)
    return start, end


def kst_night_bounds(moment: datetime) -> tuple[datetime, datetime]:
    """이어진 한가한 구간을 한국 시간으로."""
    start, end = off_hours_span(moment)
    return start.astimezone(KST), end.astimezone(KST)


def describe_window(moment: datetime | None = None) -> str:
    moment = moment or datetime.now(ET)
    start, end = off_hours_span(moment)
    kst_start, kst_end = start.astimezone(KST), end.astimezone(KST)
    et_start, et_end = start.astimezone(ET), end.astimezone(ET)
    et = _as_et(moment)
    dst = et.dst() not in (None, timedelta(0))
    season = "서머타임(EDT)" if dst else "표준시(EST)"
    return (
        f"NCBI 한가한 시간: 미국 동부 주말 종일, 평일 {NIGHT_START_HOUR:02d}:00–"
        f"{NIGHT_END_HOUR:02d}:00 ({season}). "
        f"지금 기준으로 이어진 구간은 {et_start:%m-%d %H:%M}–{et_end:%m-%d %H:%M} ET, "
        f"한국 시간 {kst_start:%m-%d %H:%M}–{kst_end:%m-%d %H:%M} KST 입니다."
    )


def usage_path() -> Path:
    env = os.environ.get("BLASTID_USAGE_FILE", "").strip()
    if env:
        return Path(env)
    return Path.home() / ".blastid" / "ncbi_usage.json"


class UsageLog:
    """최근 24시간에 제출한 RID(검색) 수."""

    def __init__(self, path: Path | None = None):
        self.path = Path(path) if path else usage_path()

    def _load(self) -> list[dict]:
        if not self.path.exists():
            return []
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return []
        if isinstance(data, dict):
            data = data.get("searches") or []
        return [x for x in data if isinstance(x, dict) and "ts" in x]

    def _save(self, items: list[dict]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"searches": items}
        tmp = self.path.with_suffix(self.path.suffix + ".tmp")
        tmp.write_text(json.dumps(payload), encoding="utf-8")
        tmp.replace(self.path)

    def recent(self, now: float | None = None) -> list[dict]:
        now = time.time() if now is None else now
        return [x for x in self._load() if now - float(x["ts"]) < WINDOW_SECONDS]

    def count(self, now: float | None = None) -> int:
        return len(self.recent(now))

    def add(self, rid: str = "", ts: float | None = None) -> None:
        now = time.time() if ts is None else ts
        items = self.recent(now)
        items.append({"ts": now, "rid": rid})
        self._save(items)

    def approaching_warnings(self, additional: int, now: float | None = None) -> list[str]:
        now = time.time() if now is None else now
        used = self.count(now)
        total = used + max(0, additional)
        if total < APPROACH_AT:
            return []
        msgs = [
            f"최근 24시간 BLAST 검색 {used}건, 이번 제출 {additional}건이면 합계 {total}건입니다. "
            f"24시간에 {DAILY_SOFT_LIMIT}건에 가까우면 NCBI가 느린 대기열로 옮기거나 막을 수 있습니다."
        ]
        if total > DAILY_SOFT_LIMIT:
            msgs.append(
                f"24시간 {DAILY_SOFT_LIMIT}건을 넘깁니다. "
                "NCBI 사용 지침상 느린 대기열 또는 차단 대상입니다."
            )
        return msgs


@dataclass
class SubmissionPlan:
    searches: int
    warnings: list[str] = field(default_factory=list)
    wait_until: datetime | None = None


def remote_batch_budget(used: int, n_batches: int, limit: int = DAILY_SOFT_LIMIT, backend: str = "auto") -> int:
    """앞으로 원격 RID로 보낼 수 있는 묶음 수.

    한 번 더 제출하면 24시간 합계가 ``limit`` 을 넘을 때부터는 0이다.
    ``remote`` 는 호출자가 원격만 쓰라고 한 경우이고, ``local`` 은 전부 로컬이다.
    """
    if backend == "local":
        return 0
    if backend == "remote":
        return n_batches
    return min(n_batches, max(0, int(limit) - int(used)))


def plan_submission(
    n_searches: int,
    off_hours: bool,
    now: datetime,
    usage: UsageLog | None = None,
) -> SubmissionPlan:
    usage = usage or UsageLog()
    warnings = list(usage.approaching_warnings(n_searches, now.timestamp()))
    wait_until = None
    if n_searches > MANY_SEARCHES:
        text = (
            f"이번에 {n_searches}건의 BLAST 검색(RID 묶음)을 제출합니다. "
            f"{MANY_SEARCHES}건을 넘으면 NCBI는 주말 또는 평일 미국 동부 "
            f"{NIGHT_START_HOUR:02d}:00–{NIGHT_END_HOUR:02d}:00에 실행하라고 합니다. "
            + describe_window(now)
        )
        if not is_off_hours(now):
            if off_hours:
                wait_until = next_off_hours_start(now)
                kst = wait_until.astimezone(KST)
                text += (
                    f" 한가한 시간까지 기다립니다: {wait_until:%Y-%m-%d %H:%M} ET"
                    f" / {kst:%Y-%m-%d %H:%M} KST."
                )
            else:
                text += " `--off-hours`(또는 작업창의 같은 옵션)를 켜면 그때까지 기다립니다."
        warnings.append(text)
    return SubmissionPlan(n_searches, warnings, wait_until)


def sleep_until(when: datetime, log=print, tick: float = 60.0) -> None:
    while True:
        now = datetime.now(when.tzinfo or ET)
        remain = (when - now).total_seconds()
        if remain <= 0:
            log("NCBI 한가한 시간에 들어왔습니다. 제출을 시작합니다.")
            return
        log(
            f"NCBI 한가한 시간까지 약 {remain / 3600:.1f}시간 대기 "
            f"({when.astimezone(ET):%Y-%m-%d %H:%M} ET / {when.astimezone(KST):%Y-%m-%d %H:%M} KST)"
        )
        time.sleep(min(tick, max(1.0, remain)))
