"""한가한 시간(미국 동부, 서머타임 포함)과 24시간 검색 수."""
import tempfile
import time
from datetime import datetime, timedelta
from pathlib import Path

from blastid.ncbi import EUTILS_MIN_INTERVAL, pace_eutils
from blastid.policy import (
    ET,
    KST,
    UsageLog,
    describe_window,
    is_off_hours,
    kst_night_bounds,
    next_off_hours_start,
    off_hours_span,
    plan_submission,
)

def test_weekend_window_is_one_kst_span():
    """금요일 밤부터 월요일 새벽까지는 한국 시간으로 한 구간이다."""
    friday_evening = datetime(2026, 10, 23, 20, 0, tzinfo=KST)
    assert friday_evening.weekday() == 4
    start, end = off_hours_span(friday_evening)
    assert start.astimezone(KST) == datetime(2026, 10, 24, 10, 0, tzinfo=KST)
    assert end.astimezone(KST) == datetime(2026, 10, 26, 18, 0, tzinfo=KST)
    text = describe_window(friday_evening)
    assert "10-24 10:00" in text and "10-26 18:00" in text
    # 한겨울(EST)에는 같은 금요일이 11:00–19:00 KST로 밀린다.
    winter_friday = datetime(2026, 1, 16, 20, 0, tzinfo=KST)
    assert winter_friday.weekday() == 4
    wstart, wend = off_hours_span(winter_friday)
    assert wstart.astimezone(KST).hour == 11
    assert wend.astimezone(KST).hour == 19
    assert (wend.astimezone(ET) - wstart.astimezone(ET)) > timedelta(days=2)
    # 주 중 밤은 8시간만.
    wed = datetime(2026, 7, 15, 15, 0, tzinfo=ET)
    nstart, nend = off_hours_span(wed)
    assert nstart.astimezone(ET).hour == 21
    assert nend.astimezone(ET).hour == 5
    assert (nend - nstart) == timedelta(hours=8)


def test_remote_budget_stops_at_the_limit():
    from blastid.policy import remote_batch_budget

    assert remote_batch_budget(99, 5, 100, "auto") == 1
    assert remote_batch_budget(100, 5, 100, "auto") == 0
    assert remote_batch_budget(0, 3, 100, "local") == 0
    assert remote_batch_budget(100, 3, 100, "remote") == 3
    assert remote_batch_budget(0, 150, 100, "auto") == 100


def test_off_hours_window_respects_dst():
    summer_night = datetime(2026, 7, 15, 21, 0, tzinfo=ET)
    summer_day = datetime(2026, 7, 15, 15, 0, tzinfo=ET)
    winter_night = datetime(2026, 1, 14, 21, 0, tzinfo=ET)
    winter_morning = datetime(2026, 1, 14, 4, 30, tzinfo=ET)
    winter_close = datetime(2026, 1, 14, 5, 0, tzinfo=ET)
    saturday = datetime(2026, 7, 18, 12, 0, tzinfo=ET)
    assert summer_night.weekday() < 5
    assert is_off_hours(summer_night)
    assert not is_off_hours(summer_day)
    assert is_off_hours(winter_night)
    assert is_off_hours(winter_morning)
    assert not is_off_hours(winter_close)
    assert saturday.weekday() == 5
    assert is_off_hours(saturday)
    nxt = next_off_hours_start(summer_day)
    assert nxt.astimezone(ET).hour == 21
    assert nxt.astimezone(ET).date() == summer_day.date()
    assert next_off_hours_start(summer_night) == summer_night
    # 21:00 ET → KST. EDT(UTC-4)는 +13시간, EST(UTC-5)는 +14시간.
    assert summer_night.astimezone(KST).hour == 10
    assert winter_night.astimezone(KST).hour == 11
    kst_start, kst_end = kst_night_bounds(summer_day)
    assert kst_start.hour == 10 and kst_end.hour == 18
    text = describe_window(summer_day)
    assert "KST" in text and "10:00" in text and "18:00" in text
    winter_text = describe_window(winter_night)
    assert "11:00" in winter_text and "19:00" in winter_text


def test_plan_waits_only_when_off_hours_enabled():
    afternoon = datetime(2026, 7, 15, 15, 0, tzinfo=ET)
    usage = UsageLog(Path(tempfile.mkdtemp()) / "u.json")
    waiting = plan_submission(51, True, afternoon, usage)
    assert waiting.wait_until is not None
    assert waiting.wait_until.astimezone(ET).hour == 21
    assert any("51" in w for w in waiting.warnings)
    assert any("KST" in w for w in waiting.warnings)
    immediate = plan_submission(51, False, afternoon, usage)
    assert immediate.wait_until is None
    assert any("off-hours" in w or "한가한" in w for w in immediate.warnings)
    small = plan_submission(50, True, afternoon, usage)
    assert small.wait_until is None
    night = plan_submission(80, True, datetime(2026, 7, 15, 22, 0, tzinfo=ET), usage)
    assert night.wait_until is None


def test_usage_log_prunes_and_warns_near_100(tmp_path):
    log = UsageLog(tmp_path / "usage.json")
    log.add(rid="old", ts=time.time() - 90_000)
    assert log.count() == 0
    now = time.time()
    for i in range(79):
        log.add(rid=f"r{i}", ts=now)
    assert log.approaching_warnings(0, now) == []
    warns = log.approaching_warnings(1, now)
    assert warns and "100" in warns[0]
    for i in range(30):
        log.add(rid=f"x{i}", ts=now)
    hard = log.approaching_warnings(0, now)
    assert any("넘깁" in w for w in hard)


def test_eutils_paced_at_three_per_second():
    start = time.perf_counter()
    pace_eutils()
    pace_eutils()
    elapsed = time.perf_counter() - start
    assert elapsed + 0.05 >= EUTILS_MIN_INTERVAL
