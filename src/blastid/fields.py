"""템플릿에 넣을 수 있는 결과 필드.

숫자는 웹 BLAST 표를 엑셀에 넣던 방식 그대로다.
Query Cover·Per. Ident·겹침 일치도는 0–1 분수(셀 서식 ``0%`` / ``0.00%``).
E value 표시가 ``0.0`` 이면 0.
"""
from __future__ import annotations

from typing import Any

FIELD_HELP: dict[str, str] = {
    "no": "순번 (새 워크북을 만들 때만 기록)",
    "id": "시료 ID",
    "sequence": "BLAST에 제출한 서열",
    "scientific_name": "1위 Scientific Name",
    "max_score": "Max Score (bit score를 웹처럼 내림)",
    "total_score": "Total Score",
    "query_cover": "Query Cover (0–1, 권장 서식 0%)",
    "evalue": "E value (표시가 0.0이면 0)",
    "per_ident": "Per. Ident (0–1, 권장 서식 0.00%)",
    "acc_len": "Accession 서열 길이",
    "accession": "Accession.version",
    "lineage": "Lineage (Fungi부터)",
    "taxid": "NCBI taxid",
    "rid": "BLAST RID",
    "rid_url": "RID를 웹에서 여는 링크",
    "flag": "[이상] 또는 [애매]",
    "comment": "검토 코멘트",
    "screenshot": "스크린샷 시트로 가는 링크",
    "alternatives": "1위가 아닌 hit 요약",
    "overlap_len": "정·역방향 겹침 길이 (bp)",
    "overlap_ident": "겹침 일치도 (0–1, 권장 서식 0.00%)",
    "n_conflicts": "겹침 구간에서 염기가 다른 자리 수",
    "n_resolved": "품질 차이로 한쪽을 고른 불일치 수",
    "n_ambiguities": "합의 서열에 남은 IUPAC 모호 염기 수",
    "n_indels": "겹침 내부 indel 수",
    "merge_mode": "consensus / fallback / single",
    "merge_backend": "biopython / tracy / none",
    "trimmed_len": "트리밍(합침) 후 길이",
    "raw_len": "원본 크로마토그램 길이",
    "trim_mean_q": "트리밍 후 평균 Phred Q",
    "merge_warning": "합침 경고 문장",
    "direction": "F / R / FR",
    "qlen": "BLAST query 길이",
    "source": "NCBI 웹 core_nt RID 또는 로컬 BLAST+ 데이터베이스",
}

# 기존 워크북에 다시 쓸 때 입력 열(ID, Seq 등)은 덮어쓰지 않는다.
INPUT_FIELDS = {"no", "id", "sequence"}
LINK_FIELDS = {"rid_url", "screenshot"}


def build_field_values(
    *,
    no: int | None = None,
    sample_id: str | None = None,
    sequence: str | None = None,
    blast: dict | None = None,
    level: str | None = None,
    comment: str = "",
    merge: Any = None,
) -> dict[str, Any]:
    """한 행에 들어갈 필드 → 값. 스크린샷 링크 문자열은 나중에 채운다."""
    from blastid.scoring import alternatives_text, flag_label

    hits = (blast or {}).get("hits") or []
    top = hits[0] if hits else None
    values: dict[str, Any] = {
        "no": no,
        "id": sample_id,
        "sequence": sequence,
        "flag": flag_label(level),
        "comment": comment or None,
        "screenshot": None,
        "trimmed_len": len(sequence) if sequence else None,
        "qlen": (blast or {}).get("qlen"),
        "source": (blast or {}).get("source"),
    }
    if top is None and blast is not None:
        values.update(
            scientific_name="No significant similarity found",
            max_score=None,
            total_score=None,
            query_cover=None,
            evalue=None,
            per_ident=None,
            acc_len=None,
            accession=None,
            lineage=None,
            taxid=None,
            rid=blast.get("rid"),
            rid_url=blast.get("rid_url"),
            alternatives=None,
        )
    elif top is not None:
        values.update(
            scientific_name=top["name"],
            max_score=top["max"],
            total_score=top["total"],
            query_cover=top["qc"],
            evalue=top["evalue"],
            per_ident=top["pi"],
            acc_len=top["acc_len"],
            accession=top["accession"],
            lineage=(blast or {}).get("lineage") or "",
            taxid=top.get("taxid"),
            rid=blast.get("rid"),
            rid_url=blast.get("rid_url"),
            alternatives=alternatives_text(hits) or None,
        )
    if merge is not None:
        values.update(merge.field_values())
        if merge.trimmed_len:
            values["trimmed_len"] = merge.trimmed_len
    return values
