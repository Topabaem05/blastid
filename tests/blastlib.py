"""테스트용 BLAST JSON2_S / XML 묶음. 네트워크를 쓰지 않는다."""
from __future__ import annotations


def json_hit(spec: dict) -> dict:
    acc = spec["acc"]
    bare = acc.split(".")[0]
    return {
        "len": spec["length"],
        "description": [
            {
                "id": f"gi|1|gb|{acc}|",
                "accession": bare,
                "sciname": spec["name"],
                "taxid": spec["taxid"],
                "title": spec.get("title") or spec["name"],
            }
        ],
        "hsps": [
            {
                "bit_score": h["bit"],
                "score": 1,
                "evalue": h["evalue"],
                "identity": h["ident"],
                "align_len": h["alen"],
                "query_from": h["qfrom"],
                "query_to": h["qto"],
            }
            for h in spec["hsps"]
        ],
    }


def blast_json(queries: list[tuple]) -> dict:
    """queries: (label, seq, [hit spec, ...])"""
    reports = []
    for label, seq, hits in queries:
        reports.append(
            {
                "report": {
                    "program": "blastn",
                    "version": "BLASTN 2.16.0+",
                    "search_target": {"db": "core_nt"},
                    "params": {
                        "expect": 0.05,
                        "filter": "L;m;",
                        "sc_match": 1,
                        "sc_mismatch": -2,
                    },
                    "results": {
                        "search": {
                            "query_title": label,
                            "query_len": len(seq),
                            "hits": [json_hit(h) for h in hits],
                        }
                    },
                }
            }
        )
    return {"BlastOutput2": reports}


def blast_xml(queries: list[tuple]) -> str:
    """각 쿼리의 1위 hit(첫 spec)만 XML로. 검증 파서는 1위 hit 만 본다."""
    parts = [
        "<?xml version='1.0'?>",
        "<BlastOutput><BlastOutput_iterations>",
    ]
    for _label, seq, hits in queries:
        parts.append("<Iteration>")
        parts.append(f"<Iteration_query-len>{len(seq)}</Iteration_query-len>")
        parts.append("<Iteration_hits>")
        if hits:
            h = hits[0]
            acc = h["acc"]
            bare = acc.split(".")[0]
            parts.append("<Hit>")
            parts.append(f"<Hit_id>gi|1|gb|{acc}|</Hit_id>")
            parts.append(f"<Hit_accession>{bare}</Hit_accession>")
            parts.append(f"<Hit_len>{h['length']}</Hit_len>")
            parts.append("<Hit_hsps>")
            for hp in h["hsps"]:
                parts.append("<Hsp>")
                parts.append(f"<Hsp_bit-score>{hp['bit']}</Hsp_bit-score>")
                parts.append(f"<Hsp_evalue>{hp['evalue']}</Hsp_evalue>")
                parts.append(f"<Hsp_query-from>{hp['qfrom']}</Hsp_query-from>")
                parts.append(f"<Hsp_query-to>{hp['qto']}</Hsp_query-to>")
                parts.append(f"<Hsp_identity>{hp['ident']}</Hsp_identity>")
                parts.append(f"<Hsp_align-len>{hp['alen']}</Hsp_align-len>")
                parts.append("</Hsp>")
            parts.append("</Hit_hsps></Hit>")
        parts.append("</Iteration_hits></Iteration>")
    parts.append("</BlastOutput_iterations></BlastOutput>")
    return "".join(parts)
