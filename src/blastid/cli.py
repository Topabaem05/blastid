"""blastid 명령줄."""
from __future__ import annotations

import argparse
import sys

from blastid import __version__
from blastid.pipeline import execute
from blastid.template import init_template


def _blast_args(sp):
    sp.add_argument("--out-dir", help="로그·원본 결과·이미지 폴더 (기본 ./blast_runs/날짜_시간)")
    sp.add_argument("--template", help="출력 템플릿 YAML (기본: 패키지에 들어 있는 default.yaml)")
    sp.add_argument("--sheet", help="시트 이름 (기본은 템플릿의 sheet)")
    sp.add_argument("--batch-size", type=int, default=10, help="한 RID에 넣을 서열 수 (기본 10)")
    sp.add_argument("--min-ident", type=float, default=90.0, help="이 값 미만 Per. Ident 면 [이상] (기본 90)")
    sp.add_argument("--min-cover", type=float, default=80.0, help="이 값 미만 Query Cover 면 [이상] (기본 80)")
    sp.add_argument(
        "--alt-margin",
        type=float,
        default=2.0,
        help="다른 속 후보로 볼 Ident 차이(%%p, 기본 2)",
    )
    sp.add_argument("--no-verify", action="store_true", help="독립 검증(XML·Taxonomy 재조회) 생략")
    sp.add_argument("-o", "--output", help="결과 엑셀 경로")
    sp.add_argument(
        "--reuse-rids",
        help="새로 제출하지 않고 기존 RID를 묶음 순서대로 사용 (쉼표 구분, 약 36시간 이내)",
    )
    sp.add_argument("--dry-run", action="store_true", help="대상만 보여주고 NCBI에 제출하지 않음")
    sp.add_argument("--min-seq-len", type=int, default=50, help="BLAST에 넣을 최소 서열 길이 (기본 50)")


def _excel_args(sp):
    sp.add_argument("--rows", help="처리할 행. 예: 2-6 또는 3,24,62-70")
    sp.add_argument("--bold-only", action="store_true", help="ID 칸이 굵은 글씨인 행만")
    sp.add_argument("--only-empty", action="store_true", help="Scientific Name 이 비어 있는 행만")
    sp.add_argument("--in-place", action="store_true", help="입력 엑셀을 덮어씀 (기본은 새 파일)")
    sp.add_argument(
        "--clear-old-comments",
        action="store_true",
        help="문제가 없어진 행의 기존 검토 코멘트·스크린샷 링크를 지움",
    )
    sp.add_argument("--id-header", help="ID 열 헤더를 직접 지정")
    sp.add_argument("--seq-header", help="서열 열 헤더를 직접 지정")


def _merge_args(sp):
    sp.add_argument(
        "--pair-regex",
        action="append",
        help="파일명(확장자 제외)에서 시료와 방향을 나누는 정규식. "
        "(?P<id>...) 와 (?P<dir>...) 가 필요. 지정하면 기본 규칙을 대체. 여러 번 쓸 수 있음",
    )
    sp.add_argument("--pairs", help="정·역방향 짝 CSV (열: forward, reverse, 선택 id)")
    sp.add_argument("--no-merge", action="store_true", help="짝을 맞추지 않고 리드마다 BLAST")
    sp.add_argument("--merger", choices=["biopython", "tracy"], default="biopython",
                    help="합의 서열 백엔드 (기본 biopython, tracy는 PATH에 tracy가 있을 때)")
    sp.add_argument("--tracy-trim", type=int, default=4, help="tracy consensus -t 트리밍 강도 1–9 (기본 4)")
    sp.add_argument("--trimmer", choices=["mott", "biopython"], default="mott")
    sp.add_argument("--trim-cutoff", type=float, default=0.05, help="Mott 오류확률 cutoff (기본 0.05)")
    sp.add_argument("--min-len", type=int, default=150, help="트리밍 후 BLAST 최소 길이 (기본 150)")
    sp.add_argument("--min-q", type=float, default=25, help="평균 Q 경고 기준 (기본 25)")
    sp.add_argument("--min-overlap", type=int, default=50, help="합침에 필요한 최소 겹침 bp (기본 50)")
    sp.add_argument("--min-overlap-ident", type=float, default=90.0, help="겹침 일치도 %% 미만이면 단일 리드 (기본 90)")
    sp.add_argument("--conflict-q", type=float, default=5, help="불일치를 고를 Phred 차이 (기본 5)")
    sp.add_argument("--indel-q", type=float, default=15, help="겹침 안 indel 염기를 남길 최소 Q (기본 15)")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="blastid",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "지의류/균류 ITS Sanger 서열을 NCBI BLAST(웹과 같은 core_nt megablast)로 동정하고 "
            "엑셀에 웹 표와 같은 값으로 기록합니다."
        ),
        epilog="""예시:
  blastid run 시료.xlsx -o 결과.xlsx
  blastid run ./ab1폴더 -o 결과.xlsx
  blastid run seqs.fasta
  blastid run 시료.xlsx --rows 2-10 --only-empty --dry-run
  blastid ab1 ./ab1폴더 --pairs pairs.csv          # 합치기만
  blastid init-template -o ./내템플릿
  blastid verify 결과.xlsx --log blast_runs/.../log.json

NCBI_EMAIL, NCBI_API_KEY 는 환경변수가 있을 때만 전송합니다. 없어도 동작합니다.
""",
    )
    p.add_argument("--version", action="version", version=f"blastid {__version__}")
    sub = p.add_subparsers(dest="cmd")

    run = sub.add_parser("run", help="엑셀·FASTA·ab1 폴더를 자동으로 알아내 끝까지 BLAST")
    run.add_argument("inputs", nargs="+")
    _blast_args(run)
    _excel_args(run)
    _merge_args(run)

    excel = sub.add_parser("excel", help="엑셀의 서열 열을 BLAST해 결과 열을 채움")
    excel.add_argument("input")
    _blast_args(excel)
    _excel_args(excel)

    ab1 = sub.add_parser("ab1", help=".ab1 품질 트리밍·정역 합침 (BLAST는 --blast)")
    ab1.add_argument("inputs", nargs="+")
    ab1.add_argument("--blast", action="store_true", help="합친 서열을 바로 BLAST")
    _blast_args(ab1)
    _merge_args(ab1)

    fasta = sub.add_parser("fasta", help="FASTA를 BLAST해 템플릿 엑셀을 만듦")
    fasta.add_argument("input")
    _blast_args(fasta)

    verify = sub.add_parser("verify", help="log.json 의 RID로 엑셀 값을 다시 검증")
    verify.add_argument("input")
    verify.add_argument("--log", required=True, dest="log_path")
    verify.add_argument("--sheet")
    verify.add_argument("--template")
    verify.add_argument("--out-dir")

    init = sub.add_parser("init-template", help="기본 YAML·xlsx 템플릿을 만들어 줌")
    init.add_argument("-o", "--output", default="blastid_template", help="만들 폴더 (기본 ./blastid_template)")
    return p


def _ns_from(args, kind, inputs):
    args.kind = kind
    args.inputs = inputs
    return args


def main(argv=None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    else:
        argv = list(argv)
    known = {"run", "excel", "ab1", "fasta", "verify", "init-template"}
    if argv and argv[0] not in known and not argv[0].startswith("-"):
        argv = ["run", *argv]
    parser = build_parser()
    if not argv:
        parser.print_help()
        return
    args = parser.parse_args(argv)
    if args.cmd is None:
        parser.print_help()
        return
    if args.cmd == "init-template":
        from pathlib import Path

        yml, xlsx = init_template(Path(args.output))
        print(f"템플릿을 썼습니다.\n  {yml}\n  {xlsx}")
        return
    if args.cmd == "verify":
        _cmd_verify(args)
        return
    if args.cmd == "excel":
        execute(_ns_from(args, "excel", [args.input]))
    elif args.cmd == "fasta":
        execute(_ns_from(args, "fasta", [args.input]))
    elif args.cmd == "ab1":
        execute(_ns_from(args, "ab1", list(args.inputs)))
    elif args.cmd == "run":
        execute(_ns_from(args, "run", list(args.inputs)))
    else:
        parser.error(f"알 수 없는 명령 {args.cmd}")


def _cmd_verify(args) -> None:
    import json
    from pathlib import Path

    import openpyxl

    from blastid.ncbi import NCBI
    from blastid.pipeline import default_outdir, make_logger, verify
    from blastid.template import load_template
    from blastid.workbook import ensure_fields, map_fields, pick_sheet

    outdir = default_outdir(args.out_dir)
    log = make_logger(outdir)
    template = load_template(args.template)
    if args.sheet:
        template.sheet = args.sheet
    recs = json.loads(Path(args.log_path).read_text(encoding="utf-8"))
    wb = openpyxl.load_workbook(args.input)
    ws = pick_sheet(wb, template.sheet, log)
    cols = map_fields(ws, template)
    cols = ensure_fields(ws, template, cols, log)
    entries = [
        dict(
            row=r["row"],
            rid=r["rid"],
            query_index=r["query_index"],
            seq_len=r.get("seq_len") or r.get("qlen"),
            taxid=r.get("taxid"),
            name_json=r.get("name"),
        )
        for r in recs
    ]
    probs = verify(NCBI(log), ws, template, cols, entries, log)
    log(
        f"검증 결과: {len(entries) - len(probs)}/{len(entries)}행 일치"
        + (f", 불일치 {sorted(probs)}" if probs else "")
    )
    raise SystemExit(1 if probs else 0)


if __name__ == "__main__":
    main()
