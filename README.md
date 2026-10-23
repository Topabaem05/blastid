# blastid

지의류·균류 ITS Sanger 서열을 NCBI BLAST로 동정하고, 웹 BLAST Descriptions 표의 **1위 hit** 을 엑셀에 그대로 옮기는 도구입니다.

엑셀·CSV 서열 열, FASTA, `.ab1` 폴더(정방향·역방향 합침 포함)를 한 명령으로 처리합니다. Excel 작업창에서도 같은 서버로 돌릴 수 있습니다. 설정은 웹사이트의 blastn 기본값(core_nt, megablast)과 같고, 숫자를 넣는 방식도 웹 표와 같습니다. 기록한 뒤에는 같은 RID를 XML로 다시 받아 셀마다 대조합니다.

## 설치

```bash
pip install -e ".[dev]"
blastid --help
```

Python 3.10 이상이 필요합니다. 의존성은 Biopython, openpyxl, PyYAML, matplotlib, Pillow, FastAPI, uvicorn 입니다.

NCBI에 연락 이메일을 알리고 싶으면 환경변수만 설정합니다. **메일을 안 넣어도 동작하고, 코드에 메일이 박혀 있지 않습니다.** `NCBI_EMAIL` / `NCBI_API_KEY` 가 있을 때만 요청에 실립니다. API 키는 Taxonomy(E-utilities) 조회에만 쓰입니다.

```bash
export NCBI_EMAIL=you@example.com   # 선택
export NCBI_API_KEY=...             # 선택. https://www.ncbi.nlm.nih.gov/account/
```

## 한 명령으로 돌리기

입력이 엑셀인지, FASTA인지, `.ab1` 인지 자동으로 구분합니다.

```bash
# 엑셀: ID·Seq 헤더를 찾아 BLAST하고 결과 열을 채운 새 파일
blastid run 시료.xlsx -o 결과.xlsx

# 크로마토그램 폴더: Mott 트리밍 → 정·역 합침 → BLAST → 템플릿 엑셀
blastid run ./ab1폴더 -o 결과.xlsx

# FASTA·CSV (utf-8, BOM, cp949)
blastid run seqs.fasta -o 결과.xlsx
blastid run 시료.csv -o 결과.xlsx

# 대상만 확인 (제출하지 않음)
blastid run 시료.xlsx --rows 2-10 --only-empty --dry-run
```

`blastid 시료.xlsx` 처럼 하위 명령 없이 경로만 줘도 `run` 과 같습니다.

## 엑셀

기본 템플릿은 아래 머리글을 가정합니다. 기존 파일에 쓸 때는 **열 글자(A, K…)가 아니라 머리글 문구**로 칸을 찾습니다. 순서가 바뀌어 있어도 됩니다.

`No`, `Collection No.`, `ID(Voucher)`, `Lichen name(1차 가동정)`, `Substratum`, `과제`, `ELF`, `납품`, `DNA`, `Seq`, `Scientific Name`, `Max Score`, `Total Score`, `Query Cover`, `E value`, `Per. Ident`, `Acc. Len`, `Accession`, `Lineage`, `검토 코멘트 / 스크린샷`, `스크린샷 (링크)`

```bash
blastid excel 시료.xlsx --rows 62-102 --bold-only -o 결과.xlsx
blastid excel 시료.xlsx --only-empty          # Scientific Name 이 빈 행만
blastid excel 시료.xlsx --rows 3,24,62-70 --dry-run
blastid excel 시료.xlsx --in-place            # 원본을 덮어쓸 때만
blastid excel 시료.xlsx --id-header "시료번호" --seq-header "ITS"
```

- `--bold-only`: ID 칸이 굵은 행만
- `--only-empty`: Scientific Name 이 비어 있는 행만
- `--clear-old-comments`: 이번에 문제가 없어진 행의 예전 검토 코멘트를 지움. 이 옵션이 없으면, blastid 가 단 메모(작성자 blastid)만 지우고 사람이 적어 둔 코멘트는 남깁니다.
- 서열이 아니거나 50bp 미만(`--min-seq-len`)인 칸은 건너뜁니다.
- 결과 파일이 없으면 `blast_runs/<날짜_시간>/<이름>_blast.xlsx` 에 저장합니다. 입력 파일을 덮어쓰려면 `--in-place` 가 필요합니다.

## .ab1 정·역방향 합치기

파일 이름(확장자 제외)으로 짝을 찾습니다. 기본으로 알아듣는 예:

| 정방향 | 역방향 |
|---|---|
| `ID_F`, `ID-F`, `ID.F`, `ID_fwd`, `ID_forward` | `ID_R`, `ID-R`, `ID_rev`, `ID_reverse` |
| `ID_ITS1F`, `ID_ITS5`, `ID_ITS1`, `ID_LR1` | `ID_ITS4`, `ID_ITS4R`, `ID_NL4` |

규칙은 바꿀 수 있습니다. 정규식에는 `(?P<id>...)` 와 `(?P<dir>...)` 가 있어야 하고, `--pair-regex` 를 주면 기본 규칙 대신 그 식만 씁니다.

```bash
blastid run ./ab1 --pair-regex '^(?P<id>.+?)_(?P<dir>F|R)$'
```

이름 규칙이 안 맞으면 짝 CSV 를 줍니다. 경로는 CSV 파일이 있는 폴더 기준입니다.

```csv
id,forward,reverse
Spec001,traces/weird_left.ab1,traces/weird_right.ab1
```

```bash
blastid run ./ab1 --pairs pairs.csv -o 결과.xlsx
```

처리 순서:

1. 각 리드를 Mott 알고리즘으로 품질 트리밍합니다 (`--trim-cutoff`, 기본 0.05). 품질값이 전부 0인 오래된 장비 파일은 자르지 않습니다. `--trimmer biopython` 이면 SeqIO `abi-trim` 을 씁니다.
2. 역방향 리드를 상보·역전합니다.
3. Biopython `PairwiseAligner` 로 끝 갭은 공짜, 내부 갭은 벌점인 overlap 정렬을 합니다 (match 2, mismatch −5, gap open −3, extend −1).
4. 두 염기가 다르면 Phred 차이가 `--conflict-q`(기본 5) 이상일 때 높은 쪽을 고르고, 아니면 IUPAC 모호 부호를 남깁니다. 겹침 한가운데의 indel 은 `--indel-q`(기본 15) 이상일 때만 염기를 남깁니다.
5. 겹침이 `--min-overlap`(기본 50bp) 보다 짧거나, A/C/G/T 가 비교된 자리의 일치도가 `--min-overlap-ident`(기본 90%) 보다 낮으면 합치기를 실패로 보고 **평균 Q 가 더 높은 리드 하나**를 씁니다. 코멘트는 `[애매]` 로 표시됩니다.
6. 짝이 없는 리드는 그대로 BLAST 합니다. 트리밍 후 `--min-len`(기본 150bp) 미만은 BLAST 하지 않습니다.

합침 통계(겹침 길이, 일치도, 충돌, 모호 염기, 사용한 백엔드)는 `merge_qc.csv` 와 로그에 남습니다. 엑셀에도 넣으려면 템플릿에 `overlap_len`, `overlap_ident` 같은 필드를 추가하면 됩니다.

`blastid ab1 ./폴더` 는 합치고 FASTA·QC 표만 만듭니다. BLAST까지 하려면 `--blast` 를 붙이거나 `blastid run` 을 쓰세요.

### tracy (선택)

[gear-genomics/tracy](https://github.com/gear-genomics/tracy) 의 `consensus` 를 쓸 수 있습니다. 정적 바이너리를 PATH 에 두면 됩니다. 저장소에는 바이너리를 넣지 않습니다.

```bash
blastid run ./ab1 --merger tracy -o 결과.xlsx
```

tracy 는 트레이스에서 염기를 다시 부르고 방향도 스스로 고릅니다. 기본 백엔드(biopython)는 크로마토그램에 들어 있는 basecall 과 품질값을 씁니다. tracy 가 없거나 실패하거나 겹침 기준을 못 넘기면 biopython 결과(또는 더 좋은 단일 리드)로 돌아갑니다.

## FASTA

```bash
blastid fasta seqs.fasta -o 결과.xlsx
```

헤더의 첫 토큰이 ID가 됩니다. 결과 엑셀은 템플릿으로 새로 만듭니다.

## 출력 템플릿

어느 필드가 어느 머리글·서식으로 들어갈지는 YAML이 정합니다. xlsx 는 그 머리글이 이미 스타일된 빈 시트입니다.

```bash
blastid init-template -o ./내템플릿
# ./내템플릿/template.yaml
# ./내템플릿/template.xlsx
blastid run 시료.xlsx --template ./내템플릿/template.yaml -o 결과.xlsx
```

기존 엑셀에 쓸 때는 YAML의 `header` 또는 `aliases` 와 같은 머리글을 찾아 값을 넣습니다. 템플릿에만 있고 시트에 없는 결과 열은 오른쪽에 추가합니다. ID·Seq 같은 입력 열은 덮어쓰지 않습니다.

새 파일을 만들 때(ab1, FASTA)는 템플릿 xlsx 를 복사한 뒤 2행부터 채웁니다.

필드 예: `id`, `sequence`, `scientific_name`, `max_score`, `total_score`, `query_cover`, `evalue`, `per_ident`, `acc_len`, `accession`, `lineage`, `taxid`, `rid`, `rid_url`, `flag`, `comment`, `screenshot`, `alternatives`, `overlap_len`, `overlap_ident`, `n_conflicts`, `n_resolved`, `n_ambiguities`, `n_indels`, `merge_mode`, `merge_backend`, `trimmed_len`, `raw_len`, `trim_mean_q`, `merge_warning`, `direction`, `qlen`.

Query Cover·Per. Ident·겹침 일치도는 **0–1 분수**로 저장합니다. 권장 서식은 Query Cover `0%`, Per. Ident `0.00%` 입니다. E value 표시가 `0.0` 이면 숫자 0을 넣습니다.

서식(기본 템플릿): 머리글 진한 파랑(`1F4E78`), 틀 고정 `D2`, 자동 필터, 짝수 행 줄무늬, 결과 칸 연두, 문제 행은 빨강/노랑, 데이터 행 높이 30.

## 값 계산과 검증

| 단계 | 방법 |
|---|---|
| 제출·수신 | NCBI BLAST URL API `Blast.cgi` `CMD=Put` / `CMD=Get`, `FORMAT_TYPE=JSON2_S` |
| 설정 | `core_nt`, megablast, expect 0.05, word 28, reward/penalty 1/−2, linear gap, filter `L;m;`, hit 100. RID의 `params` 가 다르면 경고 |
| Max / Total Score | bit score 를 웹처럼 99.9 초과에서 소수점 버림 |
| Query Cover | HSP query 구간의 합집합 ÷ query 길이, 반올림(0.5는 올림). 100%를 넘으면 99% |
| Per. Ident | bit score 가 가장 높은 HSP의 identity ÷ align_len, 소수 둘째 자리 |
| E value | 웹과 같은 자릿수 문자열. `0.0` 은 0 |
| Lineage | E-utilities taxonomy. `Fungi` 부터 |
| 검증 | 저장한 엑셀을 다시 열고, 같은 RID를 **XML로 다시 받아 다른 파서**로 계산하고 taxonomy 도 다시 받아 셀과 비교 |

서열은 기본 10개씩 한 RID로 보냅니다. 제출 간격은 10초 이상, 같은 RID 상태 확인은 처음 20초(또는 RTOE) 뒤부터 1분에 한 번입니다.

각 행의 RID와 웹 링크가 `blast_runs/<날짜_시간>/log.csv`, `log.json`, `run.log` 에 남습니다.

```text
https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=<RID>
```

NCBI는 RID를 약 36시간 보관합니다. 그 안에 `--reuse-rids RID1,RID2` 로 다시 제출하지 않고 엑셀만 다시 만들 수 있습니다. RID 개수는 묶음 수와 같아야 합니다.

```bash
blastid verify 결과.xlsx --log blast_runs/<날짜_시간>/log.json
```

## 문제 행

- **[이상] (빨강)**: 1위 Per. Ident < 90% 또는 Query Cover < 80% (`--min-ident`, `--min-cover`), 또는 hit 가 없음
- **[애매] (노랑)**:
  - 1위와 Max Score·Ident 가 같은 다른 종
  - 90%·Cover 80% 를 넘는 다른 속 hit 가 1위와 Ident 차이가 `--alt-margin`(기본 2%p) 이내
  - 1위가 sp. / aff. / cf. / uncultured 처럼 종이 안 적힌 서열
  - 정·역 합침이 실패해 단일 리드로 대체했거나, 합의 서열의 모호 염기가 3개 이상

문제 행은 검토 코멘트(셀 메모 포함)와 `스크린샷` 시트의 상위 15 hit 표 이미지, 그 이미지로 가는 링크를 갖습니다.

## 로그

`blast_runs/<날짜_시간>/`

- `log.csv`, `log.json` — RID, 웹 링크, 표시값, 판정, 검증, 상위 15 hit, 합침 통계
- `raw/<RID>.json` — NCBI 원본
- `shots/*.png` — 문제 행 표
- `run.log`
- ab1 인 경우 `merged.fasta`, `merged.fastq`, `ab1_qc.csv`, `merge_qc.csv`

## CSV

`.csv` 도 엑셀과 같이 첫 줄을 머리글로 읽습니다. 인코딩은 UTF-8(BOM 포함)을 먼저 시도하고, 아니면 CP949입니다. 구분자 `,` `;` 탭을 구분합니다. 서식과 그림은 CSV에 담을 수 없어서 결과는 `.xlsx` 로 저장합니다.

```bash
blastid run 시료.csv --id-header ID --seq-header Seq -o 결과.xlsx
```

## Excel 작업창

브라우저에서 NCBI를 직접 치면 CORS와 호출 제한에 걸리므로, 작업창은 이 컴퓨터의 `blastid serve` 로만 말합니다. 서버가 BLAST를 하고 결과를 머리글 이름으로 돌려주면, 작업창이 시트에 쓰고 `[이상]`/`[애매]` 색과 첫 행 고정을 적용합니다.

```bash
blastid serve
# https://localhost:8787/taskpane.html
# 매니페스트: src/blastid/addin/manifest.xml  (주소는 https://localhost:8787)
```

`addin/` 은 패키지 데이터로 들어갑니다. `pip install`(편집 설치가 아닌 경우)로 받아 저장소 밖에서 `blastid serve` 를 실행해도 `/taskpane.html` 을 제공합니다. 그 파일이 없으면 404로 넘기지 않고, 시작할 때 한국어 오류로 멈춥니다.

처음 실행하면 `~/.blastid/dev-certs/localhost.crt` 를 만듭니다. 이미 `npx office-addin-dev-certs install` 로 만들어 둔 인증서가 있으면 그것을 씁니다. Excel이 인증서를 신뢰해야 작업창이 열립니다. 브라우저에서 `https://localhost:8787/taskpane.html` 을 한 번 열어 경고를 확인한 뒤, 인증서를 신뢰할 루트에 넣으세요.

작업창에서 고를 수 있는 것:

- 현재 시트 또는 선택한 행, 행 범위(`2-10,15`), 굵은 ID만, Scientific Name 이 빈 행만
- ID·서열 머리글, 기본 템플릿 또는 YAML 파일, 최소 Ident 90 / Cover 80 / 다른 속 차이 2, 묶음 크기
- 이메일(기본 빈 칸, 비어 있으면 요청에 넣지 않음), 50건을 넘길 때 한가한 시간까지 기다리기
- 실행 위치(자동 / NCBI 웹 / 로컬 BLAST+), 로컬 DB, RID 최대 대기(분, 기본 60). 최근 24시간 NCBI 검색 수도 같이 보여 줍니다.

설정은 Office 로밍 설정에 저장하고, 그게 없으면 브라우저 `localStorage` 에 남깁니다. 진행 목록에는 RID와 NCBI 결과 페이지 링크가 나옵니다.

사이드로드:

- **Windows Excel**: 삽입 → 추가 기능 → 내 추가 기능 → 내 추가 기능 업로드 → `src/blastid/addin/manifest.xml`. 홈 탭에 **BLAST 동정** 단추가 생깁니다.
- **Mac Excel**: 같은 업로드 메뉴가 있으면 그것을 씁니다. 없으면 `manifest.xml` 을 `~/Library/Containers/com.microsoft.Excel/Data/Documents/wef/` 에 넣고 Excel을 다시 엽니다.
- **웹 Excel** (excel.office.com): 삽입 → 추가 기능 → 추가 기능 업로드 → `manifest.xml`. 이 컴퓨터의 브라우저가 `https://localhost:8787` 을 열 수 있어야 합니다. 서버는 `127.0.0.1` 만 받습니다.

`yo office` 로 만드는 작업창과 같은 Office.js 매니페스트(XML, 홈 탭 단추, 작업창)입니다. 페이지는 webpack 개발 서버 대신 `blastid serve` 가 직접 제공합니다.

작업창 경로의 독립 XML 재검증은 하지 않습니다. 셀을 다시 읽어 RID XML과 맞추는 검증은 `blastid run` / `blastid verify` 가 합니다.

## NCBI 사용 지침

[BLAST 개발자 안내](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html)와 [NCBI 정책](https://www.ncbi.nlm.nih.gov/home/about/policies/)을 따릅니다.

- 모든 BLAST·E-utilities 요청에 `tool=blastid` 를 붙입니다. `email` 은 `NCBI_EMAIL` 또는 작업창에 입력한 주소가 있을 때만 붙입니다.
- 제출 간격은 10초 이상, 같은 RID 확인은 1분에 한 번 이하, E-utilities는 초당 3회 이하입니다.
- 검색 1건은 새로 제출한 RID 1회입니다. `--reuse-rids` 로 이미 있는 RID를 다시 열 때는 로그에 `기존 RID 재사용` 이라고 적고, 24시간 횟수에 넣지 않습니다. 한 번에 50건을 넘기면 경고합니다. `--off-hours` 또는 작업창의 같은 옵션이 켜져 있고 지금이 한가한 시간이 아니면, 주말 또는 평일 미국 동부 21:00–05:00(서머타임 반영)까지 기다립니다. 금요일 밤·주말·월요일 새벽은 하나의 구간으로 안내합니다. 서머타임에는 토요일 10:00 KST부터 월요일 18:00 KST까지이고, 표준시에는 토요일 11:00 KST부터 월요일 19:00 KST까지입니다. 주중 밤은 8시간(서머타임 기준 다음날 10:00–18:00 KST)입니다.
- RID가 `WAITING` 으로 남으면 `--max-wait`(기본 60분, 작업창의 같은 항목)까지 기다립니다. 그 안에 READY가 아니면 한국어 오류와 RID 링크를 남기고 멈춥니다. 같은 RID는 나중에 `--reuse-rids` 로 이어서 볼 수 있습니다. `--backend auto` 이면 그 묶음만 로컬 BLAST+로 넘길 수 있고, BLAST+나 DB가 없으면 이미 끝난 행만 저장한 뒤 멈춥니다.
- 최근 24시간 검색 수는 `~/.blastid/ncbi_usage.json` 에 남습니다. 이번에 더하면 80건 이상이면 100건에 가깝다고 경고하고, 100건을 넘으면 NCBI가 느리게 하거나 막을 수 있다고 알립니다. 로컬 BLAST+ 검색은 이 횟수에 넣지 않습니다.
- 아래 고지는 작업창과 실제 제출 로그에 표시합니다.

> BLAST는 NCBI의 공개 서비스입니다. 미국 정부는 NCBI가 제공하는 문서와 소프트웨어의 정확성·완전성·유용성을 보증하지 않으며 법적 책임을 지지 않습니다. NCBI는 특정 상업 제품·방법·서비스를 보증하거나 권하지 않습니다. 웹사이트에 실린 저자의 견해가 미국 정부의 입장을 뜻하지는 않으며, 광고나 제품 보증에 쓸 수 없습니다.
> https://www.ncbi.nlm.nih.gov/home/about/policies/
> https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html

## 로컬 BLAST+ (24시간 100건)

원격 검색은 `~/.blastid/ncbi_usage.json` 에 RID 하나씩 남습니다. 다음 RID를 내면 최근 24시간이 100건(`--remote-limit`)을 넘게 되면, 그 묶음부터는 이 컴퓨터의 NCBI BLAST+로 돌립니다. 주소를 바꾸거나 한도를 피하지 않습니다. `--backend remote` 는 경고만 하고 계속 NCBI에 보내고, `--backend local` 은 처음부터 로컬만 씁니다. 기본은 `auto` 입니다.

기본 데이터베이스는 NCBI `ITS_RefSeq_Fungi` 이고, 학명·taxid용 `taxdb` 를 같이 받습니다. 다른 NCBI DB(`ITS_eukaryote_sequences` 등)나 직접 만든 FASTA(`makeblastdb`)도 `--local-db` 로 지정합니다. 로컬 megablast 설정은 웹과 같습니다. word 28, reward 1 / penalty -2, gap 0/0, dust·소문자 마스크, evalue 0.05, 대상 100개. 표의 숫자는 웹과 같은 규칙으로 계산하고, taxid가 있으면 lineage는 지금까지처럼 E-utilities(`tool=blastid`, 초당 3회 이하)로 가져옵니다.

로컬 행의 `검색 출처` 는 `로컬 BLAST+ (ITS_RefSeq_Fungi)` 처럼 적히고, 결과 칸은 웹 결과의 연두색과 다른 파란색입니다. 코멘트에는 core_nt와 DB가 달라 hit이 다를 수 있다는 문장이 붙습니다. `log.json` 의 `search_backend` 가 `local` 이면 `blastid verify` 는 RID를 받지 않고 같은 로컬 DB로 다시 검색합니다.

BLAST+나 DB가 없으면 그 시점에서 멈추고, 이미 끝난 행은 결과 파일에 남깁니다. 이어서 하려면 `blastid setup-local` 뒤에 그 파일을 `--only-empty` 로 다시 실행합니다.

```bash
blastid setup-local                          # BLAST+ , ITS_RefSeq_Fungi, taxdb
blastid setup-local --db ITS_eukaryote_sequences
blastid setup-local --fasta unite.fasta --name unite
blastid run 시료.xlsx --backend auto --local-db ITS_RefSeq_Fungi
```

2026-10-23에 받은 `ITS_RefSeq_Fungi.tar.gz` 는 68,941,961바이트(약 66MB)입니다. 풀어 본 구성 파일 합은 303,702,490바이트(약 290MB)입니다. NCBI가 적은 bytes-total 9,103,975(서열 20,216개)는 핵산 인덱스만이고, 같은 압축 안에 `taxdb.btd`, `taxdb.bti`, `taxonomy4blast.sqlite3` 가 같이 들어 있습니다. 별도 `taxdb.tar.gz` 는 65,544,367바이트이고 NCBI가 적은 용량은 294,519,351바이트(약 281MB)입니다. `setup-local` 이 taxdb를 다시 받아도 같은 파일이면 이 용량을 덮어씁니다. BLAST+ 리눅스 x64 실행 파일 묶음은 약 282MB입니다.

## 테스트

```bash
pytest
```

NCBI 응답은 흉내 내고, `.ab1` 은 합성 크로마토그램을 씁니다. 실제 NCBI를 보려면:

```bash
BLASTID_LIVE_XLSX=/path/to/sample.xlsx pytest -m live
```

실서열이 들어 있는 표는 저장소에 넣지 않습니다.

## 한계

- core_nt 가 갱신되면 1위 hit·퍼센트·accession 이 날짜에 따라 달라질 수 있습니다. 그때의 기준은 그 RID를 웹에서 연 결과입니다.
- 1위 이름은 BLAST가 보여 준 이름입니다. Taxonomy 의 현재 학명과 다르면 로그에만 적습니다.
- 합침 일치도는 겹침 구간에서 양쪽이 모두 A/C/G/T 인 자리만 셉니다. N 은 분모에 넣지 않습니다.
- tracy 를 쓰면 basecall 이 크로마토그램 파일에 들어 있는 서열과 달라질 수 있습니다.
- Excel 작업창은 `https://localhost:8787` 인증서를 신뢰해야 열립니다. 서버를 끄면 작업창도 멈춥니다.
- 작업창은 시트에 결과를 쓰지만, RID XML로 셀을 다시 대조하는 검증은 명령줄 `blastid run` 이 합니다.
- 50건을 넘겨 `--off-hours` 로 기다릴 때는 그 프로세스가 켜져 있어야 합니다. 대기는 몇 시간이 될 수 있습니다.
