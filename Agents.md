# Agents.md

이 문서는 이 저장소를 다루는 에이전트가 `run_pipe.sh` 기준의 실제 실행 파이프라인을 빠르게 이해하도록 정리한 메모다.

## 분석 범위

- 기본 원칙: `.gitignore`에 걸리는 파일과 디렉터리는 코드 이해 대상에서 제외한다.
- 예외: 사용자 요청상 `run_pipe.sh`는 파이프라인 진입점으로 직접 참고했다.
- 따라서 아래 설명은 주로 추적 가능한 비무시 Python 코드와 `run_pipe.sh`의 호출 관계를 기준으로 작성되었다.

무시 대상으로 간주한 대표 패턴:

- `20_acc_pipe`, `30_skype_pipe*`, `PanDepth`,`backup*`, `depth/`
- `*.pkl`, `*.txt`, `*.ipynb`, `*.h5`, `*.png`
- `*.sh` 전반

## 엔트리포인트

실제 실행은 아래 명령을 기준으로 본다.

```bash
conda activate skype
bash run_pipe.sh <CELL_LINE> [OUTPUT_PREFIX]
```

`run_pipe.sh`의 동작 요약:

- `THREAD=64` 고정
- `UTG="utg"` 고정
- 입력은 대부분 `/home/hyunwoo/ACCtools-pipeline/90_skype_run/<CELL_LINE>/...` 아래에서 읽는다
- 출력 prefix를 주지 않으면 `30_skype_pipe/<paf basename>_<HH_MM_SS>` 형식으로 생성한다
- 레퍼런스 BED/FAI는 저장소의 `public_data/`를 사용한다
- `22_save_matrix.py` 호출 시 `--not_use_nclose_weight`를 넘기므로 기본 파이프라인에서는 nclose weight를 쓰지 않는다

## 실제 실행 순서

`run_pipe.sh`는 현재 아래 9단계를 모두 실행한다.

| 단계 | 스크립트 | 핵심 역할 | 주요 산출물 |
|---|---|---|---|
| 0 | `00_depth_norm.py` | 암 샘플 depth를 CHM13 depth 기준으로 정규화 | `*_normalized.win.stat.gz` |
| 1 | `02_Build_Breakend_Graph_Limited.py` | PAF 전처리, `.ppc.paf` 생성, telomere/NClose 기반 그래프 구축, path 탐색 | `<PREFIX>/00_raw`, `.ppc.paf`, 각종 `pkl/txt` |
| 2 | `11_Ref_Outlier_Contig_Modify.py` | type4/type2 outlier 이벤트를 depth 측정용 PAF 쌍으로 분리 | `<PREFIX>/11_ref_ratio_outliers/front_jump`, `back_jump` |
| 3 | `21_run_depth.py` | path/outlier별 PAF를 분할하고 PanDepth 실행용 입력 준비 및 실행 | `<PREFIX>/21_pat_depth`, outlier별 `*.paf`, `contig_pat_vec_data.pkl` |
| 4 | `22_save_matrix.py` | path depth를 행렬화하여 NNLS 입력 생성 | `matrix.h5`, `23_input.pkl`, `tar_chr_data.pkl`, `B.npy` |
| 5 | `23_run_nnls.py` | 기본 NNLS 가중치 계산 | `weight.npy`, `predict_B.npy`, `A_idx_list.pkl` |
| 6 | `24_cluster_weight.py` | Julia 기반 추가 필터링/클러스터링 가중치 계산 | `weight_filter.npy`, `predict_B_filter.npy`, `weight_cluster.npy`, `predict_B_cluster.npy` |
| 7 | `30_virtual_sky.py` | Virtual SKY 그림 생성 | `virtual_sky*.pdf`, `virtual_sky*.png`, `tot_loc_list.pkl` |
| 8 | `31_depth_analysis.py` | Circos, BED, VCF 등 최종 결과 생성 | `total_cov*.pdf/png`, `SKYPE_result*.bed`, `SV_call_result*.vcf` |

## 단계별 상세

### 0. `00_depth_norm.py`

입력:

- 샘플 window depth: `01_depth/<CELL_LINE>.win.stat.gz`
- 레퍼런스 depth: `deps/SKYPE/public_data/CHM13.win.stat.gz`
- censat BED: `public_data/chm13v2.0_censat_v2.1.m.bed`

동작:

- `chrM` 제거
- ref meandepth의 중앙값으로 정규화 offset 계산
- censat 구간과 `chrY`는 보정하지 않음
- `meandepth` 기반으로 `totaldepth`를 다시 계산

출력:

- 원본 depth 파일 옆에 `<CELL_LINE>_normalized.win.stat.gz`

### 1. `02_Build_Breakend_Graph_Limited.py`

이 저장소의 핵심 단계다. 실제로 뒤 단계가 의존하는 대부분의 중간 파일을 여기서 만든다.

입력:

- primary contig align PAF: `<CELL_LINE>.ctg.aln.paf`
- alt unitig align PAF: `<CELL_LINE>.utg.aln.paf`
- 원본 PAF들: `--orignal_paf_loc ...ctg.paf ...utg.paf`
- raw read BAM: `01_depth/<CELL_LINE>.bam`
- normalized depth
- reference/telomere/repeat/censat annotation

핵심 동작:

- contig 전처리 후 extended PAF 형태의 `.ppc.paf` 생성
- telomere 연결 contig 탐색
- NClose 계산 및 compress 정보 생성
- 필요 시 내부적으로 `03_Anal_bam.py`를 호출해 read-level evidence로 translocation NClose 방향을 보정
- breakend graph 구성
- telomere-to-telomere path 탐색
- 경로 수가 너무 많으면 제한값을 낮춰 재시도
- inversion-only 그래프에서 짧은 cycle을 찾아 ecDNA 후보를 저장

주요 산출물:

- `<primary paf>.ppc.paf`
- `<PREFIX>/telomere_connected_list.txt`
- `<PREFIX>/telomere_connected_list_readable.txt`
- `<PREFIX>/00_raw/<chrpair>/<n>.index.txt`
- `<PREFIX>/path_data.pkl`
- `<PREFIX>/path_di_data.pkl`
- `<PREFIX>/paf_file_path.pkl`
- `<PREFIX>/report.txt`
- `<PREFIX>/nclose_chunk_data.pkl`
- `<PREFIX>/ecdna_circuit_data.pkl`
- `<PREFIX>/conjoined_type4_ins_del.pkl`

주의:

- `run_pipe.sh`에서 `--verbose`를 항상 넘기므로 `00_raw` 인덱스 파일 생성이 실제 파이프라인의 일부다.
- 이 단계가 실패하면 대부분의 후속 단계는 진행 불가다.

## Breakend Graph 구현 메모

사용자가 정리한 논문 초안을 코드에 대응시키면, `02_Build_Breakend_Graph_Limited.py`의 핵심은 "breakend abstraction + read evidence 보정 + reverse-complement를 반영한 이중 그래프"로 요약된다.

### 1. Abstraction

코드는 PAF의 각 alignment chunk를 하나의 node로 본다.

- 하나의 contig/unitig는 chunk의 순서열로 표현된다.
- preprocessing 뒤 `extract_bnd_contig()`가 breakend를 포함하는 contig를 우선 추린다.
- `extract_nclose_node()`가 각 breakend contig에서 nclose pair를 만든다.
- `nclose_calc()`는 이 pair들을 후처리하면서 압축과 필터링을 수행한다.

구현상 중요한 점:

- nclose pair는 `(s, e)`로 저장되며 순서가 의미를 가진다.
- `get_corr_dir()`는 pair가 압축 과정에서 반대 순서로 관측되더라도 strand를 기준 방향으로 다시 정규화한다.
- `nclose_start_compress`, `nclose_end_compress`는 시작점과 끝점이 가까운 pair를 묶기 위한 압축 정보다.
- 압축 전/후 pair 목록은 각각 `all_nclose_nodes_list.txt`, `compressed_nclose_nodes_list.txt`로 남는다.

### 2. BAM analysis가 그래프에 들어오는 방식

논문 초안의 BAM analysis는 구현에서 `02` 단계 내부 substep으로 들어간다.

흐름:

- `nclose_calc()`가 translocation 성격의 nclose만 추려 `03_anal_bam_input.pkl`을 만든다.
- 그 다음 `03_Anal_bam.py`가 read BAM을 직접 읽어 orientation별 지지 read 수를 센다.
- read evidence가 없고 depth도 깨끗하면 기존 nclose를 삭제할 수 있다.
- 반대 방향 orientation이 충분히 지지되면 reverse-direction nclose를 추가할 수 있다.
- reciprocal translocation으로 보이면 `nclose2cov.pkl`에 coverage를 기록한다.

즉, 최종 breakend graph는 preprocessing 결과만으로 고정되지 않고, read-level evidence를 통해 방향성이 보강된 nclose 집합을 입력으로 삼는다.

### 3. Reverse-Complement를 반영한 이중 그래프

핵심 함수는 `initialize_bnd_graph(contig_data, nclose_nodes, telo_contig)`이다.

이 구현은 생물학적 node 하나를 그대로 vertex 하나로 쓰지 않고, "그 node를 어느 방향으로 읽고 있는가"를 분리해서 이중 그래프로 만든다.

코드 상수:

- `DIR_FOR = 1`, `DIR_BAK = 0`
- `DIR_OUT = 2`, `DIR_IN = 3`

의미:

- nclose node `x`는 `(DIR_FOR, x)`와 `(DIR_BAK, x)` 두 상태를 가진다.
- 이는 같은 chunk라도 unitig를 forward로 읽는 경우와 reverse-complement 방향으로 읽는 경우를 분리한 것이다.
- telomere에 연결된 chunk는 legacy 구현에서 `(DIR_OUT, x)`, `(DIR_IN, x)`로 표현된다.
- telomere endpoint 자체는 layer를 따로 갖지 않는다. base graph에서는 `chr1f`, `chr1b` 같은 문자열 노드이고, path graph에서는 `(chr1f, chr_change, dir_change)`처럼 확장된다.
- `telomere_connected_list.txt`에는 `(layer, node_idx, dist)`가 저장되지만, current legacy graph에서는 이 layer가 terminal state로 이중화되어 유지되지는 않는다.

`chr_correlation_maker()`의 역할:

- `chr1f`, `chr1b` 같은 telomere endpoint 이름과 내부 synthetic index를 매핑한다.
- 여기서 `f`/`b`는 chromosome front/back telomere endpoint를 뜻한다.
- path 탐색 시 source/target은 이 telomere endpoint 문자열을 사용한다.

논문 표기와 코드의 대응:

- 논문의 `(x, 0)`, `(x, 1)`은 코드에서 대체로 `(DIR_BAK, x)`, `(DIR_FOR, x)`에 대응한다.
- telomere node `t`도 실제 구현에서는 "외부 endpoint 문자열 노드"와 "그 endpoint에 연결된 실제 chunk node의 in/out 상태"로 나뉜다.
- 다만 terminal layer는 `DIR_OUT`/`DIR_IN`으로만 표현되므로, 내부 `DIR_FOR`/`DIR_BAK`처럼 완전히 이중화되어 있지는 않다.

### 4. NClose 내부 edge

논문의 규칙 a)는 코드에서 아래 두 줄로 구현된다.

- `(DIR_FOR, s) -> (DIR_FOR, e)`
- `(DIR_BAK, e) -> (DIR_BAK, s)`

해석:

- forward 복사본에서는 pair의 시작점에서 끝점으로 진행한다.
- reverse-complement 복사본에서는 끝점에서 시작점으로 되돌아간다.
- 즉, 하나의 nclose pair를 넣을 때 항상 그 reverse-complement 짝 edge도 함께 들어간다.

### 5. NClose 사이 edge

논문의 규칙 b)는 코드에서 `i_ind`, `j_ind`의 8개 경우 분기로 구현된다.

```text
f+ / b+ / f- / b-
```

여기서:

- `f`와 `b`는 해당 node가 nclose pair의 첫 번째 끝인지 두 번째 끝인지
- `+`와 `-`는 chunk의 `CTG_DIR`

이 분기들은 논문에서 정의한 다음 조건을 코드 형태로 푼 것이다.

- 두 node가 같은 reference chromosome에 있어야 함
- 두 node에서 기대하는 reference position 변화 방향이 같아야 함
- 실제 reference 좌표의 선후 관계가 그 기대 방향과 모순되지 않아야 함

구현 디테일:

- 코드는 `i_ind/j_ind` 조합마다 `CHR_STR`, `CHR_END` 비교를 직접 적어서 조건을 푼다.
- 논문의 `incr`와 `sgn`를 함수형으로 계산하지 않고, 각 경우의 수를 수동 분기한 형태다.
- edge를 하나 추가할 때 항상 reverse-complement 짝 edge도 같이 추가한다.

예를 들어, 어떤 경우에는

- `(DIR_BAK, x) -> (DIR_FOR, y)`
- `(DIR_BAK, y) -> (DIR_FOR, x)`

처럼 들어가고, 다른 경우에는

- `(DIR_FOR, x) -> (DIR_FOR, y)`
- `(DIR_BAK, y) -> (DIR_BAK, x)`

처럼 들어간다.

첫 번째 예시는 cross-layer edge다.

- 대표적으로 `f+` 와 `f-`를 잇는 경우
- 대표적으로 `b+` 와 `b-`를 잇는 경우
- 즉, 두 endpoint가 같은 chromosome에 있더라도 strand와 endpoint 타입까지 같이 보면 같은 layer에서 바로 잇는 것은 모순이고, 한쪽은 reverse-complement layer에서 나가서 다른 쪽의 forward layer로 들어가야 하는 경우다

코드 관점에서는 "reference 좌표의 증가/감소 방향을 맞추려면 layer를 바꿔야 하는 경우"라고 이해하면 된다.

두 번째 예시는 same-layer edge다.

- 같은 nclose pair 내부에서는 항상 `(DIR_FOR, s) -> (DIR_FOR, e)`와 `(DIR_BAK, e) -> (DIR_BAK, s)`를 넣는다
- 서로 다른 nclose pair 사이라도 `b+` 와 `f+`
- 또는 `f-` 와 `b-`
- 처럼 strand와 endpoint 타입 조합상 같은 reading direction을 유지한 채 연결할 수 있으면 same-layer edge가 들어간다

즉:

- cross-layer edge는 "이 연결을 따라가면 읽기 방향 복사본을 한 번 바꿔야 한다"는 뜻
- same-layer edge는 "이 연결을 따라가도 현재 읽기 방향 복사본을 유지할 수 있다"는 뜻

그리고 왜 항상 두 줄씩 같이 들어가느냐가 중요하다.

- 하나는 현재 방향에서 본 연결
- 다른 하나는 그 연결을 reverse-complement로 읽었을 때의 대응 연결

이 짝 추가가 바로 구현상의 "reverse-comp.를 고려한 이중 그래프"다.

### 6. Telomere 관련 edge

논문의 규칙 c), d)는 구현에서 두 단계로 나뉜다.

첫 번째:

- `edge_optimization()`이 telomere endpoint마다 실제로 붙일 contig chunk를 고른다.
- 결과는 `telomere_connected_list.txt`에 저장된다.

두 번째:

- `initialize_bnd_graph()`가 이 정보를 사용해 telomere-to-nclose, telomere-to-telomere edge를 만든다.

구현 포인트:

- `edge_optimization()`이 저장하는 `(layer, node_idx, dist)`의 `layer`는 telomere-connected 후보를 고를 때 사용된다.
- 하지만 `initialize_bnd_graph()`는 최종 breakend graph를 만들 때 이 layer를 terminal state로 분리해 보존하지 않고, generic `(DIR_OUT, idx)`, `(DIR_IN, idx)`를 쓴다.
- telomere-connected chunk와 nclose node의 chromosome이 같아야 한다.
- telomere chunk와 nclose chunk의 상대 위치가 단조 증가인지 감소인지 `dir1/dir2`로 본다.
- `curr_contig[CTG_DIR]`와 telomere label의 `f/b`를 합쳐 `t_ind`를 만들고, `k_ind`와 조합해서 어떤 복사본끼리 이어야 하는지 결정한다.
- 현재 로컬 수정본은 여기에 `is_telo_dir_consistent()`를 추가해, telomere suffix가 요구하는 접근 방향과 맞는 경우만 telomere-to-nclose edge를 허용한다.
- 염색체 양 끝 telomere끼리는 결국 `(DIR_OUT, node_a) -> (DIR_IN, node_b)` 꼴로 연결된다.

즉, telomere는 "외부 endpoint 문자열 + terminal adapter 상태"로 모델링되지만, terminal 자체는 legacy `DIR_OUT/DIR_IN` 구조를 유지한다.

### 7. Path 탐색용 상태 확장

`initialize_bnd_graph()`가 만든 이중 그래프는 바로 탐색하지 않고, `make_graph()`에서 다시 상태 공간으로 확장된다.

확장 뒤 node:

- telomere endpoint는 `(chr1f, chr_change, dir_change)` 같은 상태
- 내부 breakend node는 `((DIR_FOR|DIR_BAK|DIR_OUT|DIR_IN), idx, chr_change, dir_change)` 형태

의미:

- `chr_change`는 inter-chromosomal 이동 횟수
- `dir_change`는 same-contig inversion 성격의 방향 전환 횟수

따라서 실제 탐색 그래프는 "내부 breakend에 대해서는 reverse-complement 이중 그래프" 위에 "재배열 복잡도 제한 상태"를 한 번 더 덧씌운 구조라고 보는 것이 맞다. 다만 terminal은 legacy adapter state라 내부 node만큼 세분화되어 있지는 않다.

### 8. `00_raw/*.paf` 가상 염색체 파일 읽는 법

`02_Build_Breakend_Graph_Limited.py`는 `--verbose`일 때 path 하나마다

- `00_raw/<chrpair>/<n>.paf`
- `00_raw/<chrpair>/<n>.index.txt`

를 함께 저장한다.

해석 원칙:

- `*.paf`는 사람이 읽기 쉬운 chunk 나열이다.
- `*.index.txt`는 같은 path를 그래프 node 상태로 기록한 파일이다.
- 방향 해석은 `*.paf`만 보면 부족할 수 있으므로, reverse-complement 여부는 반드시 대응되는 `*.index.txt`와 같이 봐야 한다.

#### Telomere node 해석

telomere endpoint 이름의 `f`/`b`는 contig의 `+` 방향 기준으로 붙는다.

- contig의 앞쪽이 telomere 쪽이면 `f`
- contig의 뒤쪽이 telomere 쪽이면 `b`

따라서 contig 자체가 `-` strand라면 실제 reference에서 보이는 앞/뒤 해석은 반대로 뒤집힌다.

즉, `chr2b` 같은 표기는 "2번 염색체의 back telomere endpoint에 연결된 상태"를 뜻하지만, 실제 chunk가 `+`인지 `-`인지에 따라 path를 읽을 때의 local orientation은 달라질 수 있다.

#### Reverse-complement 연결 해석

가상 염색체 파일에서 같은 contig의 chunk 순서가 뒤집혀 보이면, 그 구간은 reverse-complement로 연결된 것으로 해석해야 한다.

- 이 경우 그래프 상태로는 `DIR_BAK`
- 즉, 해당 contig를 backward copy에서 읽고 있다는 뜻이다

반대로 같은 contig의 chunk가 원래 contig 좌표 순서대로 진행하면 `DIR_FOR`로 읽고 있는 것이다.

중요한 점:

- `*.paf`의 줄 순서만으로는 항상 확정할 수 없다
- 최종 판단은 대응되는 `*.index.txt`의 첫 번째 컬럼을 봐야 한다

`*.index.txt`의 첫 번째 값 해석:

- `1` = `DIR_FOR`
- `0` = `DIR_BAK`
- `2` = `DIR_OUT`
- `3` = `DIR_IN`

따라서 path 내부 contig chunk가 `0`으로 기록되어 있으면 그 chunk는 reverse-complement layer에서 읽힌 것이다.

다만 start/end chunk는 generic terminal out/in state로 나타나므로, endpoint 방향은 `CTG_TELCON`과 같이 해석해야 한다.

- `DIR_OUT`는 telomere에서 시작해 첫 chunk 쪽으로 나가는 terminal adapter 상태다
- `DIR_IN`은 마지막 chunk에서 telomere로 닫히는 terminal adapter 상태다
- terminal state 자체에는 `FOR/BAK`가 없으므로, endpoint orientation은 반드시 `CTG_TELCON`, `CTG_DIR`, 그리고 path 전파 결과를 같이 봐야 한다.

telomere 위치 판정 규칙:

- `CTG_DIR='+'`일 때 `CTG_TELCON=*f`면 contig front, `*b`면 contig back에 telomere가 있다
- `CTG_DIR='-'`일 때는 반대로 `*f`면 contig back, `*b`면 contig front에 telomere가 있다

만약 endpoint chunk의 `CTG_TELCON=0`이면, 그 chunk의 telomere label은 `*.index.txt`의 시작/끝 tuple line에 적힌 telomere 이름을 implied label로 사용한다.

이 규칙으로 먼저 start chunk의 local orientation을 정하고, 그 다음 middle chunk는 `DIR_FOR/DIR_BAK`로 전파해서 읽어야 한다.

가상 염색체 consistency 체크 순서:

1. 시작 telomere와 첫 chunk의 `CTG_TELCON`이 맞는지 본다.
2. 첫 chunk의 `CTG_TELCON`과 `CTG_DIR`로 telomere가 contig front/back 중 어디에 있는지 정한다.
3. 시작 state가 `DIR_OUT`이면 telomere 쪽 끝에서 내부로 읽으므로, 현재 가상 염색체에서 첫 chunk가 `+`로 놓이는지 `-`로 놓이는지 결정한다.
4. middle chunk는 `DIR_FOR`면 현재 orientation 유지, `DIR_BAK`면 reverse-comp로 읽는다고 해석한다.
5. 마지막 chunk의 `CTG_TELCON=0`이면 tuple line의 terminal telomere 이름으로 보정한다.
6. 마지막 chunk는 `DIR_IN`과 그 telomere label을 같이 봐서, 현재 orientation에서 telomere가 끝점에 와야 한다.
7. 마지막 telomere가 현재 orientation과 맞지 않으면 그 path는 모순이다.

#### 정리: experimental terminal 이중화는 `legacy6`에 백업됨

이 대화 중 한때 terminal state까지 `FOR/BAK`로 이중화하는 실험을 했고, 그 결과를 기준으로 설명한 문장이 남아 있었다.

- 그 실험본은 현재 mainline이 아니다.
- 관련 Python 백업은 `/home/saycorn/00_build_graph/legacy6`에 보관했다.
- 현재 mainline `02_Build_Breakend_Graph_Limited.py`는 legacy terminal model, 즉 `DIR_OUT=2`, `DIR_IN=3`을 쓴다.

따라서 future agent는 현재 repo의 `*.index.txt`를 해석할 때 `2/3/4/5` terminal model을 가정하면 안 된다. start/end는 `2/3`만 쓰고, tuple line도 `(chr_label, chr_change, dir_change)` 형식으로 읽는 것이 현재 코드와 맞다.

#### 예시: `U2OS_20_58_22/00_raw/chr2b_chr15b/1`

원본 파일 경로:

- `/home/hyunwoo/60g_skype/30_skype_pipe/U2OS_20_58_22/00_raw/chr2b_chr15b/1.paf`
- `/home/hyunwoo/60g_skype/30_skype_pipe/U2OS_20_58_22/00_raw/chr2b_chr15b/1.index.txt`

외부 파일이 없어져도 같은 예시를 재현할 수 있도록 내용을 여기에 그대로 적는다.

`1.paf`

```text
('chr2b', 0, 0)
ptg000473l	27864	0	9209	+	chr2	242696752	135664279	135673488	60	3	2247	2247	0	0	chr2b	0	0	0	+	chr2	0.5936
utg058613l	31782	0	11797	-	chr2	242696752	31270	43066	60	1	5820	5821	0	0	0	0	0	0	-	chr15	1.70292
utg058613l	31782	11797	31782	-	chr15	99753195	17445128	17465114	0	1	5820	5821	0	0	0	chr15	rin	rin	-	chr15	1.70293
utg075486l	13248	931	13248	-	chr15	99753195	13034119	13046488	1	3	6228	6228	0	0	chr15b	chr15	rin	rin	-	chr15	1.89366
('chr15b', 1, 0)
[('chr2', 135642218), ('chr15', 4430995)]
```

`1.index.txt`

```text
('chr2b', 0, 0)
2	2247	0	0
1	5820	0	0
1	5821	1	0
3	6228	1	0
('chr15b', 1, 0)
[('chr2', 135642218), ('chr15', 4430995)]
```

이 path의 index는 아래 상태를 뜻한다.

- `('chr2b', 0, 0)`에서 시작
- `(legacy DIR_OUT=2, 2247, 0, 0)`
- `(DIR_FOR, 5820, 0, 0)`
- `(DIR_FOR, 5821, 1, 0)`
- `(legacy DIR_IN=3, 6228, 1, 0)`
- `('chr15b', 1, 0)`에서 종료

해석:

- 시작 telomere는 `chr2b`다
- 첫 chunk `2247`의 `CTG_TELCON=chr2b`, `CTG_DIR='+'`이므로 telomere는 이 contig의 back 쪽에 있다
- 가상 염색체는 telomere에서 시작해야 하므로, 현재 path에서 첫 chunk의 local orientation은 `2-`로 읽어야 한다
- 이후 middle state는 `DIR_FOR -> DIR_FOR`이므로 현재 orientation을 유지한 채 `2- => 15-`로 전파된다
- 마지막 chunk `6228`은 `CTG_TELCON=chr15b`, `CTG_DIR='-'`이므로 이 경우 telomere는 contig front 쪽에 있다
- 그런데 현재 orientation은 `15-`로 전파된 상태이므로, `15b` telomere가 가상 염색체의 끝점에 와야 하는 조건과 충돌한다

즉, 이 예시는:

- middle chunk만 보면 `DIR_BAK`가 없어서 단순해 보이지만
- endpoint orientation까지 같이 보면 start에서 정해진 current orientation과 마지막 `chr15b`가 모순이다

따라서 이 path는 "일관된 telomere-to-telomere virtual chromosome" 예시가 아니라, telomere endpoint 해석 규칙을 제대로 적용하면 모순을 검출해야 하는 반례로 보는 것이 맞다.

위 예시는 현재 mainline과 같은 legacy 출력 형식이다. start/end는 `2/3`을 쓰고, tuple line도 `(chr_label, chr_change, dir_change)` 형식이다.

다만 이 legacy 구조에서는 terminal layer가 flatten되어 있으므로, `chr2b_chr15b/1` 같은 explicit telomere anchoring 반례는 여전히 주의해서 읽어야 한다.

### 9. 코드 독해 시 바로 보면 좋은 함수

breakend 구현을 따라갈 때는 아래 순서가 가장 효율적이다.

1. `nclose_calc()`
2. `03_Anal_bam.py`
3. `initialize_bnd_graph()`
4. `make_graph()`
5. `run_graph_pipeline()`

핵심 해석 규칙:

- `nclose_nodes`는 그래프의 breakend anchor 집합이다.
- `DIR_FOR/DIR_BAK`는 reverse-complement를 반영한 doubled node state다.
- telomere는 `chr*f|chr*b` 외부 노드와 `(DIR_OUT|DIR_IN, idx)` 상태 노드가 함께 써진다.
- path는 "telomere endpoint -> breakend chain -> opposite telomere endpoint"를 찾는 과정이다.
- 최종적으로 path 하나는 breakend를 포함한 mutant chromosome 후보 하나로 해석된다.

### 10. U2OS 런타임 메모

`bash run_pipe.sh U2OS`가 3분을 넘는다고 바로 그래프 loop 문제로 단정하면 안 된다.

- `02_Build_Breakend_Graph_Limited.py`는 내부에서 `03_Anal_bam.py`를 호출한다.
- 실제 U2OS 재검증에서 `02` 단독 실행은 `22:39:22` 시작, `22:49:04` 성공 로그 종료로 약 9분 42초가 걸렸다.
- 이 시간의 상당 부분은 BAM analysis가 사용했다.

즉, U2OS에서 3분 초과는 "이상 징후일 수도 있지만, 그것만으로는 그래프 무한루프의 증거가 아니다"라고 문서화해 두는 것이 맞다.

### 2. `11_Ref_Outlier_Contig_Modify.py`

이 스크립트는 이름과 달리 `.ppc.paf` 자체를 다시 쓰지 않는다.

동작:

- `02`에서 만든 `.ppc.paf`와 원본 PAF를 읽는다
- `<PREFIX>/conjoined_type4_ins_del.pkl`을 참조한다
- type4 및 conjoined type2 이벤트를 depth 측정용 PAF로 분리한다

출력:

- `<PREFIX>/11_ref_ratio_outliers/front_jump/*.paf`
- `<PREFIX>/11_ref_ratio_outliers/back_jump/*.paf`
- 각 이벤트에 대응하는 `*_base.paf`

### 3. `21_run_depth.py`

동작:

- `path_data.pkl`과 `00_raw/*.index.txt`를 이용해 path별 최종 PAF를 분리 생성
- ecDNA/type2 관련 outlier PAF도 생성
- `pandepth`를 path별 PAF와 outlier PAF에 대해 각각 실행
- path 깊이 벡터 인덱스를 정리해 다음 단계 입력으로 저장

주요 산출물:

- `<PREFIX>/21_pat_depth/*.paf`
- `<PREFIX>/21_pat_depth/*.win.stat.gz`
- `<PREFIX>/11_ref_ratio_outliers/ecdna/*.paf`
- `<PREFIX>/11_ref_ratio_outliers/type2_ins/*.paf`
- outlier 폴더 아래 `*.win.stat.gz`
- `<PREFIX>/contig_pat_vec_data.pkl`

외부 의존성:

- `./PanDepth/bin/pandepth` 바이너리

### 4. `22_save_matrix.py`

동작:

- 샘플 normalized depth와 path별 depth를 같은 좌표계의 벡터로 변환
- censat/high-depth 구간을 별도로 분리
- path, outlier, ecDNA를 모두 열(column)로 하는 설계 행렬을 생성
- 기본 chromosome path를 `tar_chr_data.pkl`로 저장

출력:

- `<PREFIX>/matrix.h5`
- `<PREFIX>/23_input.pkl`
- `<PREFIX>/tar_chr_data.pkl`
- `<PREFIX>/B.npy`

주의:

- 기본 파이프라인은 `--not_use_nclose_weight`를 사용하므로 nclose weight는 꺼진 상태다.

### 5. `23_run_nnls.py`

동작:

- `matrix.h5`에서 `A`, `B`, `A_fail`을 읽는다
- `skglm`의 `GeneralizedLinearEstimator + AndersonCD + PositiveConstraint`로 NNLS를 푼다
- nclose를 가설검정 기반으로 부분 필터링하면서 다시 적합한다

출력:

- `<PREFIX>/weight.npy`
- `<PREFIX>/predict_B.npy`
- `<PREFIX>/A_idx_list.pkl`

외부 의존성:

- `skglm`
- `h5py`

### 6. `24_cluster_weight.py`

동작:

- `23`의 결과를 바탕으로 Julia solver로 추가 필터링과 path cluster 재선정을 수행한다

실행 조건:

- `report.txt`의 총 path 수가 `HARD_PATH_COUNT_BASELINE` 이하일 때만 의미 있게 실행된다
- 조건을 넘으면 초반에 `exit(0)` 하므로 filter/cluster 산출물이 생기지 않는다

출력:

- `<PREFIX>/weight_filter.npy`
- `<PREFIX>/predict_B_filter.npy`
- `<PREFIX>/weight_cluster.npy`
- `<PREFIX>/predict_B_cluster.npy`

외부 의존성:

- `juliacall`
- `anderson_nnls.jl`

### 7. `30_virtual_sky.py`

동작:

- `weight.npy` 기반으로 기본 Virtual SKY를 그림
- `24`가 유효하게 돌 수 있는 경우 `_filter`, `_cluster` 버전도 추가 생성
- path 위치 목록을 `tot_loc_list.pkl`로 저장해 `31`이 재사용한다

출력:

- `<PREFIX>/virtual_sky.pdf`
- `<PREFIX>/virtual_sky.png`
- 조건부: `virtual_sky_filter.*`, `virtual_sky_cluster.*`
- `<PREFIX>/tot_loc_list.pkl`

### 8. `31_depth_analysis.py`

동작:

- Circos plot 생성
- breakend / insertion / deletion / amplicon을 BED로 출력
- VCF 스타일 SV 결과를 출력
- 기본 결과 외에 `_filter`, `_cluster` 버전도 조건부 생성

출력:

- `<PREFIX>/total_cov.pdf`
- `<PREFIX>/total_cov.png`
- 조건부: `total_cov_filter.*`, `total_cov_cluster.*`
- `<PREFIX>/SKYPE_result.bed`
- `<PREFIX>/SV_call_result.vcf`
- 조건부: `_filter`, `_cluster` 변형
- `<PREFIX>/cn_data.pkl`

마지막 정리:

- 스크립트 끝에서 `<PREFIX>/matrix.h5`를 삭제한다

## 실제 의존 관계

가장 중요한 파일 흐름만 적으면 아래와 같다.

```text
run_pipe.sh
  -> 00_depth_norm.py
     -> <CELL_LINE>_normalized.win.stat.gz
  -> 02_Build_Breakend_Graph_Limited.py
     -> .ppc.paf
     -> path_data.pkl / path_di_data.pkl / nclose_chunk_data.pkl / report.txt / 00_raw/*
     -> conjoined_type4_ins_del.pkl / ecdna_circuit_data.pkl
  -> 11_Ref_Outlier_Contig_Modify.py
     -> 11_ref_ratio_outliers/front_jump/*, back_jump/*
  -> 21_run_depth.py
     -> 21_pat_depth/*.win.stat.gz
     -> contig_pat_vec_data.pkl
  -> 22_save_matrix.py
     -> matrix.h5 / 23_input.pkl / tar_chr_data.pkl / B.npy
  -> 23_run_nnls.py
     -> weight.npy / predict_B.npy
  -> 24_cluster_weight.py
     -> weight_filter.npy / weight_cluster.npy
  -> 30_virtual_sky.py
     -> virtual_sky*.png/pdf
  -> 31_depth_analysis.py
     -> total_cov*.png/pdf / SKYPE_result*.bed / SV_call_result*.vcf
```

## 에이전트 작업 시 우선순위

파이프라인 동작을 바꾸는 수정은 보통 아래 순서로 보면 된다.

1. `run_pipe.sh`
2. `02_Build_Breakend_Graph_Limited.py`
3. `21_run_depth.py`
4. `22_save_matrix.py`
5. `23_run_nnls.py`
6. `24_cluster_weight.py`
7. `30_virtual_sky.py`
8. `31_depth_analysis.py`

실무적으로 중요한 점:

- 가장 큰 fan-out은 `02_Build_Breakend_Graph_Limited.py`다.
- `11_Ref_Outlier_Contig_Modify.py`는 보조 이벤트 PAF를 만드는 단계이지 `.ppc.paf` 수정 단계가 아니다.
- `30_virtual_sky.py`와 `31_depth_analysis.py`는 `24`의 실행 조건에 따라 `_filter`/`_cluster` 산출물 존재 여부가 달라진다.
- 최종 검증 시 `31_depth_analysis.py`가 `matrix.h5`를 삭제한다는 점을 염두에 둬야 한다.
