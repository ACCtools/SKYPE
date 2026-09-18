# Outputs and weight interpretation

[README](../README.md) · [Usage](usage.md) · [Methods](methods.md)

Unless a section says otherwise, this document describes native assembly
output. [VCF input](#vcf-input-outputs) and full-assembly input retain their
separate reporting rules.

## Choosing a result file

| File | Unit of a row or record | Purpose |
| --- | --- | --- |
| `SV_call_result.vcf` | A BND mate record or ordinary symbolic DEL/DUP | Variant export. A BND adjacency has two reciprocal records. |
| `SV_call_result.bnd_weights.tsv` | One structure's contribution to an emitted BND | Reproduce BND weights after PATH_SPLIT and exact geometry merging. |
| `SKYPE_result.bed` | A NClose endpoint/span or a structure summary span | Inspect genomic intervals with an explicit weight scope. |
| `nclose_report.tsv` | One original Type 1/2/4 NClose | Inspect the full NClose total and preprocessing history. |
| `structure_report.tsv` | One structure | Inspect its own coefficient and source. |
| `structure_nclose_usage.tsv` | One structure–original NClose association | Trace occurrence counts and contributions, including zero-weight structures. |
| `total_cov.png`, `total_cov.pdf` | Whole-genome plot | Compare observed/fitted depth and inspect displayed junctions/structure summaries. |
| `cn_data.pkl` | Four lists of normalized support values | Analyze CN distributions using the units described below. |

Reports include zero-weight entries. Native variant displays apply the strict
`CN > 0.1` threshold after aggregation; compound summary displays apply it to
the structure's own CN. Positive support in a report does not by itself mean
that an event passes the display threshold or establishes biological accuracy.

## Weight accounting

Let `w[s]` be a structure's original weight, `count[s, j]` its use count for
original NClose `j`, and `N` the reporting normalization unit:

```text
N = median read depth excluding chrM / 2
structure_CN[s] = w[s] / N
NClose_CN[j] = sum(w[s] × count[s, j] for all structures s) / N
```

Paths, independent Type 4, MERGE_TYPE4 and AMP use their own fitted matrix-column
coefficients. VIRTUAL_INV uses its qualified raw-read-derived expected depth
as its own weight. Centromere-only features have no NClose membership.

Different structures add. Within a path, different traversal positions add,
while duplicate descriptions of the same occurrence at one position are counted
once. An already aggregated NClose value is never used as a new structure
coefficient. Accounting preserves original NClose IDs; it does not introduce
coordinate clustering.

### Worked example

The weights below are already normalized by `N`:

| Structure | Own CN | Original NClose uses |
| --- | ---:| --- |
| Path P | 3 | A once; Type 4 event T once |
| Independent Type 4 feature | 1 | T once |
| AMP M | 2 | A once; B once |
| MERGE_TYPE4 G | 1 | B once; C once |
| VIRTUAL_INV V | 0.5 | A once; B once |

The NClose totals are A = `5.5`, B = `3.5`, C = `1`, and T = `4`.
T receives `3` from the path and `1` from its independent feature. Its total
`4` is not added again to the path's contribution. A count of two in one
structure would contribute twice that structure's weight.

### Virtual-inversion contributions

The two constituent NCloses each receive the virtual structure's weight times
their use count. The weight is not divided by two. If a constituent also occurs
in a path, AMP or another structure, those contributions add to its virtual
contribution. Two chromosome-arm display intervals for one virtual pair are
views of one structure, not two independent weights.

Current qualification requires a depth-balanced record, NClose VAF greater than
`0.2` on both sides, and the existing point-spanning evidence or qualifying
same-chromosome contiguous-span path. The spanning-path support gate remains
separate from the final output-weight gate. Virtual estimates are recorded
separately from fitted-model contributions and do not change `weight.npy` or
`predict_B.npy`.

### AMP and reference depth

Structure CN, NClose CN and segment depth are different quantities. Suppose an
AMP has structure CN `2`, uses A and B once each, and its depth feature traverses
a reference interval twice. Its contribution is `2` to A and `2` to B, while
its contribution to that interval's depth is `4 × N` in raw-depth units.

The AMP summary reports structure CN `2`; there is no additional ×2 display
normalization. Multiple reference traversals are already represented in its
depth vector. These CN/support values are not read VAFs or independently
validated absolute copy numbers.

## Output representations

### Original NCloses and compound summaries

AMP, MERGE_TYPE4 and VIRTUAL_INV export their original constituent NClose BNDs
in the native VCF. A depth-model bridge between unitigs is not exported as an
additional inferred junction. Ordinary unmerged Type 4 events remain symbolic
DEL/DUP calls.

BED and Circos also retain the compound structure's summary span. That summary
uses the structure's own weight; its constituent NClose views use their summed
support. A summary and its constituent views are different descriptions of
related evidence and must not be summed as independent events.

### PATH_SPLIT

The existing path-specific handling of internal chromosome transitions can
represent a parent NClose with several BNDs. The shared model records which
path occurrences produce those child projections. Only those occurrences move
from the parent display to the corresponding children; any other structure's
contribution to the parent remains there.

This transfer applies to VCF, BED and Circos. The original NClose's full sum
remains in `nclose_report.tsv` and the original-NClose CN distributions.
PATH_SPLIT does not add new original NClose IDs or implement a general
primitive-junction decomposition of every NClose.

Membership is recorded regardless of each path's coefficient. For example,
two paths with CN `0.06` that produce the same split can jointly pass the
`> 0.1` gate. A parent can disappear from the displayed call set after its
contributions move to children, while its original NClose report remains positive.

### Exact VCF geometry merging

VCF output merges exactly matching endpoint/retained-side pairs and aggregates
their contributions before its output threshold. Distinct original NCloses
remain distinct in the accounting and original-NClose reports. Their IDs are
retained together on the merged VCF call. Consequently, VCF junction counts
and weights need not match the number of rows or one individual row's weight
in the original-NClose report.

## Native VCF fields

Each BND adjacency has two reciprocal `MATEID` records with the same support.
Do not count the two mate records as two independent junctions. BND coordinates
convert native 0-based boundaries to 1-based retained anchor bases: a boundary
at `p` uses anchor `p` for the retained left side and `p + 1` for the retained
right side. Records are sorted by reference contig order and position.

| INFO field | Meaning |
| --- | --- |
| `SVTYPE`, `END`, `SVLEN` | Variant type and reference span information where applicable. |
| `WEIGHT` | Final normalized support. For BNDs, equals `MODEL_WEIGHT + VIRTUAL_WEIGHT`. |
| `WEIGHT_METHOD` | `STRUCTURE_SUM` for current native calls. |
| `NCLOSE_IDS` | Original IDs joining to `nclose_report.tsv`; available on BND and ordinary Type 4 records. |
| `NCLOSE_KEYS` | BND source node-index pairs, encoded as `first:second`. |
| `PARENT_IDS` | Contributing structure IDs joining to `structure_report.tsv`. |
| `PARENT_WEIGHTS` | Each structure's own normalized weight, in `PARENT_IDS` order. |
| `PARENT_MULTIPLICITY` | Occurrence counts for this exact adjacency within each parent, in the same order. |
| `MODEL_WEIGHT` | Sum of fitted-column occurrence contributions. |
| `VIRTUAL_WEIGHT` | Sum of virtual-structure occurrence contributions. |
| `MODEL_FEATURE_COUNT` | Number of distinct positive fitted columns contributing to the BND. |
| `MODEL_OCCURRENCE_COUNT` | Contributing occurrence count over those fitted columns. |
| `SVCLASS` | BND source classes, such as `NCLOSE`, `AMPLICON`, `MERGED_TYPE4`, `VIRTUAL_INV`, `PATH_SPLIT`. |
| `CTG_NAME` | Supporting contig names or the ordinary Type 4 source label. |
| `STRANDS`, `MATEID` | BND orientation and reciprocal mate linkage. |
| `BP_STEP_DEPTH_RATIO_B` | Signed observed-depth step divided by exported support in raw-depth units. |
| `BP_STEP_DEPTH_RATIO_PREDICT_B` | Corresponding signed step in fitted depth. |

The structure/occurrence fields describe BND contributions; ordinary symbolic
Type 4 records carry their original NClose ID and summed support without that
full BND field set. Depth-ratio values are ordered local/mate for BNDs and
POS/END for symbolic records. A dot denotes an unavailable ratio. These depth
diagnostics do not by themselves establish junction accuracy.

`SV_call_result.bnd_weights.tsv` contains `bnd_id`, `feature_index`,
`source_occurrence`, `occurrence_count`, `feature_weight_N`, `contribution_N`,
`weight_method`, `structure_id`, `nclose_id` and `structure_kind`.
`source_occurrence` identifies the original NClose display or its PATH_SPLIT
projection. Summing `contribution_N` by `bnd_id` reproduces each BND's `WEIGHT`.
Virtual rows have `feature_index=.`.

## TSV reports and IDs

### Structure report

`structure_report.tsv` has one row per structure, including zero-weight
structures:

| Column | Meaning |
| --- | --- |
| `structure_id` | Unique structure identifier in this result/model. |
| `kind` | `PATH`, `TYPE4`, `MERGE_TYPE4`, `AMP`, `VIRTUAL_INV` or `CENTROMERE`. |
| `feature_index` | Zero-based matrix-column index; `.` for a virtual structure. |
| `raw_weight`, `weight_N` | Own weight in raw-depth and normalized units. |
| `source` | Original feature/path location or raw virtual-pair reference; a logical path location need not be a materialized PAF file. |
| `legacy_ids` | Existing compound/report identifiers retained for tracing. |
| `nclose_occurrence_count` | Sum of the structure's original NClose use counts. |

`SKYPE.STRUCTURE.<number>` identifies a fitted column, with the number equal
to its zero-based index plus one. Virtual structures use
`SKYPE.STRUCTURE.VIRTUAL_INV.<pair_id>`. Rebuilding the candidate matrix can
change fitted-column IDs; these are not globally stable genomic identifiers.

### Structure–NClose usage

`structure_nclose_usage.tsv` contains:

| Columns | Meaning |
| --- | --- |
| `structure_id`, `kind`, `feature_index` | The contributing structure. |
| `nclose_key`, `nclose_id` | The original NClose identity. |
| `occurrence_count` | Number of uses in this structure. |
| `structure_weight`, `structure_weight_N` | Own structure weight before multiplying by the count. |
| `contribution`, `contribution_N` | Count × own weight in raw and normalized units. |

Summing `contribution_N` by `nclose_id` reproduces `nclose_report.tsv` totals.
This table describes original membership; the BND contribution TSV describes
the projected/exported representation.

### Original NClose report

`nclose_report.tsv` contains original Type 1/2/4 entries. Compound MERGE_TYPE4
summary rows belong in the structure report, not this table.

| Columns | Meaning |
| --- | --- |
| `nclose_id` | Original `SKYPE.nclose.<number>` identifier. |
| `start_chr`, `start_pos`, `start_dir`, `end_chr`, `end_pos`, `end_dir` | Catalog coordinates and orientation. BND rows retain anchor-start metadata; join by ID rather than assuming these equal VCF anchor positions. |
| `nclose_cn` | Full nonnegative original-NClose sum, including zero and below-threshold values. |
| `nclose_cn_reason` | `SUPPORTED` for positive support, otherwise `ZERO_SUPPORT`. |
| `preprocessing_reason` | Original exclusion/history status or restored-constituent provenance, independent of the current sum. |
| `kind`, `nclose_key` | Internal representation (`bnd` or `indel`) and original key. |
| `model_cn`, `virtual_cn` | Fitted and virtual parts of `nclose_cn`. |
| `nclose_filter`, `nclose_cluster` and their reason columns | Retained compatibility columns: `NA` / `NOT_RUN_RAW_NNLS` in this native workflow. |

An original path exclusion does not hide support supplied by another structure.
For example, a NClose with `FILTERED_02_NO_ELIGIBLE_PATH` history can have positive
AMP support. Native reports show that sum and retain the history separately;
they do not replace the CN with `-1`.

Existing original NClose IDs are preserved when reconstructing a result's
accounting table. Missing compound constituents receive new IDs, and IDs formerly
used for compound summaries remain reserved as structure provenance. Numbering
can therefore contain gaps. Use `NCLOSE_IDS` and `PARENT_IDS` to trace regenerated
VCF calls; sequential `SKYPE.BND.*` IDs can change when the call set changes.

## BED, Circos and CN lists

The BED columns are `#chrom`, `cordst`, `cordnd`, `type`, `weight (N)`,
`nclose_id`, `weight_scope`, and `source_id`.

- `weight_scope=NCLOSE` means projected NClose support. Original BND views have
  endpoint anchor intervals; PATH_SPLIT views have breakpoint point intervals.
  Ordinary Type 4 rows have reference spans.
- `weight_scope=STRUCTURE` means the structure's own normalized coefficient or
  virtual estimate. AMP, MERGE_TYPE4, virtual-inversion and centromere summaries
  use this scope. `source_id` joins to the structure report.

BED retains both constituent views and compound summaries. `nclose_id` can list
multiple constituent IDs on a summary, and is `.` for a centromere-only feature.

Circos uses the same display event list as BED. NClose links include all
structure contributions after PATH_SPLIT transfers; AMP contributions are not
subtracted from those links. Compound summary links use their own structure
weights. Observed and predicted depth tracks remain those of the fitted model;
the additive virtual output convention does not refit either track.

`cn_data.pkl` stores this tuple:

```text
(inversion_NClose_CN, translocation_NClose_CN, Type4_NClose_CN, new_telomere_CN)
```

The first three lists contain each qualifying original NClose once, using its
full sum and the `> 0.1` gate. They contain neither compound summary entries nor
PATH_SPLIT children. The telomere list retains its existing telomere-selection
rules and uses path coefficients. The lists contain numbers only; use the TSV
reports when IDs and coordinates are needed.

## VCF-input outputs

VCF input writes `SV_benchmark_result.vcf` instead of the native VCF/BED. It
preserves input records and adds `SKYPE_CN` and `SKYPE_STATUS`; side-specific
measurements also receive `SKYPE_CN_DETAIL` and `SKYPE_STATUS_DETAIL`. Skipped
records retain annotations describing their skipped status.

Parser diagnostics include `vcf_mode_summary.json`, `vcf_mode_summary.tsv`,
`vcf_mode_skipped_records.tsv` and `vcf_mode_orientation_mismatches.tsv`.
See ACCtools' [caller compatibility notes](https://github.com/ACCtools/ACCtools-pipeline#compatible-vcf-inputs)
for input-format details.

This mode retains its existing CN aggregation, NClose-report format and
preprocessing-status behavior. The native `STRUCTURE_SUM` reporting changes
and native structure reports do not apply to it.

## Model and intermediate files

| File or directory | Contents |
| --- | --- |
| `01_nclose_data.pkl`, `nclose_nodes.pkl` | Graph handoff and canonical NClose pairs. |
| `skype_options.json` | Saved preprocessing/graph settings. |
| `stage01_nclose_summary.json`, `stage01_nclose_rejections.tsv` | Ordered preprocessing counts and first-rejection reasons. |
| `censat_endpoint_summary.json`, `censat_endpoint_candidates.tsv` | CEN-SAT partition and endpoint-consistency results. |
| `path_data.pkl`, `contig_pat_vec_data.pkl` | Candidate path traversals and component mappings. |
| `conjoined_type4_ins_del.pkl`, `ecdna_circuit_data.pkl` | Compound Type 4 and AMP constituents. |
| `nclose_event_catalog.pkl`, `nclose_path_usage.pkl`, `nclose_filter_status.pkl` | Original event tracking and preprocessing provenance; native output reconstruction leaves these intact. |
| `tot_loc_list.pkl`, `23_input.pkl` | Matrix-column identities and depth-coordinate metadata. |
| `B.npy`, `weight.npy`, `predict_B.npy` | Observed depth, raw fitted coefficients and reconstructed depth. |
| `structure_nclose_model.pkl` | Versioned native structure membership, original NClose metadata, source signature, and stage-31 weights/projections. |
| `raw_translocation_result.pkl`, `raw_translocation_read_counts.tsv` | Paired-junction evidence used to qualify and estimate virtual structures. |
| `24_raw_rescue/` | D queries, read/OLC evidence, candidates, summaries, installation state and pre-rescue snapshots. |

The shared structure model is built at stage 22 without display thresholds.
Stage 31 joins final coefficients, adds qualified virtual structures, records
PATH_SPLIT membership and saves the normalization unit. Its model membership is
checked against the source signature on reload. Regenerating stage 22 rebuilds
the table for the current matrix columns.

## Regenerating existing results

Use a [stage-31 restart](usage.md#restarts) with the original sample, reference,
inputs and result directory. If the structure model is absent, it is reconstructed
from saved path, node, catalog and circuit artifacts. This does not require
refitting the depth model, but does require those original provenance artifacts.

When updating older native outputs, account for these format/meaning changes:

| Item | Current behavior |
| --- | --- |
| VCF `PARENT_IDS` | Structure IDs; original NClose IDs are in `NCLOSE_IDS`. |
| Virtual support | Added to fitted support for the same original NClose; shown separately in `VIRTUAL_WEIGHT`. |
| AMP summary CN | Own structure CN, without the former extra ×2. |
| NClose report | Actual CN and separate history; compound summaries live in the structure report. |
| BED | Adds `weight_scope` and `source_id`; numeric weights retain full output precision. |
| Circos/CN distributions | Full NClose contributions; CN lists exclude compound summaries and split children. |
| Output selection | Contributions are summed before thresholding; low-weight paths can collectively produce visible PATH_SPLIT calls. |

These changes can alter exported weights, adjacency sets and sequential BND IDs
while leaving the fitted coefficients and predicted depth unchanged. They
describe output accounting and representation, not evidence of improved
breakpoint accuracy by themselves.
