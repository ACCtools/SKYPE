# Outputs and weight interpretation

[README](../README.md) · [Usage](usage.md) · [Methods](methods.md)

Unless a section says otherwise, this document describes native assembly
output. [VCF input](#vcf-input-outputs) and full-assembly input retain their
separate reporting rules.

## Choosing a result file

| File | Unit of a row or record | Purpose |
| --- | --- | --- |
| `SV_call_result.vcf` | A BND mate record or ordinary symbolic DEL/DUP | Variant export. A BND adjacency has two reciprocal records. |
| `SV_call_result.bnd_weights.tsv` | One structure/source NClose's contribution to an emitted BND | Reproduce BND weights after exact geometry merging. |
| `SV_call_result.nclose_bnds.tsv` | One source NClose junction | Trace alignment pairs to merged VCF BND IDs, including unexported geometry. |
| `SKYPE_result.bed` | A NClose endpoint/span or a structure summary span | Inspect genomic intervals with an explicit weight scope. |
| `nclose_report.tsv` | One Type 1/2/4 NClose after exact BND identity reuse | Inspect the full NClose total and representative preprocessing history. |
| `nclose_sources.tsv` | One original source node pair/event | Trace source IDs, unitigs and canonical NClose identities. |
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
coefficient. Exact BND geometry reuses an existing NClose ID, preferring the
compressed-graph entry. There is no distance tolerance: a different boundary,
chromosome or retained side remains a different native identity.

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

AMP, MERGE_TYPE4 and VIRTUAL_INV first resolve to their original constituent
NCloses. The native VCF then projects each eligible source NClose into its
adjacent-alignment BNDs. A depth-model bridge between unitigs is not exported as
an additional inferred junction. Ordinary unmerged Type 4 events remain symbolic
DEL/DUP calls.

BED and Circos also retain the compound structure's summary span. That summary
uses the structure's own weight; its constituent NClose views use their summed
support. A summary and its constituent views are different descriptions of
related evidence and must not be summed as independent events.

### VCF-only NClose decomposition

BED, Circos, `nclose_report.tsv` and CN distributions use the original NClose
endpoints and full structure contributions. The native VCF separately expands
each source NClose's inclusive node interval, sorts its retained alignments by
query coordinates and exports adjacent junctions. It includes CEN-SAT endpoint
unitigs, `subtelomere_cut_contig_*`, `telomere_middle_cut_contig_*` and constituent
NCloses of VIRTUAL_INV. Stage-24 read and OLC rescue retain their original pair,
identified by `rescue_method` or `raw_rescue_read_*`/`raw_rescue_olc_*` owners.
Explicit synthetic debug pairs also retain their specified geometry.

Neighbors on the same chromosome and strand with equal signed reference and
query gaps are normal continuations and are omitted. All other retained
adjacencies are projected without a new size or MAPQ cutoff. Source alignments
come from `contig_data`, not depth PAFs containing inferred reference bridges.
Only the selected NClose interval is expanded, not its entire parent unitig.

For a BND, each structure/source contribution is:

```text
contribution_N = structure_weight_N × source NClose occurrences × BND occurrences in that source
```

For example, an `A -> X -> B` source with CN `2` gives `A-X` and `X-B` CN `2`
each. If a separate raw rescue uses `A-B` with CN `1`, only the outer `A-B`
receives that `1`. Repeated exact junctions within a source retain their true
multiplicity. No projected weight changes the original NClose or depth model.

Contributions are summed before the VCF `> 0.1` threshold, across all source
NCloses sharing an exact BND. Two different source NCloses at CN `0.06` can
therefore produce a visible shared primitive junction even when neither NClose
is individually visible in BED.

### Exact native BND identity reuse and VCF geometry merging

When AMP, MERGE_TYPE4, VIRTUAL_INV or a legacy catalog supplies a source pair
with exactly the same unordered endpoint/retained-side pair as an existing BND,
native accounting reuses that NClose. It prefers compressed-graph IDs, then
the first registered matching source. Restoring an exact alias does not allocate
a new ID. Ordinary symbolic Type 4 events retain their existing identities.

Structure weights are unchanged. Counts from different physical occurrences
add to the shared NClose. Source counts preserve the contribution of each
original pair without moving support to other junctions. The NClose report, BED
and CN lists apply their threshold after identity aggregation.
`nclose_sources.tsv` preserves original source metadata.

VCF merges projected geometry only when both endpoint chromosomes, exact
0-based boundary coordinates and retained sides match, allowing endpoint-order
reversal. There is no distance tolerance: a 1 bp difference, different mate or
different retained side keeps a separate BND. Projection happens per original
source, so exact outer aliases with different interiors do not share their
internal junction contributions. `NCLOSE_IDS` contains shared NClose identities,
while `NCLOSE_KEYS` keeps only contributing original source pairs.

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
| `NCLOSE_IDS` | Shared native IDs joining to `nclose_report.tsv`; available on BND and ordinary Type 4 records. |
| `NCLOSE_KEYS` | BND source node-index pairs, encoded as `first:second`. |
| `PARENT_IDS` | Contributing structure IDs joining to `structure_report.tsv`. |
| `PARENT_WEIGHTS` | Each structure's own normalized weight, in `PARENT_IDS` order. |
| `PARENT_MULTIPLICITY` | Occurrence counts for this exact adjacency within each parent, in the same order. |
| `MODEL_WEIGHT` | Sum of fitted-column occurrence contributions. |
| `VIRTUAL_WEIGHT` | Sum of virtual-structure occurrence contributions. |
| `MODEL_FEATURE_COUNT` | Number of distinct positive fitted columns contributing to the BND. |
| `MODEL_OCCURRENCE_COUNT` | Contributing occurrence count over those fitted columns. |
| `SVCLASS` | BND source classes, such as `NCLOSE`, `AMPLICON`, `MERGED_TYPE4`, `VIRTUAL_INV`. |
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
`occurrence_count`, `feature_weight_N`, `contribution_N`, `weight_method`,
`structure_id`, `nclose_id`, `structure_kind`, `source_nclose_key`,
`nclose_occurrence_count`, `bnd_occurrences_per_nclose` and `junction_indices`.
`nclose_id` identifies the original NClose contributing to the emitted BND.
Rows distinguish each contributing structure/source pair. `occurrence_count`
equals `nclose_occurrence_count * bnd_occurrences_per_nclose`; `junction_indices`
lists the matching 1-based adjacency indices in query order, separated by `;`.
Summing `contribution_N` by `bnd_id` reproduces each BND's `WEIGHT`.
Virtual rows have `feature_index=.`.

`SV_call_result.nclose_bnds.tsv` records the topology independently of positive
support. Its columns are `nclose_id`, `source_nclose_key`, `junction_index`,
`node_a`, `node_b`, `mode`, `chrom_a`, `pos_a0`, `side_a`, `chrom_b`, `pos_b0`,
`side_b` and `bnd_id`. Node pairs are in query order; endpoint triples are in
canonical chromosome/coordinate/retained-side order (`L` or `R`). Modes are
`ADJACENT_ALIGNMENT`, `RAW_RESCUE_OUTER` and `SYNTHETIC_OUTER`. The last two
contain a single preserved outer pair. A `.` BND ID means that geometry did not
pass the VCF weight threshold. A source with zero weight may still map to an
emitted geometry supported by another source; only the contribution TSV records
positive support. Normal continuations are absent. BND ID `SKYPE.BND.k` joins to
the two VCF records `SKYPE.BND.k_1` and `SKYPE.BND.k_2`.

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
| `source_nclose_keys`, `source_occurrence_counts` | Semicolon-separated original keys and occurrence counts, in matching order. Their counts sum to `occurrence_count`. |

Summing `contribution_N` by `nclose_id` reproduces `nclose_report.tsv` totals.
This table describes original membership; the BND contribution TSV describes
the exported BNDs after geometry merging and output thresholding.

### Original NClose report

`nclose_report.tsv` contains Type 1/2/4 entries with exact native BND aliases
sharing one row. Compound MERGE_TYPE4 summary rows belong in the structure report.
Catalog coordinates and preprocessing history describe the selected representative.

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

Existing representative NClose IDs are preserved when reconstructing a result's
accounting table. Missing constituents receive new IDs only when their exact BND
geometry is new. Legacy aliases' IDs and IDs formerly used for compound summaries
remain reserved, so numbering can contain gaps. Use `NCLOSE_IDS`, `NCLOSE_KEYS`
and `PARENT_IDS` to trace regenerated VCF calls; sequential `SKYPE.BND.*` IDs can
change when the call set changes.

`nclose_sources.tsv` joins each `source_nclose_key` and legacy `source_nclose_id`
to the shared `nclose_id` and `canonical_nclose_key`, with `is_alias`, supporting
`contig_names` and discovery `source`. A newly restored alias has no separately
allocated source ID (`.`); its node pair identifies it. This is source provenance,
not an additional variant count or a change to preprocessing acceptance rules.

## BED, Circos and CN lists

The BED columns are `#chrom`, `cordst`, `cordnd`, `type`, `weight (N)`,
`nclose_id`, `weight_scope`, and `source_id`.

- `weight_scope=NCLOSE` means the full original NClose support. BND views have
  endpoint anchor intervals. Ordinary Type 4 rows have reference spans.
- `weight_scope=STRUCTURE` means the structure's own normalized coefficient or
  virtual estimate. AMP, MERGE_TYPE4, virtual-inversion and centromere summaries
  use this scope. `source_id` joins to the structure report.

BED retains both constituent views and compound summaries. `nclose_id` can list
multiple constituent IDs on a summary, and is `.` for a centromere-only feature.

Circos uses the same display event list as BED. NClose links include all
structure contributions; AMP contributions are not subtracted from those links.
Compound summary links use their own structure weights. Observed and predicted
depth tracks remain those of the fitted model;
the additive virtual output convention does not refit either track.

`cn_data.pkl` stores this tuple:

```text
(inversion_NClose_CN, translocation_NClose_CN, Type4_NClose_CN, new_telomere_CN)
```

The first three lists contain each qualifying original NClose once, using its
full sum and the `> 0.1` gate. They exclude compound summary entries. The telomere
list retains its existing telomere-selection rules and uses path coefficients.
The lists contain numbers only; use the TSV reports when IDs and coordinates
are needed.

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
| `structure_nclose_model.pkl` | Versioned native structure membership, original NClose metadata, source signature, and stage-31 weights. |
| `raw_translocation_result.pkl`, `raw_translocation_read_counts.tsv` | Paired-junction evidence used to qualify and estimate virtual structures. |
| `24_raw_rescue/` | D queries, read/OLC evidence, candidates, summaries, installation state and pre-rescue snapshots. |

The shared structure model is built at stage 22 without display thresholds.
Stage 31 joins final coefficients, adds qualified virtual structures and saves
the normalization unit. Its model membership is checked against the source
signature on reload. Regenerating stage 22 rebuilds the table for the current
matrix columns.

## Regenerating existing results

Use a [stage-31 restart](usage.md#restarts) with the original sample, reference,
inputs and result directory. If the structure model is absent, it is reconstructed
from saved path, node, catalog and circuit artifacts. This does not require
refitting the depth model, but does require those original provenance artifacts.

Structure-model version 3 removes the former `PATH_SPLIT` representation and its saved
membership fields. Stage 31 rebuilds older cached structure models from the
original intermediate artifacts, retaining the fitted coefficients and predicted
depth. Former split children disappear; their path contributions remain on the
original NClose. The current VCF-only projection is calculated on export and
does not add split children or membership fields to the cached structure model.
It works with existing version-3 models and original node artifacts. Exported
BND counts, coordinates, weights and sequential IDs can change. Existing result
files are updated only when stage 31 is rerun.

The BND contribution TSV no longer has the redundant `source_occurrence` column;
use `nclose_id` instead. The other columns retain their order. Update consumers
that select the removed column or address later columns by position.
`source_occurrence_counts` in `structure_nclose_usage.tsv` remains unchanged.
The VCF projection appends four source/junction columns to `bnd_weights.tsv`
and adds `nclose_bnds.tsv`; use column names when reading these files.

When updating older native outputs, account for these format/meaning changes:

| Item | Current behavior |
| --- | --- |
| VCF `PARENT_IDS` | Structure IDs; original NClose IDs are in `NCLOSE_IDS`. |
| Virtual support | Added to fitted support for the same original NClose; shown separately in `VIRTUAL_WEIGHT`. |
| AMP summary CN | Own structure CN, without the former extra ×2. |
| NClose report | Actual CN and separate history; compound summaries live in the structure report. |
| BED | Adds `weight_scope` and `source_id`; numeric weights retain full output precision. |
| Circos/CN distributions | Full NClose contributions; CN lists exclude compound summaries. |
| Output selection | Contributions are summed before thresholding; low-weight paths can collectively produce a visible original NClose. |
| Exact native BND aliases | Reuse a shared NClose ID/report row; original source mappings are in `nclose_sources.tsv`. Version-1 and version-2 cached structure models are rebuilt as version 3. |
| VCF-only decomposition | Projects each eligible source NClose into adjacent junctions, then merges exact endpoint/side pairs; read and OLC rescue keep their outer pair. |

These changes can alter exported weights, adjacency sets and sequential BND IDs
while leaving the fitted coefficients and predicted depth unchanged. They
describe output accounting and representation, not evidence of improved
breakpoint accuracy by themselves.
