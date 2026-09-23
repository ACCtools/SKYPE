# Methods

[README](../README.md) · [Usage](usage.md) · [Outputs](outputs.md)

This document describes the native assembly workflow. VCF input supplies its
candidate junctions from a VCF, and complete-assembly input uses a separate
workflow described in [usage](usage.md#input-workflows).

## NCloses and structures

An NClose is the original junction unit tracked through SKYPE's assembly
processing and depth model. Its alignment chain can contain internal pieces;
the depth model, BED, plots and NClose reports retain its original endpoints.
Only the native VCF projects eligible source NCloses into primitive BNDs between
adjacent retained alignments in query order. It includes CEN-SAT endpoint-route,
telomere-derived and virtual-inversion NCloses, as well as read and OLC rescue.
Same-chromosome/same-strand neighbors with equal reference and query gaps are
normal continuations and are omitted; no new size or MAPQ cutoff is
applied by this export. Ordinary Type 4 DEL/DUP and VCF-input annotation retain
their existing representations.
Distinct source pairs with exactly equal BND endpoints and retained sides reuse
one native accounting identity, preferring an existing compressed-graph NClose.
Source-specific alignments and occurrence counts remain separate.
VCF projection uses those source counts before merging exactly equal BND
endpoint/retained-side pairs. A structure's contribution to one BND is its weight
times its source-NClose count times that BND's count within the source chain.
Projection never creates a junction between different constituent unitigs or
feeds split events back into the graph, matrix or fitted depth.
This does not change the graph's broader spatial compression or AMP eligibility.

| Term | Role |
| --- | --- |
| Type 1 NClose | Junction connecting different reference chromosomes. |
| Type 2 NClose | Intrachromosomal inversion-like junction. |
| Type 4 NClose | Indel-like reference-span change, represented by a deletion or duplication depth feature. |
| Path | A candidate chromosome traversal through the breakend graph; it can use several NCloses. |
| Independent Type 4 | A separate fitted feature for a Type 4 event, in addition to any paths that use that event. |
| MERGE_TYPE4 | A compound depth feature formed from Type 2 NCloses; its constituent NCloses remain the original junctions. |
| AMP | An ecDNA/amplicon circuit feature using its two constituent NClose pairs. |
| VIRTUAL_INV | A qualified paired-junction structure with a raw-read-derived depth estimate. Its two original NCloses are used in native output. |

Type 1, Type 2 and Type 4 are all NCloses. A compound structure is a model of
how NCloses are used, not an additional original NClose. A reference bridge
used to explain depth is not, by itself, evidence for an extra junction.

## Pipeline overview

ACCtools prepares references, assemblies, alignments and read depth, then calls
[pipeline.py](../pipeline.py). The native orchestrator owns stage commands,
completion checks, option routing and partial restarts.

| Stage | Main responsibility |
| --- | --- |
| 00 | Depth normalization when a reference-depth correction is supplied. |
| 01 and preprocessing helpers | Unitig/telomere processing, NClose classification and filtering, raw-read validation and compound candidate discovery. |
| 10 | Breakend graph construction and candidate chromosome-path search. |
| 11 | Reference-ratio outlier and Type 4 feature preparation. |
| 21 | Candidate path/component alignments and their depth vectors. |
| 22 | The observed-depth system, column metadata and structure–NClose membership. |
| 23 | One raw non-negative least-squares fit over the candidate matrix columns. |
| 24 | Optional raw-read rescue, followed by one refit if new candidates are retained. |
| 31 | Variant output, structure/NClose reports, CN lists and the depth plot. |

Native fitting has no normal-chromosome prior, post-NNLS NClose-filtering pass,
or post-NNLS clustering pass. Candidate preprocessing and graph filters still
apply before fitting. The optional rescue refit is a second model fit after
candidate augmentation, not a post-fit weight-filtering stage.

## Assembly preprocessing

[01_Preprocess_NClose.py](../01_Preprocess_NClose.py) owns the preprocessing
sequence and ecDNA discovery. In assembly mode it reads the `--alt` unitig PAF
and the last `--original-paf-loc` value. The positional primary PAF is retained
for command/file-name compatibility; the resulting primary-named `*.ppc.paf`
contains unitig-only rows.

The graph handoff, `01_nclose_data.pkl`, contains exactly `contig_data`,
`nclose_nodes` and `telo_contig`. Downstream stages read canonical NClose pairs
from `nclose_nodes.pkl`. Ordered stage counts and first-rejection reasons are
saved in `stage01_nclose_summary.json` and `stage01_nclose_rejections.tsv`.

### CEN-SAT endpoint processing

Before trimming, raw unitigs are partitioned using their two outer aligned
terminal bases against all CEN-SAT BED intervals, interpreted as 0-based,
half-open. This route assignment is fixed even if later trimming changes an
endpoint label.

Non-CEN-SAT-pair unitigs follow the standard preprocessing, terminal-pair
selection, coordinate/direction clustering and filters. One terminal pair is
retained per accepted Type 1/2 unitig.

CEN-SAT-pair unitigs are processed independently by
[censat_endpoints.py](../censat_endpoints.py):

1. Same-chromosome, same-strand endpoint pairs are excluded.
2. Raw primary and secondary alignments covering at least 50% of each original
   end chunk must agree on that end's chromosome and strand.
3. Passing ends are traced through `xi:Z:P_<index>` to raw PAF query intervals,
   extracted without reverse complementation, and realigned to the whole reference.
4. Query offsets are restored and the same original-chunk 50% consistency test
   is applied again. `A_` traces, missing traces and inconsistent traces are rejected;
   missing qualifying alignments do not count as agreement.

The endpoint realignment uses:

```text
minimap2 --cs -x asm20 --no-long-join -r2k -K10G -N 5000 -p 0.5
```

Accepted unitigs keep their original endpoint chunks and internal alignment
rows. They bypass the standard route's trimming, clustering and filters, and
are combined with that route's results before downstream graph processing.
Their graph breakpoints use junction-facing coordinates. CEN-SAT-locus pair
deduplication and the former BOTH-CEN-SAT terminal filter are not applied to
this independent route.

The orchestrator prepares `<sample>.utg.censat_endpoints/` beside the source
PAFs and passes `--censat-endpoints-dir` to stage 01. This directory is required
for direct assembly-stage invocation. See [cache behavior](usage.md#caches)
for reuse and invalidation rules.

### Anchor and telomere selection

After merging the two routes, the default `--nclose-min-ref-span 1000` removes
an assembly pair if either endpoint has `ref_end - ref_start <= 1000` bp.
Setting the option to `0` disables that cutoff. Candidates validated by the
independent CEN-SAT route (`origin == "censat_endpoint_consistency"`) are exempt;
a CEN-SAT label alone does not grant an exemption. Stage-24 rescue candidates
are unaffected by this stage-01 option.

Native and VCF-input preprocessing remove new telomere connections whose
telomere-facing aligned base lies inside chromosome-end CEN-SAT. This includes
fallback connections. Reference telomere anchors are retained. The separate
full-assembly workflow does not enter this filtering logic.

## Compound features and graph paths

### MERGE_TYPE4 construction

Stage 01 builds `conjoined_type4_ins_del.pkl` from same-chromosome Type 2
NCloses that survive standard NClose filtering, before the independent CEN-SAT
endpoint merge. For oriented NCloses `A → B` and `C → D`, the inner directions
must agree. Candidates then require either:

- `distance(B, C) < 3 Mb`; or
- two inversion junctions with both `distance(A, C) < 3 Mb` and
  `distance(B, D) < 3 Mb`.

Distances are gaps between reference alignment intervals; overlapping intervals
have distance zero. Reverse-complement combinations are also checked when the
NClose's reference-start span is at least 100 kb. In either layout the `B → C`
bridge must walk forward: along the strand B is walked on, C's breakend must lie
at or ahead of B's. Mirrored reverse-complement combinations of the same two
junctions pass the unsigned distance and direction checks, and stage 21 would
fill their bridge backwards. The outer `abs(reference_ratio - 1) > 0.1`
condition applies, with no minimum estimated outer span.

Stages 11 and 21 preserve internal alignment pieces, so a long inversion can
retain its middle depth while the signed correction feature represents flanking
deletions. The saved insertion/deletion circuits are indexed with insertion
circuits first, then deletion circuits.

Every row in a stage-11 or stage-21 depth PAF is an alignment piece that a path
selected, so it is written as `tp:A:P`. alignasm keeps minimap2's `tp:A:S` tag
on query pieces that no primary row covers, and PanDepth skips such rows by
default.

### Graph traversal

Stage 10 permits two visits to a reference chromosome-end telomere anchor,
including a synthetic `virtual_telo_*` anchor for a missing reference end.
Qualification uses the reference telomere annotation and matching terminal label;
an internal new telomere does not qualify from a `TELCON` label alone. Ordinary
NClose nodes cannot be revisited, and a third reference-end-anchor visit is rejected.

Same-terminal paths such as `chr1f → chr1f` are searched when qualifying anchors
exist and require reference-end anchors at both ends. Graph filters and the
configured search limits apply to the candidate paths.

## Depth fitting and structure weights

Stage 00 smooths the log CHM13 correction factors with a Gaussian sigma of one
100 kb bin, reflect edges and truncation at four sigma. It preserves the chrY
and large-CEN-SAT exemptions and does not smooth across coordinate gaps or
chromosome boundaries. Sample depth itself is not convolved.

Stages 21–22 construct a reference-depth vector for each candidate feature.
Paths, independent Type 4, MERGE_TYPE4, AMP and centromere fragments have
separate fitted coefficients. Stage 23 solves the raw NNLS system; predicted
depth is the sum of feature-depth vectors multiplied by their coefficients.

The structure–NClose table records how many times each feature uses each
original NClose. Stage 31 combines that table with final coefficients and
qualified VIRTUAL_INV estimates. It does not feed aggregated NClose totals
back into the NNLS coefficients. Centromere-only features have empty NClose
counts, while telomere support uses path-end counts and those same coefficients.
See [weight accounting](outputs.md#weight-accounting) for the reporting formula.

## One raw-read rescue pass

[24_raw_nclose_rescue.py](../24_raw_nclose_rescue.py) examines differences
between valid observed and predicted depth bins. D detection uses exact dynamic
programming with penalty `96` and a minimum of `8` bins per segment. Local
evidence uses up to `20` bins per side; full adjacent-segment evidence is used
when local evidence is insufficient.

### Query regions

The default `adaptive` query mode fits step and flat-ramp-flat location models
to observed depth in contexts of 1.5, 2 and 3 Mb per side, clipped to adjacent
DP segments. Each plateau needs at least five bins. Models within `6` of the
best standardized SSE + parameter-count × log(bin-count) define a location
envelope. Its union across contexts and the original D point is padded by
200 kb and clipped to the valid run; overlapping queries are merged.

This envelope describes location uncertainty, not a calibrated confidence
interval. `fixed` mode uses the original point ±200 kb queries.
`query_models.json` records the models and query coordinates. Extraction caches
are reused only when query coordinates match.

### Read and OLC evidence

Both methods use reliable terminal alignments and native NClose compression.
At least one breakend coordinate must lie in the query region; read overlap
alone is insufficient. Short internal repeat alignments cannot become selected
terminal anchors. Large CIGAR indels split alignments into pieces while keeping
their original query coordinates and exact per-piece cs strings; indel lengths
are not treated as base-alignment errors.

Default gates include MAPQ `20`, anchor length `500` bp, alignment identity
`0.9`, and `3` distinct supporting molecules. There is no additional
outward-continuity-length gate or rescue VAF cutoff. An insertion at one
unchanged reference coordinate does not add a reference-depth edge.

- `read` uses complete raw molecules, primary/SA discovery and whole-reference
  realignment, then selects each molecule's two reliable outer alignments.
- `olc` uses [local_olc.py](../local_olc.py): canonical minimizer seeds, banded
  overlap alignment, an oriented graph, transitive reduction, paths ending at
  branches and alignment-based consensus. Fully covered reads remain evidence;
  tips and bubbles are not pruned. Unitigs are realigned to the whole reference.
  Both terminal junctions require raw-molecule support, which can come from
  different molecules for a longer path. Paired outer-read support is reported
  separately; assembly input-read count is not junction support.

Standalone same-chromosome, same-strand DEL/DUP candidates must span at least
100 kb, inclusive. Smaller candidates are reported as `below_min_indel_span`.
The separate alignment-splitting threshold is 1 kb, so short internal indels can
still support a compound NClose. Accepted standalone DEL/DUP chains enter as
Type 4 with their internal pieces retained. DEL and DUP are compared separately
with existing Type 4/catalog events using the stage-11 10 kb endpoint tolerance.
Other chains use the BND discovery/compression route.

After compression retains a new NClose, the orchestrator installs the augmented
PAF/catalog and reruns stages 10–23 once. `24_raw_rescue/round.json` prevents a
completed pass from adding another round on restart. Stage 31 reports the final
fit. The rescue candidate table records `handoff_type` and `indel_event_type`;
its summary separates `added_type4_count` and `added_bnd_count`.

## Extending preprocessing

[nclose_preprocess.py](../nclose_preprocess.py) exposes ordered
`ContigPipelineStage` and `NClosePipelineStage` extension points.
`default_stage01_pipeline()` returns the active sequence.
`with_contig_stage(...)`, `with_nclose_stage(...)` and
`make_nclose_filter_stage(...)` allow placement relative to named boundaries
without changing the stage-01 driver. Options and restart behavior are described
in [usage](usage.md#native-options).
