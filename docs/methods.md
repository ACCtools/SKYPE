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

Native BND and symbolic DEL/DUP records use missing VCF `QUAL` (`.`). SKYPE does
not fit a calibrated variant-error probability, so a fixed Phred score would
overstate confidence. `WEIGHT` describes conditional model dosage; it is not a
probability that the junction or its chromosome assignment is correct. Mapping
quality, local molecule evidence, full-chain linkage and normal-sample exposure
are separate forms of evidence. An unmeasured quantity must not be read as zero.

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

Reference connections between unitigs follow their actual traversal strands,
including when their alignments overlap or one contains the other. A unitig's
own observed alignment chain may still change strand at an NClose. Telomere
entries have matching reverse exits, and both terminal anchors contribute to
the path depth.

Compound Type 2 candidates are tested in both unitig orders. Their reported
DEL/DUP span and breakend directions use the walked outer junction boundaries,
including reverse-complement traversals. The label is separate from the depth
baseline sign: overlapping duplication anchors subtract their reference union
once to retain the extra copy in the overlap. The catalog's `depth_base_sign`
records this arithmetic for stage 22; original alignment strands remain in PAF.
ecDNA length counts the walked unitig pieces and reference bridges, with overlap
counted once. Reference distances jumped across an NClose are not circle length.

Runs created before these path-geometry corrections need to restart at stage 01
to rebuild candidates, paths, depth columns and their report metadata together.

## Pipeline overview

ACCtools prepares references, assemblies, alignments and read depth, then calls
[pipeline.py](../pipeline.py). The native orchestrator owns stage commands,
completion checks, option routing and partial restarts.

| Stage | Main responsibility |
| --- | --- |
| 00 | Depth normalization when a reference-depth correction is supplied, per-sample GC correction and the B-allele-fraction depth per copy (1N). |
| 01 and preprocessing helpers | Unitig/telomere processing, NClose classification and filtering, raw-read validation and compound candidate discovery. |
| 10 | Breakend graph construction and candidate chromosome-path search. |
| 11 | Reference-ratio outlier and Type 4 feature preparation. |
| 21 | Candidate path/component alignments and their depth vectors. |
| 22 | The observed-depth system, column metadata and structure–NClose membership. |
| 23 | One raw non-negative least-squares fit over the candidate matrix columns. |
| 24 | Optional raw-read rescue, followed by one refit if new candidates are retained. |
| 25 | Karyotype mode only: greedy NClose-filter and karyotype-cluster refits. |
| 30 | Virtual SKY diagram and ISCN-style karyotype for each fit. |
| 31 | Variant output, structure/NClose reports, CN lists and the depth plot for each fit. |

Native fitting has no normal-chromosome prior. Candidate preprocessing and graph
filters apply before fitting, and the optional rescue refit is a second model
fit after candidate augmentation. The unsuffixed outputs always come from the
raw fit; the stage-25 refits only add the `_filter` and `_cluster` result sets.

## Assembly preprocessing

[01_Preprocess_NClose.py](../01_Preprocess_NClose.py) owns the preprocessing
sequence and ecDNA discovery. In assembly mode it reads the `--alt` unitig PAF
and the last `--original-paf-loc` value. The positional primary PAF is retained
for command/file-name compatibility; the resulting primary-named `*.ppc.paf`
contains unitig-only rows.

Despite its filename, `*.ppc.paf` is an internal node table, not standard PAF.
Its first nine fields retain the query/reference coordinate layout, but the
following fields hold SKYPE metadata: field 10 (1-based) is the original MAPQ,
field 11 is the node type, and fields 12–13 are the source's start/end node
indices. In particular, field 12 is **not** PAF MAPQ. Use raw or selected
alignment PAFs for standard PAF tooling and the saved `contig_data` schema for
node-level analysis.

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

1. Same-chromosome endpoint pairs are excluded, including opposite strands.
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

`realigned_consistent` means that the reported qualifying alignments agree on
chromosome and strand under these alignment options. It does not establish
unique sequence origin, exact breakpoint accuracy, or independent biological
validation. Alternative mappings may depend on the alignment options and query
context. Reusing the same assembly sequence in the two checks is not an
independent molecular assay.

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

Touching reference anchors are interpreted in the walked direction. An
inward-facing contact whose entry and exit coincide contributes zero depth;
the two complete anchor alignments must not be added as reference sequence.
Ordinary consecutive contacts retain both alignments. This rule applies to
reference bridges in chromosome paths and ecDNA circuits, while a unitig's
observed NClose connection remains a separate kind of adjacency. A corrected
zero-depth feature supplies no dosage evidence and does not establish absence
of its constituent junctions. Rebuild from stage 21 to apply this depth-rendering
correction to an otherwise unchanged saved graph.

### Graph traversal

Stage 01's ecDNA graph and stage 10's NClose graph distinguish positive interval
overlap from boundary contact when connecting different NClose sources. Only
positive overlap can bypass the coordinate ordering check; ordinary consecutive
touching anchors satisfy the ordering check and retain both anchors. Inward
contacts with zero retained span are excluded during graph construction. The
ecDNA length gate also rejects any circuit containing a nonpositive reference
bridge, even if its other components give it a positive total length.

Both graph builders independently assert that each generated reference edge
stays on the same chromosome and walked strand and retains positive reference
span. A violation stops construction before enumeration and reports node indices,
source names, coordinates, PAF strands, walked strands, and span. Same-source
NClose jumps are handled separately and are exempt from this check. Restart from
stage 01 to regenerate saved ecDNA candidates and downstream results.

Stage 10 permits two visits to a reference chromosome-end telomere anchor,
including a synthetic `virtual_telo_*` anchor for a missing reference end.
Qualification uses the reference telomere annotation and matching terminal label;
an internal new telomere does not qualify from a `TELCON` label alone. Ordinary
NClose nodes cannot be revisited, and a third reference-end-anchor visit is rejected.

Same-terminal paths such as `chr1f → chr1f` are searched when qualifying anchors
exist and require reference-end anchors at both ends. Graph filters and the
configured search limits apply to the candidate paths.

Path deduplication uses the exact ordered source-node walk, including each
node's traversal direction and the terminal labels. Expanded graph counters
are search state and are not part of this identity. Equal chromosome-total
lengths or equal depth predictions do not justify merging different junctions,
source chains or their ordering. This replaces the former 50 kb tolerance on
chromosome-total lengths, which could discard different deletion loci or
alternative junction combinations depending on enumeration order.

`path_search_diagnostics.json` records the retained graph's search limits and
terminal pairs that reached the per-pair path limit. Reaching that limit means
exhaustive enumeration is unproven; it does not establish that additional
paths exist or that omitted paths have no biological support. Runs made with
the former path compression need to restart at stage 10 or earlier before
comparing downstream fits or structural outputs.
The existing default budgets can be changed explicitly with
`--per-pair-path-limit` and `--total-path-limit` through `--option_skype`.
For controlled comparisons, fix `--limit-combinations` as well so a change in
candidate count cannot silently change the graph's permitted rearrangement
complexity through the automatic budget fallback.

## Depth fitting and structure weights

Stage 00 smooths the log CHM13 correction factors with a Gaussian sigma of one
100 kb bin, reflect edges and truncation at four sigma. It preserves the chrY
and large-CEN-SAT exemptions and does not smooth across coordinate gaps or
chromosome boundaries. Sample depth itself is not convolved.

Stage 00 then removes the sample's own GC response, the part the CHM13 reference
library does not share ([depth_gc.py](../depth_gc.py)). Log2 depth is modeled as a
segment level plus a curve of the 100 kb window GC
(`public_data/chm13v2.0_gc.win.tsv.gz`): segments come from an exact optimal
partition of the corrected depth, refined over rounds, and the curve is the
median per 0.01 GC bin with at least 30 windows from 10 segments, flat beyond the
supported GC range. The curve is shrunk by its agreement between the two halves of
the chromosomes and applied to every window except chrY and large CEN-SAT
exemptions; the median depth of all windows is kept. Too few training windows, a
failed or unreproducible fit, or a GC table of another reference leave the depth
uncorrected with the reason in `depth_scale.json`; `--no-gc-correction` skips it.

The depth per copy (1N) comes from depth and B-allele fractions
([baf_ploidy.py](../baf_ploidy.py)), not from the assumption that the median window
holds two copies. Heterozygous sites are counted from the sample BAM in the middle
20 kb of each fitted window, outside segmental duplications and CEN-SAT. Windows
are segmented by depth, het-site rate and allele fraction; segments that are
balanced (1:1), one third, one quarter or credibly LOH anchor the scale, and 1N is
the largest depth scale at which they sit on whole copy numbers (the fewest copies
consistent with the anchors). A whole-genome doubling therefore stays invisible.
The 1N, its status (`supported`, `ambiguous_scale`, `scale_sensitive`, `poor_fit`,
`insufficient_anchors`, `insufficient_data`) and diagnostics are written to
`depth_scale.json`. Stages 30 and 31 use this 1N for every copy number they show or
report whenever it exists, whatever its status, and median depth / 2 otherwise; all
cuts on small values (0.1N, 0.2N, 5% of the median) stay in median depth / 2. A
`depth_scale.json` whose row count, median or sum of depths differs from the depth
table a stage reads is ignored, and stage 00 removes the old file before it runs. The
depth viewer (`skype_view.py`) uses the same 1N for its 1N and 4N rings and CN values
when the depths it shows (`B.npy`) are those of the recorded table.
`--no-baf-n` skips the estimate.

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

## Karyotype-mode refits

[25_cluster_weight.py](../25_cluster_weight.py) runs in karyotype mode, the
default unless variant mode is requested, the input is a VCF or the graph has
more than 1000 NClose nodes. Both refits use the Julia NNLS solver on the
stage-22 matrix; the filter refit starts from the raw-fit coefficients.

The filter refit considers cross-chromosome NCloses that carry at least 1% of
the median depth and have at least one non-CEN-SAT side without a depth step
between its 1 Mb flanks (a step needs a ratio above 1.2 and a difference of at
least 5). For each candidate the
matrix is refitted without the paths that use it. Let the error envelope of a
chromosome be the largest absolute cumulative difference between predicted and
observed depth along it. Candidates whose single removal keeps every envelope
within 1.5 times the raw fit are tried in order of that ratio. A removal is kept
when, together with the removals already kept, no envelope exceeds twice the
raw fit. chrY is ignored when it lacks coverage.

The cluster refit starts from the filter fit. It summarizes paths above 0.1N
and the pure telomere-to-telomere path of each chromosome as chromosome and
strand segments. A path joins the first cluster whose representative shares at
least 80% of the longer length after the best shift (90% when the terminal
segments differ). Each cluster above 0.2N in total contributes its path with
the fewest NCloses. Type 4 events longer than 1 Mb and ecDNA above 0.2N and all
centromere fragments are kept as well, restricted to columns the filter kept. After the refit, coefficients at or below
0.1N are set to zero.

Stage 25 records a reason in `nclose_filter_status.pkl` for every source NClose
a fit removes (`FILTERED_24_GREEDY_DEPTH`, `FILTERED_24_COFILTERED_PATH`,
`FILTERED_24_CLUSTER_SELECTION`, `FILTERED_24_CLUSTER_MIN_WEIGHT`), including
merged Type 4/ecDNA constituents outside the stage-23 catalog. A source passes
while any column carrying it stays in the fit, even if preprocessing removed its
own paths. When the 0.1N floor removes the only positive support, the source
keeps the minimum-weight reason even if an exact-zero carrier stays selected.
Stage 31 combines the sources sharing a canonical NClose: positive support,
including qualified virtual support, is `PASS`; otherwise a minimum-weight cut
wins, then a source that kept a carrier (`0 / PASS`), then the source removed at
the latest stage. Preprocessing exclusions remain separate history.

[30_virtual_sky.py](../30_virtual_sky.py) draws each fit's paths above 0.2N,
Type 4 events over 10 Mb, virtual inversions and centromere fragments, and writes
the matching `karyotype.txt` labels.

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
At least one primitive breakend coordinate must lie in the query region; read overlap
alone is insufficient. Short internal repeat alignments cannot become selected
terminal anchors. Large CIGAR indels split alignments into pieces while keeping
their original query coordinates and exact per-piece cs strings; indel lengths
are not treated as base-alignment errors.

Default gates include MAPQ `20`, anchor length `500` bp, alignment identity
`0.9`, and `3` distinct supporting molecules. There is no additional
outward-continuity-length gate or rescue VAF cutoff. An insertion at one
unchanged reference coordinate does not add a reference-depth edge.
This exclusion applies to an unmapped insertion. A mapped internal template
adds explicit reference adjacencies and is retained even when its outer anchors
are adjacent or collinear. Reference-continuation fragments outside the first
and last junction are trimmed before selecting the source's outer coordinates.
Candidate admission keeps the 1 kb query/reference discrepancy requirement for
same-strand adjacencies; a collection of smaller discrepancies alone does not
become a compound rescue. Chromosome or orientation changes remain eligible.

- `read` uses complete raw molecules, primary/SA discovery and whole-reference
  realignment, then groups molecules by their complete ordered primitive
  junction chain. Equal outer anchors alone cannot pool different interiors.
  Reverse-complement chains are equivalent; exact reference-continuation
  fragments do not add junctions.
- `olc` uses [local_olc.py](../local_olc.py): canonical minimizer seeds, banded
  overlap alignment, an oriented graph, transitive reduction, paths ending at
  branches and alignment-based consensus. Fully covered reads remain evidence;
  tips and bubbles are not pruned. Unitigs are realigned to the whole reference.
  Every primitive junction, including internal junctions, requires distinct
  raw-molecule support. Different molecules can support different junctions of
  a longer path; OLC layout supplies its inferred linkage. This does not by
  itself establish unique chromosome placement or a uniquely resolved path.
  Per-junction read names/counts and complete-chain molecule support are
  reported separately from outer-pair and assembly input-read counts.
  Passing the per-junction gate establishes local junction evidence; it does
  not verify clone-specific linkage or repeated-junction multiplicity.

Standalone single-junction, same-chromosome, same-strand DEL/DUP candidates must span at least
100 kb, inclusive. Smaller candidates are reported as `below_min_indel_span`.
The separate alignment-splitting threshold is 1 kb, so short internal indels can
still support a compound NClose. Accepted standalone DEL/DUP chains enter as
Type 4. Compound sources retain their primitive junctions through the BND
route, including same-chromosome, same-strand templates. DEL and DUP are compared separately
with existing Type 4/catalog events using the stage-11 10 kb endpoint tolerance.
Other chains use the BND discovery/compression route.
For compound sources, compression also requires equivalent ordered primitive
chains within the existing compression tolerance, so equal outer anchors
cannot merge different internal templates.

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
