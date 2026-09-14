# Graph-based cancer structural-variant caller

[![Tests](https://github.com/ACCtools/SKYPE/actions/workflows/tests.yml/badge.svg)](https://github.com/ACCtools/SKYPE/actions/workflows/tests.yml)

**SKYPE** discovers structural variants from assembly alignments and assigns
copy-number support from observed read depth. It can also import an existing
VCF and annotate its records with SKYPE depth support.

The native pipeline:

- classifies and compresses assembly segments around NClose junctions;
- builds a breakend graph and enumerates candidate chromosome paths;
- constructs one observed-depth matrix for all candidate paths;
- fits every matrix column with one raw non-negative least-squares solve per pass;
- optionally discovers missing NCloses from depth-triggered raw-read regions and refits once; and
- emits copy-number-weighted variants, an NClose report, and a coverage plot.

There is no karyotype/variant mode switch, normal-chromosome prior, post-NNLS
NClose filtering, or clustering stage.

### Stage 01/10 pipeline

ACCtools `SKYPE.py` prepares references and assembly alignments, then calls
`pipeline.py` once for native assembly and VCF-input runs. `pipeline.py` owns
normalization, CEN-SAT endpoint preparation, stages 01/10/11/21/22/23/24/31,
completion checks, and partial restarts. Existing ACCtools options and output
locations are preserved; `--print_args` still prints the individual stage
commands. Run `python deps/SKYPE/pipeline.py --help` from the workspace root
for direct invocation with prepared input paths. Full-assembly runs continue
to use `full_assembly_pipeline.py`.

`01_Preprocess_NClose.py` now owns contig/telomere preprocessing, NClose
post-processing, ecDNA discovery, and the exact three-field
`01_nclose_data.pkl` handoff consumed by `10_Graph_Find_Paths.py`. VCF input is
usable through the split stages. In assembly mode the positional primary PAF
is retained only for command/file-name compatibility: stage 01 reads the
`--alt` unitig PAF and the last `--original-paf-loc` value, and writes
unitig-only rows under the primary-named `*.ppc.paf`.

Before any trimming, raw unitigs are partitioned by their two outer aligned
terminal bases against **all** CEN-SAT BED intervals (0-based, half-open).
Non-CEN-SAT-pair unitigs follow the existing preprocessing, one terminal pair
per retained type-1/2 unitig, and coordinate/direction clustering and filters.
Route membership is fixed even if later trimming changes endpoint labels.

CEN-SAT-pair unitigs follow `censat_endpoints.py` independently. Same-chromosome,
same-strand endpoints are excluded. Both remaining ends must have only matching
chromosome/strand assignments among raw P+S alignments covering at least 50%
of each original end chunk. Passing ends are traced through `xi:Z:P_<index>`
to raw PAF query intervals, extracted
without reverse complementation, and realigned to the whole reference with
`minimap2 --cs -x asm20 --no-long-join -r2k -K10G -N 5000 -p 0.5`.
The same original-chunk 50% test is repeated after restoring query offsets.
`A_`, missing, or inconsistent source traces are excluded and reported.
Missing qualifying alignments never count as consistent evidence.

Accepted unitigs retain their original endpoint chunks and internal rows;
graph breakpoints use the existing junction-facing coordinate convention.
They bypass legacy stage-01 trimming, clustering, and filters, and are merged
only after the other route finishes. The old BOTH-CEN-SAT terminal/MAPQ/strand
filter and CEN-SAT-locus pair deduplication have been removed. Stage 10 and
later graph/path processing remain in effect for the merged candidates.

SKYPE `pipeline.py` prepares `<sample>.utg.censat_endpoints/` beside the source
PAFs and passes `--censat-endpoints-dir` to stage 01 (required for direct
assembly-stage invocation). This cache contains the partition, candidate
metadata, source-coordinate manifest, extracted FASTA, and realignment PAF.
It is reused across SKYPE result directories when assembly/PAF/reference/BED
signatures, workflow implementation, minimap2 version, and options match.
`--skype_force` alone reuses the alignment cache; alignment `--force` rebuilds
it. VCF and full-assembly modes keep their separate input paths.

Default SKYPE and VCF-input runs exclude new telomere connections whose
telomere-facing aligned base lies inside chromosome-end CEN-SAT. Stage 01
removes these candidates before graph handoff, including fallback telomere
edges; reference telomere anchors are retained. Full-assembly input uses its
separate workflow and does not enter this filtering logic.

Downstream stages read canonical NCloses from `nclose_nodes.pkl`; the former
endpoint-compression tuple is no longer an external artifact.

Stage 01 exposes two ordered extension points in `nclose_preprocess.py`:
`ContigPipelineStage` for adding another unitig split/augmentation policy and
`NClosePipelineStage` (or `make_nclose_filter_stage`) for candidate filtering
or augmentation. `default_stage01_pipeline()` returns the active sequence, and
new stages can be placed relative to a named boundary with
`with_contig_stage(...)` or `with_nclose_stage(...)` without changing the
stage-01 driver.

Every run also writes `stage01_nclose_summary.json` with ordered stage counts
and `stage01_nclose_rejections.tsv` with the first filter/reason that removed
each candidate. SKYPE `pipeline.py` runs stages 01 and 10 for native assembly/VCF input.

Pass native pipeline options through `--option_skype` (or `--option-skype`):

```bash
bash run.sh --option_skype="--skip_bam_analysis --add_indel_graph" HCC1937
```

Stage 01 parses the string using `skype_options.py`, applies preprocessing
options, and saves the stage settings in `skype_options.json`. Stage 10 loads
its graph options from that file. A new stage-01 run replaces saved settings;
a stage-10 restart reuses them. Explicit stage-10 CLI options override saved
values. A direct stage-10 `--option_skype` replaces the saved graph options for
that invocation; preprocessing changes require rerunning stage 01.
The graph handoff remains the same three-field `01_nclose_data.pkl`.

Stage 10 permits two visits to a reference chromosome-end telomere anchor,
including a synthetic `virtual_telo_*` anchor for a missing reference end.
It uses the reference telomere annotation plus the matching terminal label;
an internal newly detected telomere does not qualify merely by having a
`TELCON` label. Ordinary NClose nodes keep the no-revisit rule, and a third
visit to a chromosome-end anchor is rejected. Same-terminal paths (for example
`chr1f -> chr1f`) are searched when such anchors exist and require reference-end
anchors at both ends. Other path filters and search limits are unchanged.

### Stage 24: one raw-read rescue pass

Stage 00 smooths the **log CHM13 correction factors**, using a Gaussian sigma
of one 100 kb bin, reflect edges and truncation at four sigma. It preserves the
existing chrY/large-CEN-SAT exemptions and never smooths across coordinate gaps
or chromosome boundaries. Sample depth itself is not convolved.

`24_raw_nclose_rescue.py` applies D detection to the valid stage-23 observed and
predicted depth bins: exact DP (penalty 96, minimum 8 bins), local evidence over
up to 20 bins per side, and full adjacent-segment evidence only when local
evidence is insufficient. D candidate boundaries and evidence gates are retained.
The default `--query-mode adaptive` fits step and flat-ramp-flat location models
to the observed depth in 1.5, 2, and 3 Mb contexts per side, clipped to the two
adjacent DP segments. Each model keeps at least five bins per plateau. Models
within 6 of the best standardized SSE + parameter-count * log(bin-count) define
a transition/location envelope; its union across contexts and the original D
point is padded by 200 kb and clipped to the valid run. Overlapping queries are
merged. This envelope describes location uncertainty, not a calibrated confidence
interval. `--query-mode fixed` reproduces the original point +/-200 kb queries.
`query_models.json` records the models, envelope and old/new query coordinates.
Read/OLC extraction caches are reused only when the query coordinates match.

The two methods use the same query regions, reliable terminal alignment
selection, and native NClose compression. At least one **breakend coordinate**
must lie inside the query; overlap of a read with the query alone is insufficient.
Short internal repeat alignments cannot become the selected terminal anchors.
Large CIGAR insertions/deletions split an alignment into separate pieces before
terminal selection. Their lengths are not counted as base-alignment errors;
the original query coordinates and exact per-piece cs strings are retained.
An insertion at one unchanged reference coordinate does not add a reference-depth edge.
Default alignment gates are MAPQ 20, anchor 500 bp, and three distinct supporting
molecules. There is no extra outward-continuity-length gate or VAF cutoff.
Standalone same-chromosome, same-strand DUP/DEL candidates must span
at least 100 kb, matching `VCF_TYPE4_MIN_SPAN`. Smaller candidates are recorded
as `below_min_indel_span` before compression and are not installed. The 1 kb
alignment-splitting threshold remains separate, so short internal indels can
still support a compound NClose. Exactly 100 kb passes the addition threshold.

- `read` uses complete raw molecules, primary/SA discovery and whole-reference
  realignment, then selects each molecule's two reliable outer alignments.
- `olc` uses the custom assembler in `local_olc.py`: canonical minimizer seeds,
  banded overlap alignment, an oriented overlap graph, transitive reduction,
  paths that stop at branches, and an alignment-based consensus. No external
  assembler/overlapper is used. Fully covered reads remain consensus evidence;
  tips and bubbles are not pruned. Unitigs are realigned to the whole reference.
  Both terminal junctions require raw-molecule support; for a longer path this
  need not be the same molecule at both ends. The paired outer-read support is
  reported separately, and assembly input-read count is never junction support.

Native assembly runs default to raw-read-only rescue (`read`), without local
OLC assembly. OLC remains an explicit comparison option; VCF input defaults to
`off`. The workspace `run.sh` accepts `--raw-rescue-method` directly (this takes
precedence over the method in `--option_skype`). Select a method and tune
stage-24 arguments through the native orchestrator:

```bash
bash run.sh HCC1937 skype/HCC1937_read
bash run.sh --raw-rescue-method read HCC1937 skype/HCC1937_read_explicit
bash run.sh --raw-rescue-method olc HCC1937 skype/HCC1937_olc
bash run.sh --raw-rescue-method off HCC1937 skype/HCC1937_baseline
bash run.sh --raw-rescue-method olc --option_skype='--raw-rescue-options "--olc-min-overlap 1500 --olc-min-identity 0.99 --olc-jobs 4"' HCC1937
```

`--skype-start-at 24` on `run.sh` starts from a prepared stage-23 result. Compare
methods in separate copies of the same baseline. Standalone stage 24 only
prepares an augmented handoff; `pipeline.py` installs it, updates the PAF and
event catalog, and invokes stages 10--23 **once** if compression retained a new
NClose. Stage 31 then emits the final results. `24_raw_rescue/round.json` records
installation/refitting, so restarting an already completed pass cannot add a
second round. A restart reuses the recorded rescue method when none is supplied.
A fresh stage-01 run resets that state.

`24_raw_rescue/` contains the D BED/TSV/PDF, extracted reads, alignment evidence,
OLC graph/unitig details, accepted and rejected candidates, compression matches,
and `before/` snapshots. `summary.json` records the added count and installation
files. New nodes use `raw_rescue_<method>_<candidate>` names in the event catalog.
OLC parameters include seed size/window/copy limit, overlap length/identity,
overhang, indel limit and band width; `--help` lists the stage-24 options.
Query controls include `--query-context-bins`, `--query-flank-bins`,
`--query-model-delta`, and `--query-padding`.
The stage uses the existing SciPy, pysam and numba packages in the SKYPE environment.

## Results

All files are written to the SKYPE output directory.

| Stage | Output | Description |
| --- | --- | --- |
| `01_Preprocess_NClose.py` | `01_nclose_data.pkl` | Exact three-field graph handoff: `contig_data`, `nclose_nodes`, and `telo_contig`. |
| `01_Preprocess_NClose.py` | `skype_options.json` | Options routed to preprocessing and graph search. |
| `01_Preprocess_NClose.py` | `stage01_nclose_summary.json`, `stage01_nclose_rejections.tsv` | Ordered split/filter counts and first-rejection provenance. |
| `01_Preprocess_NClose.py` | `censat_endpoint_summary.json`, `censat_endpoint_candidates.tsv` | Fixed partition counts and independent endpoint consistency decisions. |
| `23_run_nnls.py` | `weight.npy`, `predict_B.npy` | Full-column raw NNLS weights and reconstructed depth from one solve. |
| `24_raw_nclose_rescue.py` | `24_raw_rescue/` | D queries, raw/OLC evidence, compressed additions and one-pass provenance. |
| `31_depth_analysis.py` | `SV_call_result.vcf` | Native `BND`, `INV`, `DEL`, and `DUP` calls with normalized copy-number support. |
| `31_depth_analysis.py` | `SKYPE_result.bed` | Simplified native breakend, indel, centromere-fragment, amplicon, and virtual-inversion calls. |
| `31_depth_analysis.py` | `nclose_report.tsv` | Per-NClose coordinates, orientation, copy-number support, and NClose-preprocessing exclusion provenance. |
| `31_depth_analysis.py` | `total_cov.png`, `total_cov.pdf` | Circos view of observed and reconstructed depth with the reported variants. |
| `31_depth_analysis.py` | `SV_benchmark_result.vcf` | Copy-number-annotated input VCF, produced instead of the native VCF/BED in VCF-input runs. |

## Running SKYPE

Run SKYPE through the **[ACCtools pipeline](https://github.com/ACCtools/ACCtools-pipeline)**,
which prepares the assembly, alignments, depth data, reference resources, and
stage arguments.

Complete-assembly input remains a separate workflow. ACCtools invokes
`full_assembly_pipeline.py`, which treats mapped FASTA records as matrix
features and retains its Virtual SKY, karyotype, and coverage outputs. The
native numbered pipeline is not entered for that input type.
