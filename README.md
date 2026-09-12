# Graph-based cancer structural-variant caller

[![Tests](https://github.com/ACCtools/SKYPE/actions/workflows/tests.yml/badge.svg)](https://github.com/ACCtools/SKYPE/actions/workflows/tests.yml)

**SKYPE** discovers structural variants from assembly alignments and assigns
copy-number support from observed read depth. It can also import an existing
VCF and annotate its records with SKYPE depth support.

The native pipeline:

- classifies and compresses assembly segments around NClose junctions;
- builds a breakend graph and enumerates candidate chromosome paths;
- constructs one observed-depth matrix for all candidate paths;
- fits every matrix column with one raw non-negative least-squares solve; and
- emits copy-number-weighted variants, an NClose report, and a coverage plot.

There is no karyotype/variant mode switch, normal-chromosome prior, post-NNLS
NClose filtering, or clustering stage.

### Stage 01/10 pipeline

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

ACCtools `SKYPE.py` prepares `<sample>.utg.censat_endpoints/` beside the source
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
each candidate. ACCtools runs stages 01 and 10 for native assembly/VCF input.

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

## Results

All files are written to the SKYPE output directory.

| Stage | Output | Description |
| --- | --- | --- |
| `01_Preprocess_NClose.py` | `01_nclose_data.pkl` | Exact three-field graph handoff: `contig_data`, `nclose_nodes`, and `telo_contig`. |
| `01_Preprocess_NClose.py` | `skype_options.json` | Options routed to preprocessing and graph search. |
| `01_Preprocess_NClose.py` | `stage01_nclose_summary.json`, `stage01_nclose_rejections.tsv` | Ordered split/filter counts and first-rejection provenance. |
| `01_Preprocess_NClose.py` | `censat_endpoint_summary.json`, `censat_endpoint_candidates.tsv` | Fixed partition counts and independent endpoint consistency decisions. |
| `23_run_nnls.py` | `weight.npy`, `predict_B.npy` | Full-column raw NNLS weights and reconstructed depth from one solve. |
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
