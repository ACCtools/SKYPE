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

### Stage-01/10 refactor status

`01_Preprocess_NClose.py` now owns contig/telomere preprocessing, NClose
post-processing, ecDNA discovery, and the exact three-field
`01_nclose_data.pkl` handoff consumed by `10_Graph_Find_Paths.py`. VCF input is
usable through the split stages. In assembly mode the positional primary PAF
is retained only for command/file-name compatibility: stage 01 reads the
`--alt` unitig PAF and the last `--original-paf-loc` value, and writes
unitig-only rows under the primary-named `*.ppc.paf`. It extracts one
path-ordered terminal NClose from each retained type-1/2 unitig and globally
clusters those pairs with the same coordinate/direction rules used by the
legacy stage 02.
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
each candidate. ACCtools runs stages 01 and 10 by default; `--legacy` selects
the unchanged combined stage 02 route.

## Results

All files are written to the SKYPE output directory.

| Stage | Output | Description |
| --- | --- | --- |
| `01_Preprocess_NClose.py` | `01_nclose_data.pkl` | Exact three-field graph handoff: `contig_data`, `nclose_nodes`, and `telo_contig`. |
| `01_Preprocess_NClose.py` | `stage01_nclose_summary.json`, `stage01_nclose_rejections.tsv` | Ordered split/filter counts and first-rejection provenance. |
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
