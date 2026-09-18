# Usage

[README](../README.md) · [Methods](methods.md) · [Outputs](outputs.md)

## Setup and command locations

Use [ACCtools-pipeline](https://github.com/ACCtools/ACCtools-pipeline) to prepare
assemblies, reference resources, alignments and depth. Follow its
[environment instructions](https://github.com/ACCtools/ACCtools-pipeline#dependency)
for the Python environment and external tools. Keep the caller and ACCtools
checkouts compatible when using features from this documentation.

```bash
git clone https://github.com/ACCtools/ACCtools-pipeline.git
cd ACCtools-pipeline
mamba activate skype
python SKYPE.py --help
```

The commands below run from that ACCtools-pipeline directory. `/path/to/...`
denotes a path to replace with your own input or output. `--dependency_loc`
selects an existing dependency directory; ACCtools installs dependencies when
no directory is supplied. That directory contains the SKYPE checkout along
with the other tools and reference resources.

The research workspace's `bash run.sh` is a convenience wrapper around
ACCtools. It is not part of the SKYPE repository, and its sample names and data
locations are specific to that workspace. Public examples here use the ACCtools
CLI directly.

## Input workflows

### Native assembly from reads

For PacBio HiFi:

```bash
python SKYPE.py run_hifi --reference hg38 -p SAMPLE -t 16 \
  /path/to/work /path/to/reads.fastq.gz
```

Place all options before the working directory for `run_hifi`: every argument
after that directory is interpreted as an input read file. Multiple read files
can follow it.

Other assembly inputs supported by ACCtools include:

```bash
# ONT R10/HQ with hifiasm
python SKYPE.py run_hifi --reference hg38 -t 16 \
  --hifiasm_args="--ont --chem-c 0" /path/to/work /path/to/ont.fastq.gz

# PacBio CLR with Flye
python SKYPE.py run_flye --reference hg38 -t 16 \
  /path/to/work pacbio-raw /path/to/clr.fastq.gz

# ONT with Flye
python SKYPE.py run_flye --reference hg38 -t 16 \
  /path/to/work nano-raw /path/to/ont.fastq.gz
```

### Existing assembly and depth

```bash
python SKYPE.py analysis --reference hg38 -p SAMPLE -t 16 \
  --skype_dir /path/to/results \
  /path/to/work /path/to/contigs.fa /path/to/unitigs.fa \
  /path/to/depth.win.stat.gz
```

This reuses the supplied assembly and PanDepth result while ACCtools prepares
the alignment and reference inputs for SKYPE. Use a working directory containing
the corresponding ACCtools read-mapping inputs: native BAM validation and
raw-read rescue require the sample's mapped reads. Supplying a depth file alone
does not provide raw-read junction evidence.

### VCF input

```bash
python SKYPE.py analysis --reference hg38 -p SAMPLE -t 16 \
  --benchmark_vcf_loc /path/to/input.vcf \
  --skype_dir /path/to/vcf_results \
  /path/to/work /path/to/contigs.fa /path/to/unitigs.fa \
  /path/to/depth.win.stat.gz
```

`--benchmark_vcf_loc` supplies the variant candidates instead of discovering
NCloses from assembly alignments. Assembly alignment still provides
telomere/neotelomere anchors, and read depth supplies the support estimates.
The option is also available on `run_hifi` and `run_flye`.

Only `FILTER=PASS` and `FILTER=.` records are evaluated by default. Replace that
set, if needed, with `--option_skype="--vcf-filter-pass PASS . Candidate"`.
Matching is exact and case-sensitive. Parsing, skipped records and orientation
diagnostics are described in [VCF-input outputs](outputs.md#vcf-input-outputs).

### Full-assembly input

```bash
python SKYPE.py analysis --reference hs1 -p SAMPLE -t 16 \
  --full_assembly /path/to/complete.fa \
  --skype_dir /path/to/full_assembly_results \
  /path/to/work /path/to/complete.fa /path/to/complete.fa \
  /path/to/depth.win.stat.gz
```

The contig/unitig positional slots are ignored in this workflow; the example
fills them with the complete FASTA. Each FASTA record becomes a matrix path.
ACCtools aligns the FASTA with minimap2 and alignasm and invokes
[full_assembly_pipeline.py](../full_assembly_pipeline.py), preserving its
Virtual SKY, karyotype and coverage outputs. It does not run the native numbered
pipeline or use the native structure-output reporting contract.

`--full_assembly` and `--benchmark_vcf_loc` are mutually exclusive. Full assembly
is also accepted by the read-based ACCtools commands, which still map reads and
calculate sample depth but skip hifiasm/Flye assembly.

### Reference and output directory

`--reference` accepts `hs1` or `hg38` and defaults to `hs1`. Assembly alignments,
depth, input VCF coordinates and reference annotations must use the same build.
For native/VCF runs, the default SKYPE result directory is `30_skype/` under the
working directory for `hs1`, or `31_skype_hg38/` for `hg38`.
The `analysis` command's `--skype_dir` selects another result directory.

## Native options

Forward preprocessing and graph settings with `--option_skype` (the alias
`--option-skype` also works):

```bash
python SKYPE.py analysis --reference hg38 \
  --option_skype="--nclose-min-ref-span 10000 --add-indel-graph" \
  /path/to/work /path/to/contigs.fa /path/to/unitigs.fa \
  /path/to/depth.win.stat.gz
```

Common settings are:

| Option inside `--option_skype` | Default | Meaning |
| --- | --- | --- |
| `--nclose-min-ref-span` | `1000` bp | Stage-01 anchor-span cutoff; `0` disables it. Independently validated CEN-SAT candidates are exempt. |
| `--exclude-nclose-list-loc` | Unset | Stage-01 user exclusion list. |
| `--check-nclose-count` | Disabled | Enable the stage-01 raw-read junction VAF filter. |
| `--nclose-count-vaf-threshold` | `0.1` | Threshold for that optional VAF filter. |
| `--skip-bam-analysis` | Disabled | Skip BAM validation of translocation candidates and the associated virtual-inversion candidate removal. |
| `--add-indel-graph` | Disabled | Add selected depth-supported Type 4 edges without increasing graph dimensions; VCF INS events are excluded. |
| `--vcf-filter-pass` | `PASS .` | Accepted input VCF filter values; ignored outside VCF mode. |
| `--limit-combinations` | Automatic | File specifying the graph-search limit pair. |
| `--verbose` | Disabled | Graph-search diagnostics. |

Do not combine `--skip-bam-analysis` and `--check-nclose-count`, since the latter
needs BAM analysis. Use ACCtools' top-level `-t`, `-d` and `--progress` options
for execution settings. Input PAFs, `--vcf-input` and reference/depth paths are
owned by the orchestrator and are not accepted as native extra options.

Stage 01 saves routed stage-01/10 settings in `skype_options.json`. A fresh
stage-01 run replaces those settings. Stage 10 normally reuses them; an explicit
nonempty `--option_skype` on a stage-10 restart replaces its graph options for
that invocation. Direct stage-10 CLI values take precedence over saved values;
a direct `--option_skype=""` selects defaults. Preprocessing changes require a
stage-01 restart. [skype_options.py](../skype_options.py) defines the accepted
options and their aliases.

## Raw-read rescue

Native assembly defaults to `read`; VCF input defaults to `off`. Select a method
through ACCtools' extra-option string:

```bash
python SKYPE.py analysis --reference hg38 \
  --option_skype="--raw-rescue-method off" \
  /path/to/work /path/to/contigs.fa /path/to/unitigs.fa \
  /path/to/depth.win.stat.gz
```

Use `read` for raw-molecule rescue or `olc` for the explicit local-assembly
comparison. Nested rescue settings use a quoted argument string:

```bash
python SKYPE.py analysis --reference hg38 \
  --option_skype='--raw-rescue-method olc --raw-rescue-options "--olc-min-overlap 1500 --olc-min-identity 0.99 --olc-jobs 4"' \
  /path/to/work /path/to/contigs.fa /path/to/unitigs.fa \
  /path/to/depth.win.stat.gz
```

Query controls include `--query-mode adaptive|fixed`, `--query-context-bins`,
`--query-flank-bins`, `--query-model-delta` and `--query-padding`.
OLC controls include seed/window/copy limits, overlap length and identity,
overhang, band, indel limit and job count. The complete CLI is available with
`python /path/to/deps/SKYPE/24_raw_nclose_rescue.py --help`.
See [rescue methods](methods.md#one-raw-read-rescue-pass) for default gates and
the one-refit behavior.

## Restarts

Native restart stages are `0, 1, 10, 11, 21, 22, 23, 24, 31`. Use the same
sample, reference, inputs and result directory as the original run. A nonzero
restart checks the required artifacts from skipped stages.

To regenerate reports and plots without refitting:

```bash
python SKYPE.py analysis --reference hg38 -p SAMPLE -t 16 \
  --skype_dir /path/to/existing_results --skype_force --skype_start_at 31 \
  /path/to/work /path/to/contigs.fa /path/to/unitigs.fa \
  /path/to/depth.win.stat.gz
```

For VCF-input results, also pass the original `--benchmark_vcf_loc`.
`--print_args` prints the SKYPE subprocess commands for inspection.

| Change or task | Restart point |
| --- | --- |
| Stage-01 options, anchor cutoff or compound-candidate construction | `1` |
| Graph options with prepared preprocessing inputs | `10` |
| Candidate matrix/column membership rebuild | `22`, with current upstream depth artifacts |
| Raw-read rescue from a prepared stage-23 baseline | `24` |
| Native output accounting/report regeneration with a completed fit | `31` |

Stage 31 reconstructs a missing `structure_nclose_model.pkl` from saved path,
catalog and circuit data. Stage 22 rebuilds the table with the matrix columns.
Output regeneration does not alter the fitted coefficients or predicted depth.
See [result regeneration](outputs.md#regenerating-existing-results) for output
format changes and required provenance.

To compare rescue methods, use separate copies of the same pre-rescue baseline.
Standalone stage 24 prepares augmented inputs; the orchestrator installs them
and refits once when new NCloses survive compression. A completed
`24_raw_rescue/round.json` prevents another rescue round on restart. An omitted
method reuses the recorded rescue method; a fresh stage-01 run resets the state.
Completed rescue snapshots are not automatically migrated when candidate rules
change; evaluate such changes from a pre-rescue baseline.

Full-assembly restarts belong to its separate workflow and accept stages
`21, 22, 23, 30, 31`; a fresh default run begins at its stage 21.

## Caches

The CEN-SAT preparation directory contains the partition, candidate metadata,
source-coordinate manifest, extracted FASTA and realignment PAF. It is reused
when assembly/PAF/reference/BED signatures, workflow implementation, minimap2
version and options match. `--skype_force` alone reuses this alignment cache;
ACCtools alignment `--preprocess_force` or the direct native CLI's
`--alignment-force` requests its rebuild.

Reference-based minimap2 alignments share indexes under
`<dependency_loc>/reference_indexes/`. Indexes are built on demand and keyed by
reference FASTA signature, preset and minimap2 version. The native orchestrator,
CEN-SAT CLI and rescue CLI accept `--reference-index-cache`; an explicit rescue
`--reference-index` inside `--raw-rescue-options` takes precedence.

## Direct native invocation

With all inputs and resources already prepared, inspect:

```bash
python /path/to/deps/SKYPE/pipeline.py --help
```

This entry point requires the PAF inputs, depth, dependency directory and
reference FASTA/FAI/BED resources normally supplied by ACCtools. It accepts
`--raw-rescue-method`, `--raw-rescue-options` and `--reference-index-cache`
directly. Flag spellings differ between the ACCtools and direct native CLIs;
use the help for the entry point you are invoking.
