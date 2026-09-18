# SKYPE

[![Tests](https://github.com/ACCtools/SKYPE/actions/workflows/tests.yml/badge.svg)](https://github.com/ACCtools/SKYPE/actions/workflows/tests.yml)

**SKYPE** is an assembly-based structural-variant caller that uses raw-read
depth and junction evidence to support its calls. It builds candidate chromosome
paths from assembly junctions, fits their depth contributions, and reports
variants with normalized support. Native runs can also rescue missing junction
candidates from raw reads in regions with unexplained depth changes.

SKYPE also supports depth annotation of an existing VCF and a separate workflow
for complete genome assemblies.

[Usage and options](docs/usage.md) · [Methods](docs/methods.md) ·
[Output formats and weight interpretation](docs/outputs.md)

## Quickstart

Run SKYPE through [ACCtools-pipeline](https://github.com/ACCtools/ACCtools-pipeline),
which prepares assemblies, alignments, depth and reference resources. First set
up the environment in its [installation guide](https://github.com/ACCtools/ACCtools-pipeline#dependency).
Then, from the ACCtools-pipeline checkout:

```bash
mamba activate skype

# PacBio HiFi reads; options precede the working directory.
python SKYPE.py run_hifi --reference hg38 -p SAMPLE -t 16 \
  /path/to/work /path/to/reads.fastq.gz

# Reuse an assembly and a PanDepth result.
python SKYPE.py analysis --reference hg38 -p SAMPLE -t 16 \
  /path/to/work /path/to/contigs.fa /path/to/unitigs.fa \
  /path/to/depth.win.stat.gz
```

Use one reference build for the assembly alignments, depth and any input VCF.
The reference choices are `hs1` and `hg38`; the default is `hs1`.
See [usage](docs/usage.md) for input requirements, ONT/Flye commands, VCF input,
raw-read rescue options and partial restarts. The research workspace's `run.sh`
is a separate convenience wrapper and is not included in this repository.

## Input workflows

| Workflow | Input and purpose | Main result |
| --- | --- | --- |
| Native assembly | Discover NCloses from assembly alignments, fit depth and optionally rescue candidates from raw reads. | `SV_call_result.vcf` |
| VCF input | Evaluate the supplied variants using SKYPE depth support; select with `--benchmark_vcf_loc`. | `SV_benchmark_result.vcf` |
| Full assembly | Use each supplied FASTA record as a matrix path; select with `--full_assembly`. | Separate coverage, Virtual SKY and karyotype outputs |

Native assembly runs default to raw-read rescue (`read`); VCF input defaults to
rescue `off`. Local OLC rescue is an explicit option. Full-assembly input uses
its own workflow. See [input workflows](docs/usage.md#input-workflows).

## NCloses and structures

An **NClose** is SKYPE's tracked assembly-junction unit. It includes Type 1,
Type 2 and Type 4 events. An NClose can retain several internal alignment
pieces; it is not necessarily a decomposition into primitive junctions.

A **structure** carries a weight and a count of how often it uses each original
NClose. Structures include chromosome paths, independent Type 4 features,
MERGE_TYPE4, AMP and VIRTUAL_INV. Several structures may use the same NClose.
See the [definitions and pipeline](docs/methods.md#ncloses-and-structures).

## Main results

Files are written to the selected SKYPE output directory.

| File | What to use it for |
| --- | --- |
| `SV_call_result.vcf` | Native variant calls: constituent BNDs for compound structures, and symbolic DEL/DUP for ordinary Type 4. |
| `SKYPE_result.bed` | NClose and structure summary views, with their weight scopes identified. |
| `nclose_report.tsv` | Original NClose totals and separate preprocessing history. |
| `structure_report.tsv` | Each native structure's own weight and provenance. |
| `structure_nclose_usage.tsv` | NClose occurrence counts and each structure's contribution. |
| `total_cov.png`, `total_cov.pdf` | Observed and fitted depth with NClose and structure links. |
| `SV_benchmark_result.vcf` | VCF-input records annotated with `SKYPE_CN` and `SKYPE_STATUS`. |

The native structure reports do not replace the VCF-input workflow's existing
reporting rules. Detailed schemas, IDs and intermediate files are documented in
[outputs](docs/outputs.md).

## Interpreting weights

Native outputs share this calculation:

```text
NClose CN = sum(structure weight × NClose occurrence count) / N
N = median read depth excluding chrM / 2
```

For example, a path with normalized weight `3` and a virtual inversion with
normalized weight `2` each using NClose A once give A a total of `5`. A virtual
inversion contributes its weight to each constituent NClose according to its
usage count; the weight is not divided between its two NCloses.

These values are support estimates, not read VAFs. Structure CN, NClose CN and
the predicted depth of a reference segment describe different quantities.
Contributions are summed before the native output threshold (`CN > 0.1`).
VCF geometry merging and the existing PATH_SPLIT representation can make its
record count differ from the original NClose count.

See [weight accounting](docs/outputs.md#weight-accounting) for AMP normalization,
virtual contributions, Type 4 path contributions and worked examples, and
[output representations](docs/outputs.md#output-representations) for PATH_SPLIT.

## Further documentation

- [Usage](docs/usage.md): setup, inputs, options, rescue, caches and restarting.
- [Methods](docs/methods.md): NClose processing, graph paths, depth fitting and raw-read rescue.
- [Outputs](docs/outputs.md): weight semantics, file schemas, VCF fields and result regeneration.
