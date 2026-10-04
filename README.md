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
The native VCF alone expands each eligible NClose into BNDs between adjacent
alignments. CEN-SAT, telomere-derived and virtual-inversion NCloses, as well as
read and OLC rescue, are included. Normal reference continuations are omitted.
The depth model, BED and plots retain original NCloses.
Native accounting reuses an existing NClose ID when another source pair has
exactly the same two BND boundaries and retained sides. Source pairs and unitigs
remain traceable in `nclose_sources.tsv`; nearby coordinates are not merged.

A **structure** carries a weight and a count of how often it uses each original
NClose. Structures include chromosome paths, independent Type 4 features,
MERGE_TYPE4, AMP and VIRTUAL_INV. Several structures may use the same NClose.
See the [definitions and pipeline](docs/methods.md#ncloses-and-structures).

## Main results

Files are written to the selected SKYPE output directory.

| File | What to use it for |
| --- | --- |
| `SV_call_result.vcf` | Native variant calls: adjacent-alignment BNDs within eligible constituent NCloses, and symbolic DEL/DUP for ordinary Type 4. |
| `SV_call_result.nclose_bnds.tsv` | Each original NClose source's junctions and their exact merged VCF BND IDs. |
| `SV_call_result.bnd_weights.tsv` | Each structure/source NClose's contribution to an exported BND. |
| `SKYPE_result.bed` | NClose and structure summary views, with their weight scopes identified. |
| `nclose_report.tsv` | NClose totals after exact BND identity reuse and representative preprocessing history. |
| `nclose_sources.tsv` | Original node pairs/unitigs and their shared NClose identities. |
| `structure_report.tsv` | Each native structure's own weight and provenance. |
| `structure_nclose_usage.tsv` | NClose occurrence counts and each structure's contribution. |
| `terminal_evidence.json` | Physical terminal contexts and explicitly assembly-derived sequence evidence, with source/reference verification states. |
| `structure_terminal_usage.tsv` | Every native PATH start and end, including zero weights, joined to its evidence context. |
| `terminal_host_inventory.tsv` | All graph terminal hosts, including hosts with no modeled PATH use. |
| `total_cov.png`, `total_cov.pdf` | Observed and fitted depth with NClose and structure links. |
| `SV_benchmark_result.vcf` | VCF-input records annotated with `SKYPE_CN` and `SKYPE_STATUS`. |

The native structure reports do not replace the VCF-input workflow's existing
reporting rules. Detailed schemas, IDs and intermediate files are documented in
[outputs](docs/outputs.md).

## Interpreting weights

Native outputs share this calculation:

```text
NClose CN = sum(structure weight × NClose occurrence count) / N
BND CN = sum(structure weight × source NClose occurrence count × BND occurrences in that source) / N
N = median read depth excluding chrM / 2
```

For example, a path with normalized weight `3` and a virtual inversion with
normalized weight `2` each using NClose A once give A a total of `5`. A virtual
inversion contributes its weight to each constituent NClose according to its
usage count; the weight is not divided between its two NCloses.
Likewise, splitting a source NClose gives each distinct junction its full source
contribution. Junctions merge only when both chromosomes, boundaries and retained
sides match exactly, including reverse-complement descriptions. Source-specific
counts keep different internal chains separate even if their outer endpoints match.

These values are support estimates, not read VAFs. Structure CN, NClose CN and
the predicted depth of a reference segment describe different quantities.
Contributions are summed before the native output threshold (`CN > 0.1`).
VCF geometry merging and reciprocal BND mate records can make its
record count differ from the original NClose count.

See [weight accounting](docs/outputs.md#weight-accounting) for AMP normalization,
virtual contributions, Type 4 path contributions and worked examples, and
[output representations](docs/outputs.md#output-representations) for original
NClose endpoints and compound summaries.

## Native terminal evidence

Stage 21 records its actual resolved terminal traversal in
`terminal_component_context.json`, including a host whose depth row was trimmed
to zero. Input hashes are captured before resolution and checked again before
publication; rebuilding stage 21 invalidates any older capture first.
Stage 31 joins this record to each PATH end and writes separate evidence
reports. Graph labels and unused direct edges do not substitute for the physical
port used in a path. Both occurrences of a terminal node in one path are retained.
The annotation does not change graph construction, PAFs, depth, weights, or calls.
It applies to the native PATH workflow; it does not annotate the separate
full-assembly-only or VCF-only workflows.

The four display classes are ordinary graph anchor, observed nonreference
telomeric extension, donor-route repeat, and unavailable sequence/ambiguous port.
Always read the separate physical-port, `assembly_evidence_state`, source-repeat,
reference-comparison, anchor-quality and internal-repeat/farther-flank fields.
An ordinary graph anchor is a model label, not a measured absence of a telomere.
A repeat extension can still be internal. An assembly boundary does not establish
a chromosome cap, healing event, somatic event, or complete chromosome path.
The usage table's conditional fitted contributions do not measure identified
telomere dosage or chromosome-end counts.

`terminal_evidence_rules.json` fixes the uniform motif and context rules.
Full tract span, canonical-covered bases, aligned-side overlap, and outward-only
coverage are distinct quantities. Strong extension descriptions require the
primary span/canonical minima on the outward partition and a complete callable
reference comparison. A clipped reference-end window cannot supply a negative
comparison. No independent-molecule evidence is assumed by these reports.

The pipeline passes its effective assembly FASTA, assembly-alignment reference,
and raw PAF to the report. Input FASTAs are read without creating indexes.
Existing runs can be annotated with `python terminal_evidence.py PREFIX` and
optional `--source-fasta`, `--reference-fasta`, `--raw-paf`, and `--source-binding`.
Missing stage-21 capture or source inputs remain explicit unknown states; this
command does not rebuild the graph, alignment, or fit. `terminal_source_inputs.json`
records requested files only and never certifies historical generation.

Historical source verification requires an upstream
`<raw_paf>.source_binding.json` (or explicit binding path) with schema
`SKYPE.assembly_alignment_source.v1`, complete fresh-generation status,
`producer_stage=raw_and_alternate_alignment_generation`, equal `inputs_before`
and `inputs_after` FASTA/reference path-SHA pairs, primary/alternate output
path-SHA pairs, matching `outputs_at_generation` recorded immediately after each
mapper finishes, and recorded generation argv. Its `reference_index_binding` must
use `SKYPE.reference_index_source.v1`, bind the exact reference, preset and
minimap2 to the index, and link the actual `-d` generation output to the published
index by equal content SHA. The gap-extraction step, when present, retains its
script identity and source/raw/query-file argument linkage. The existing
`.aln.paf.alignasm.json` must then link those exact PAF inputs to the selected PAF.
Only the upstream fresh-generation producer may issue this attestation.
Readable legacy sequences are explicitly unverified; a downstream hash or an
old existence/size/mtime cache never upgrades them to bound source evidence.

## Optional local HiFi evidence

For native assembly models, local read support can be reported separately from
conditional fitted contributions. Supply a coordinate-sorted, indexed HiFi BAM
and its uncompressed reference FASTA with `--local-hifi-bam` and
`--local-hifi-reference` inside the existing pipeline `--option_skype` string.
Add `--local-hifi-export` for separate experimental BND VCFs. Declare the known
assembly-input relationship with `--local-hifi-same-input-as-assembly yes`, `no`,
or `unknown` (the default). This declaration is not independent validation.
These options do not change fitting, weights, candidates, or the default VCF/BED.
The separate full-assembly-only and VCF-input workflows are outside this scope.

A completed native result can be annotated without rerunning the pipeline:

```bash
python local_hifi_evidence.py PREFIX \
  --local-hifi-bam sample.sorted.bam \
  --local-hifi-reference reference.fa \
  --local-hifi-same-input-as-assembly yes \
  --local-hifi-export
```

Every exact retained primitive is included before support or coefficient
selection, including candidates without a modeled carrier and unencodable
reference boundaries. Both fixed tolerances (100 bp primary, 500 bp sensitivity)
use MAPQ ≥20, 500 bp query/reference-span anchors, query gaps from −500 to 1,000 bp,
and at least three distinct read names. The small decoder and
`local_hifi_policy.json` define the complete versioned rule. One read name must
identify one physical CCS molecule. The BAM's SQ names/lengths and all available
M5 checksums are compared to the supplied FASTA; missing M5 remains explicitly
unverified. Matching lengths do not certify the historical model/reference
generation relationship. Native-node/reference length disagreements are recorded
per chromosome. Candidate endpoints on those chromosomes remain unassessable
with null support and an explicit VCF exclusion reason; other candidates retain
the same gate. Coordinates are never silently renumbered.

`PREFIX/local_hifi_evidence/` contains separate attempt directories. A successful
attempt publishes `latest.json`; failures retain a failed manifest and do not
replace that pointer. Each attempt records input and implementation hashes,
all-candidate JSON/TSV evidence, supporting read names, per-endpoint exposure
witnesses, source/structure memberships, and nearby/shared-molecule relations.
`--local-hifi-matrix PATH` can supply a retained matrix for depth-design states;
otherwise the current matrix is used when available. Missing design data remain
unavailable rather than zero. Zero design, entirely masked design and no modeled
carrier are distinct from a fitted coefficient of zero.

`LOCAL_HIFI_SUPPORTED` means the fixed local alignment gate passed. It does not
establish biological existence, unique locus origin, insertion sequence identity,
somatic origin, full-chain linkage, identified dosage, absolute CN, or VAF.
Normal, ONT and whole-chain evidence are null here. Conditional flank exposure
is neither a VAF denominator nor a detection-power estimate. Nearby candidates
and candidates sharing reads are retained as representations; they must not be
counted automatically as independent biological events.

Optional `ExperimentalEvidence.100bp.vcf` and `.500bp.vcf` contain supported
reciprocal BND pairs only. Every record has `FILTER=ExperimentalEvidence`, missing
QUAL and no sample/genotype columns. They are distinct from the default calls.
Unencodable candidates remain in the full sidecar with explicit export reasons.
Any benchmark of these files must explicitly opt into that experimental filter
and retain the full candidate denominator and unresolved cases.

## Further documentation

- [Usage](docs/usage.md): setup, inputs, options, rescue, caches and restarting.
- [Methods](docs/methods.md): NClose processing, graph paths, depth fitting and raw-read rescue.
- [Outputs](docs/outputs.md): weight semantics, file schemas, VCF fields and result regeneration.
