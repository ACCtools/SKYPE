# Graph-based cancer genome karyotyping and structural-variant analysis tool

[![Tests](https://github.com/ACCtools/SKYPE/actions/workflows/tests.yml/badge.svg)](https://github.com/ACCtools/SKYPE/actions/workflows/tests.yml)

**SKYPE** generates the following results.

- A **Circos plot** showing observed and reconstructed copy number together with rearrangement junctions, indels, amplicons, and neotelomeres
- A **Virtual SKY diagram** showing depth-weighted combinations of reconstructed cancer chromosome paths
- An **ISCN-style karyotype summary** and copy-number-weighted structural-variant results

**SKYPE** implements the following steps to generate these results.

- Derive rearrangement junctions from contig/unitig alignments, or import them from a supported VCF
- **Classify contig segments** by their mutational role and prune noisy or unnecessary segments
- **Compress contigs** around their most significant terminal nodes (NClose)
- **Build a breakend graph** by connecting NClose and telomere nodes
- Find candidate **telomere-to-telomere paths** and circular components in the breakend graph
- **Recover full chromosome paths** with normal reference segments
- Estimate the copy number of each path using non-negative least squares
- Visualize the reconstructed karyotype and report copy-number-weighted variants

## Results

All files below are written to the SKYPE output directory.

| Stage | Output | Description |
| --- | --- | --- |
| `30_virtual_sky.py` | `virtual_sky.png`, `virtual_sky.pdf` | Virtual SKY diagram of reconstructed chromosome paths, colored by chromosome of origin and labeled with normalized copy number. |
| `30_virtual_sky.py` | `karyotype.txt` | Tab-separated ISCN-style karyotype labels and their normalized depth (`N`). |
| `31_depth_analysis.py` | `total_cov.png`, `total_cov.pdf` | Circos view of observed and reconstructed chromosome depth, with copy-number-weighted breakends, inversions, indels, amplicons, and neotelomeres. |
| `31_depth_analysis.py` | `SV_call_result.vcf` | Structural-variant calls (`BND`, `INV`, `DEL`, and `DUP`) with normalized copy-number support. Produced in assembly-input modes. |
| `31_depth_analysis.py` | `SKYPE_result.bed` | Simplified BED report of breakends, deletions, duplications, centromere fragments, amplicons, and virtual inversions. Produced in assembly-input modes. |
| `31_depth_analysis.py` | `SV_benchmark_result.vcf` | Copy-number-annotated copy of the input VCF. Produced only in VCF input mode instead of `SV_call_result.vcf` and `SKYPE_result.bed`. |

The unsuffixed result files are always generated. Karyotype mode also generates `_filter` and `_cluster` versions of the Virtual SKY, karyotype, Circos, VCF, and BED results. Variant mode and VCF input mode generate only the unsuffixed versions.

## Running SKYPE

Run SKYPE through the **[ACCtools pipeline](https://github.com/ACCtools/ACCtools-pipeline)**, which prepares the assembly, alignments, depth data, reference resources, and stage arguments required by SKYPE.
