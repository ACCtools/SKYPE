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
| `24_cluster_weight.py` | `cluster_raw_zero_path_decisions.tsv` | High-quality, non-CENSAT MAPQ-60 zero-read NCloses and any carrier paths removed before compositional A/B/AB filtering and refitting. |
| `25_cluster_nclose_read_count.py` | `cluster_nclose_read_counts.tsv` | Raw-read junction support plus strand-aware 5 kb reference support sharing the actual junction anchors for BND NCloses and graph-only indels carried by displayed cluster paths. |
| `26_nclose_read_depth_plot.py` | `nclose_read_depth_comparison.tsv`, `nclose_read_depth_groups.pdf` | Compares fitted cluster depth with strand-aware one-sided 500 kb local-depth/VAF-derived raw-read depth and groups NCloses by co-carriage in displayed paths. |
| `26_nclose_read_depth_plot.py` | `nclose_read_depth_scatter.png`, `nclose_read_depth_scatter.pdf` | Direct `y=x` comparison of VAF-derived and NClose-derived raw depth after excluding every NClose with either endpoint overlapping CENSAT; remaining quality and agreement problems use distinct markers. |
| `30_virtual_sky.py` | `virtual_sky.png`, `virtual_sky.pdf` | Virtual SKY diagram of reconstructed chromosome paths, colored by chromosome of origin and labeled with normalized copy number. |
| `30_virtual_sky.py` | `virtual_sky_cluster.png`, `virtual_sky_cluster.pdf` | Cluster Virtual SKY with path weights shown as fitted raw depth (`x`) and each countable junction shown as local-coverage/VAF-derived depth (`t(...) : depth_x`); uncountable junctions are labeled `NA`. |
| `30_virtual_sky.py` | `virtual_sky_cluster_refspan_C*.png`, `virtual_sky_cluster_refspan_C*.pdf` | One chromosome-coordinate branch graph per non-singleton ambiguity component formed by at least 100 kb of continuous shared breakend-graph reference span. Each chromosome uses an independent Mb axis; route width and labels show fitted path contributions. No-NClose paths covering at least 90% of one reference chromosome are excluded. |
| `30_virtual_sky.py` | `virtual_sky_cluster_refspan_lane_C*`, `virtual_sky_cluster_refspan_profile_C*`, `virtual_sky_cluster_refspan_attrib_C*` | The same components with the cell line's own reference depth. `lane` puts a copy-number ribbon inside each chromosome lane, `profile` adds full copy-number axes per chromosome, and `attrib` turns every shared block and member footprint into observed-versus-contributed bars. Member depth is the path's own NNLS matrix column scaled by its fitted weight, so it is read at exactly the windows the solver fitted. Select layouts with `SKYPE_REFSPAN_MODES` (default `plain,lane,profile,attribution`). |
| `30_virtual_sky.py` | `reference_path_clusters.tsv`, `reference_path_cluster_overlaps.tsv` | Auditable path membership/full-reference exclusions and exact half-open reference intervals supporting every cluster edge. |
| `30_virtual_sky.py` | `reference_path_cluster_depth.tsv` | Every number behind the depth figures: observed, fitted (`A·w`), cluster-summed, other-path, and per-member copy number over each shared block and each member's whole reference footprint. |
| `30_virtual_sky.py` | `karyotype.txt` | Tab-separated ISCN-style karyotype labels and their normalized depth (`N`). |
| `31_depth_analysis.py` | `total_cov.png`, `total_cov.pdf` | Circos view of observed and reconstructed chromosome depth, with copy-number-weighted breakends, inversions, indels, amplicons, and neotelomeres. |
| `31_depth_analysis.py` | `SV_call_result.vcf` | Structural-variant calls (`BND`, `INV`, `DEL`, and `DUP`) with normalized copy-number support. Produced in assembly-input modes. |
| `31_depth_analysis.py` | `SKYPE_result.bed` | Simplified BED report of breakends, deletions, duplications, centromere fragments, amplicons, and virtual inversions. Produced in assembly-input modes. |
| `31_depth_analysis.py` | `SV_benchmark_result.vcf` | Copy-number-annotated copy of the input VCF. Produced only in VCF input mode instead of `SV_call_result.vcf` and `SKYPE_result.bed`. |

The unsuffixed result files are always generated. Karyotype mode also generates `_filter` and `_cluster` versions of the Virtual SKY, karyotype, Circos, VCF, and BED results. Variant mode and VCF input mode generate only the unsuffixed versions.

## Running SKYPE

Run SKYPE through the **[ACCtools pipeline](https://github.com/ACCtools/ACCtools-pipeline)**, which prepares the assembly, alignments, depth data, reference resources, and stage arguments required by SKYPE.

For complete-assembly input, ACCtools runs minimap2 followed by alignasm and
invokes `full_assembly_pipeline.py` once with the resulting `*.aln.paf`. That
entrypoint owns depth splitting, matrix construction, NNLS, Virtual SKY, and
coverage rendering; the numbered stage scripts are not entered. The normal and
full-assembly routes call the same reusable stage-30 and stage-31 renderers, so
their figure layout and styling stay consistent.
