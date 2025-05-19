# Graph-based cancer genome scaffolding tool

**SKYPE** generates following results.
- **Circos plot**, which shows copy number of chromosomes
- **Virtual SKY diagram**, which shows most reasonable combinations of mutated chromosome
</br></br>

**SKYPE** implements following steps to achieve above results.
- **Classify contigs** as their mutational role and prune unnecessary parts
- **Compress contigs** by their most significant terminal node (NClose)
- **Build breakend graph** by connecting NClose nodes and telomere nodes
- Find all simple **telomere-to-telomere paths of breakend graph**
- **Recover full path** with normal contigs
- Optimize copy number of each paths by using customized **SI-NNLS** algorithm
- Visualize results

## Prequisites
It is recommended to use **[ACCtools pipeline](https://github.com/ACCtools/ACCtools-pipeline)**, rather than using SKYPE alone.

Results are saved as 2 png files.
- \<SKYPE result folder\>/\<CELL_LINE_RESULT_FOLDER\>/total_cov.png : Circos plot
- \<SKYPE result folder\>/\<CELL_LINE_RESULT_FOLDER\>/virtual_sky.png : Virtual SKY diagram

## Arguments and Variables
First argument defines cell line to analyze.
Second argument(optional) defines the folder to save SKYPE result.
Change THREAD variable 

## Example
Please check `README.md` in [ACCtools pipeline](https://github.com/ACCtools/ACCtools-pipeline)
