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

SKYPE receives 
- **Result of [alignasm](https://github.com/ACCtools/alignasm)** : Aligned unitig(.r.aln.paf) and primary(.p.aln.paf) PAF file 
- **Reference read depth data** (.win.stat.gz)
  
for input.

## How to run

```bash
CELL_LINE=$1
THREAD=64 # Number of threads for multiprocessing

ALIGNASM_RESULT_FOLDER="<Alignasm result folder>" # Alignasm result folder
SKYPE_RESULT_FOLDER="<SKYPE result folder>" # Result folder name for SKYPE

PAF_LOC=$ALIGNASM_RESULT_FOLDER/$CELL_LINE.p/$CELL_LINE.p.aln.paf
PAF_UTG_LOC=$ALIGNASM_RESULT_FOLDER/$CELL_LINE.r/$CELL_LINE.r.aln.paf

if [ -n "$2" ]; then
    PREFIX="$2"
else
    PAF_LOC_BASENAME="${PAF_LOC##*/}"
    PREFIX="$SKYPE_RESULT_FOLDER/${PAF_LOC_BASENAME%%.*}_$(date +"%H_%M_%S")"
fi

TEL_BED="public_data/chm13v2.0_telomere.bed"
CHR_FAI="public_data/chm13v2.0.fa.fai"
RPT_BED="public_data/chm13v2.0_repeat.m.bed"
RCS_BED="public_data/chm13v2.0_censat_v2.1.m.bed"
CYT_BED="public_data/chm13v2.0_cytobands_allchrs.bed"

MAIN_STAT_LOC="<Reference read depth folder>/$CELL_LINE.win.stat.gz" # File location of reference read depth

PYTHON="python" # Execution option

PROGRESS="" # Use "--progress" to disable tqdm progress bar

$PYTHON 00_Contig_Preprocessing.py $PAF_LOC $TEL_BED $CHR_FAI $RPT_BED $RCS_BED $MAIN_STAT_LOC $PREFIX --alt $PAF_UTG_LOC

$PYTHON 02_Build_Breakend_Graph_Limited.py $PAF_LOC".ppc.paf" $CHR_FAI $RCS_BED $PREFIX \
--orignal_paf_loc 20_acc_pipe/$CELL_LINE.p/$CELL_LINE.p.paf 20_acc_pipe/$CELL_LINE.r/$CELL_LINE.r.paf -t $THREAD $PROGRESS

$PYTHON 11_Ref_Outlier_Contig_Modify.py $PAF_LOC $CHR_FAI $PAF_LOC".ppc.paf" $PREFIX --alt $PAF_UTG_LOC

$PYTHON 21_run_depth.py $PAF_LOC $PAF_LOC".ppc.paf" $PREFIX --alt $PAF_UTG_LOC -t $THREAD $PROGRESS

$PYTHON 22_save_matrix.py $RCS_BED $PAF_LOC".ppc.paf" $MAIN_STAT_LOC $TEL_BED $CHR_FAI $CYT_BED $PREFIX -t $THREAD $PROGRESS

$PYTHON -X juliacall-threads=$THREAD -X juliacall-handle-signals=yes \
23_run_nnls.py $PREFIX -t $THREAD

$PYTHON 30_depth_analysis.py $RCS_BED $PAF_LOC".ppc.paf" $MAIN_STAT_LOC $TEL_BED $CHR_FAI $CYT_BED $PREFIX -t $THREAD $PROGRESS

$PYTHON 31_virtual_sky.py $PAF_LOC".ppc.paf" $MAIN_STAT_LOC $TEL_BED $CHR_FAI $PREFIX $CELL_LINE
```
```
bash run.sh <CELL_LINE_NAME>

bash run.sh <CELL_LINE_NAME> <CELL_LINE_RESULT_FOLDER>
```


Results are saved as 2 png files.
- \<SKYPE result folder\>/\<CELL_LINE_RESULT_FOLDER\>/total_cov.png : Circos plot
- \<SKYPE result folder\>/\<CELL_LINE_RESULT_FOLDER\>/virtual_sky.png : Virtual SKY diagram

## Arguments and Variables
First argument defines cell line to analyze.
Second argument(optional) defines the folder to save SKYPE result.
Change THREAD variable 

## Example
```
bash run.sh U2OS_telo 
bash run.sh PC-3 mypipe/myfolder/
```

