# Graph-based cancer genome scaffolding tool

## How to run

```bash
PAF_LOC=$1
TEL_BED="public_data/chm13v2.0_telomere.bed"
CHR_FAI="public_data/chm13v2.0.fa.fai"

python 00_Contig_Preprocessing.py $PAF_LOC $TEL_BED $CHR_FAI
python 01_Contig_Build_Graph.py $PAF_LOC".ppc.paf" $TEL_BED $CHR_FAI
python 10_Graph_Find_Paths.py $PAF_LOC".ppc.paf" $PAF_LOC".ppc.paf.op.graph.txt"
```
```
bash run.sh <PAF_FILE_LOCATION>
```

This bash command lines will output every possible paths from each chromosomes' telomere.  
  
  
Paths are stored in folders. Folder's name represents each paths' startpoint and endpoint.