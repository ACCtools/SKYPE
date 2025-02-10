PAF_LOC=$1
PAF_ALT_LOC=$2

TEL_BED="public_data/chm13v2.0_telomere.bed"
CHR_FAI="public_data/chm13v2.0.fa.fai"

#python 00_Contig_Preprocessing_rewrite.py $PAF_LOC $TEL_BED $CHR_FAI
#python 00_Contig_Preprocessing_rewrite.py $PAF_LOC $TEL_BED $CHR_FAI --alt $PAF_ALT_LOC
#python 01_Contig_Build_Graph_rewrite.py $PAF_LOC".ppc.paf" $TEL_BED $CHR_FAI
python 10_Graph_Find_Paths.py $PAF_LOC".ppc.paf" $PAF_LOC".ppc.paf.op.graph.txt"
