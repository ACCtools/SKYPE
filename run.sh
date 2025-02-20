CELL_LINE=$1

PAF_LOC=20_acc_pipe/$CELL_LINE.p/$CELL_LINE.p.aln.paf
PAF_ALT_LOC=20_acc_pipe/$CELL_LINE.a/$CELL_LINE.a.aln.paf

# PAF_LOC=$1
# PAF_ALT_LOC=$2
# PREFIX=$3_$(date +"%H_%M_%S")

if [ -n "$2" ]; then
    PREFIX="$2"
else
    PAF_LOC_BASENAME="${PAF_LOC##*/}"
    PREFIX="30_skype_pipe/${PAF_LOC_BASENAME%%.*}_$(date +"%H_%M_%S")"
fi

TEL_BED="public_data/chm13v2.0_telomere.bed"
CHR_FAI="public_data/chm13v2.0.fa.fai"
RPT_BED="public_data/chm13v2.0_repeat.m.bed"
RCS_BED="public_data/chm13v2.0_censat_v2.1.m.bed"
MAIN_STAT_LOC="/home/hyunwoo/51g_cancer_denovo/51_depth_data/$CELL_LINE.win.stat.gz"

python 00_Contig_Preprocessing.py $PAF_LOC $TEL_BED $CHR_FAI $RPT_BED $RCS_BED --alt $PAF_ALT_LOC
python 01_Contig_Build_Graph.py $PAF_LOC".ppc.paf" $TEL_BED $CHR_FAI
python 04_Build_Breakend_Graph_Limited.py $PAF_LOC".ppc.paf" $CHR_FAI $PAF_LOC".ppc.paf.op.graph.txt" $PREFIX
python 12_Breakend_Graph_Build_Paths.py $PAF_LOC".ppc.paf" $PAF_LOC".ppc.paf.op.graph.txt" $PREFIX
# python 20_Fill_Path.py $PAF_LOC $PAF_LOC".ppc.paf" $PREFIX --alt $PAF_ALT_LOC
# python 21_run_depth.py $PREFIX
# python 22_depth_analysis.py $MAIN_STAT_LOC $PREFIX