#bsub -W 2:00 -M 100G -o combined.out -e combined.err ./combined.sh
module load BEDTools
FASTA="hg38.fa"
BED="outs/consensus_peak_calling/cellranger.bed"
FA_OUT="output.fa" 

PYTHON_SCRIPT="create_cisTarget_databases/create_cistarget_motif_databases.py"
MOTIF_COLLECTION="aertslab_motif_colleciton/v10nr_clust_public/singletons"
MOTIFS="motifs.txt"
OUTPUT_DIR="BPD_Control_cellranger" ##need change
THREADS=20
dos2unix $BED
bedtools getfasta -fi $FASTA -bed $BED -fo $FA_OUT

module load anaconda3
source activate scenicplus_3.11
python $PYTHON_SCRIPT -f $FA_OUT -M $MOTIF_COLLECTION -m $MOTIFS -o $OUTPUT_DIR -t $THREADS