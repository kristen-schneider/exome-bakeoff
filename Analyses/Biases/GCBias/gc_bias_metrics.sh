echo "-------"
echo $1
echo $2
echo $3
echo $4
echo "-------"
java -jar $PICARD CollectGcBiasMetrics \
I=$1 \
O=$2 \
CHART=$3 \
S=$4 \
R=/scratch/Shares/layer/ref/hg37/human_g1k_v37.fasta.gz
