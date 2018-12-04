# Set-up Conda
source ~/.bashrc

RDS=$1
NCPU=$2
ID=$3
param=$4
if [ -z $param ]; then
	param=32
fi
MAT=/lustre/scratch117/cellgen/team218/TA/TemporaryFileDir/counts_$ID.csv
OUT=/lustre/scratch117/cellgen/team218/TA/TemporaryFileDir/$ID
mkdir -p $OUT

/software/R-3.4.2/bin/Rscript ~/MAGIC/extract_count_mat.R $RDS $MAT

conda activate dca2

dca --threads=$NCPU -s 64,$param,64 $MAT $OUT

/software/R-3.4.2/bin/Rscript ~/MAGIC/add_count_mat.R $RDS $OUT/mean.tsv "dca"

# clean-up
rm $MAT
rm -r $OUT
