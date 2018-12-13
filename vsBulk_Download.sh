# Tung et al.
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77288/suppl/GSE77288_molecules-raw-single-per-sample.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77288/suppl/GSE77288_reads-raw-bulk-per-sample.txt.gz


#Kolo et al.
wget https://www.ebi.ac.uk/teichmann-srv/espresso/static/counttable_es.csv
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR523/ERR523096/ERR523096_1.fastq.gz #a2i
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR523/ERR523096/ERR523096_2.fastq.gz #a2i
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR522/ERR522876/ERR522876_1.fastq.gz # a2i
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR522/ERR522876/ERR522876_2.fastq.gz # a2i
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR523/ERR523093/ERR523093_1.fastq.gz # serum
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR523/ERR523093/ERR523093_2.fastq.gz # serum
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR523/ERR523027/ERR523027_1.fastq.gz # serum
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR523/ERR523027/ERR523027_2.fastq.gz # serum

# Quantify Kolo bulk
STAR=~/RNASeqPipeline/software/STAR-STAR_2.4.0j/bin/Linux_x86_64_static/STAR
GENOME=/lustre/scratch117/cellgen/team218/TA/STRIPED_GENOMES/STAR_Mouse
NTHREADS=1
PARAMS=~/RNASeqPipeline/2_STAR_Parameters.txt
OUTDIR=/lustre/scratch117/cellgen/team218/TA/TemporaryFileDir
FC=~/RNASeqPipeline/software/subread-1.4.6-p2-Linux-x86_64/bin/featureCounts
GTF=/lustre/scratch117/cellgen/team218/TA/genomebuilding/Mus_musculus.GRCm38.79.gtf

$STAR --runThreadN $NTHREADS --runMode alignReads --genomeDir $GENOME --readFilesCommand zcat --parametersFiles $PARAMS --outFileNamePrefix $OUTDIR/a2i_1 --outTmpDir $OUTDIR/a2i_1 --readFilesIn ERR523096_1.fastq.gz ERR523096_2.fastq.gz 

$STAR --runThreadN $NTHREADS --runMode alignReads --genomeDir $GENOME --readFilesCommand zcat --parametersFiles $PARAMS --outFileNamePrefix $OUTDIR/a2i_2 --outTmpDir $OUTDIR/a2i_2 --readFilesIn ERR522876_1.fastq.gz ERR522876_2.fastq.gz 

$STAR --runThreadN $NTHREADS --runMode alignReads --genomeDir $GENOME --readFilesCommand zcat --parametersFiles $PARAMS --outFileNamePrefix $OUTDIR/serum_1 --outTmpDir $OUTDIR/serum_1 --readFilesIn ERR523093_1.fastq.gz ERR523093_2.fastq.gz 

$STAR --runThreadN $NTHREADS --runMode alignReads --genomeDir $GENOME --readFilesCommand zcat --parametersFiles $PARAMS --outFileNamePrefix $OUTDIR/serum_2 --outTmpDir $OUTDIR/serum_2 --readFilesIn ERR523027_1.fastq.gz ERR523027_2.fastq.gz 


$FC -O -p -T $NTHREADS -a $GTF -o $OUTDIR/a2i_1.fragmentcounts $OUTDIR/a2i_1Aligned.sortedByCoord.out.bam
$FC -O -p -T $NTHREADS -a $GTF -o $OUTDIR/a2i_2.fragmentcounts $OUTDIR/a2i_2Aligned.sortedByCoord.out.bam
$FC -O -p -T $NTHREADS -a $GTF -o $OUTDIR/serum_1.fragmentcounts $OUTDIR/serum_1Aligned.sortedByCoord.out.bam
$FC -O -p -T $NTHREADS -a $GTF -o $OUTDIR/serum_2.fragmentcounts $OUTDIR/serum_2Aligned.sortedByCoord.out.bam

perl ~/MAGIC/6_Get_Expression_featureCounts.pl $OUTDIR Kolo_Bulk.txt

#bsub -R"select[mem>3000] rusage[mem=3000]" -M3000 -o count_bulk.%J.out -e count_bulk.%J.err
#bsub -R"select[mem>30000] rusage[mem=30000]" -M30000 -o map_bulk.%J.out -e map_bulk.%J.err

