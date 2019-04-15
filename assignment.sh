# assignment
mkdir ~/workdir/assignment
cd ~/workdir/assignment
source activate ngs1
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR746/SRR7462222/SRR7462222.sra  
mv SRR7462222.sra sample.sra 
sudo apt install sra-toolkit
# convert sra to one file fastq
fastq-dump --split-spot sample.sra  
# get 5 samples each of 460000 reads
gzip sample.fastq
seqkit split sample.fastq.gz -s 460000
# shuffle
seqkit shuffle sample.fastq.gz > shuffled.fastq.gz
# get 5 shuffled samples each one of 460000 reads
seqkit split shuffled.fastq.gz -s 460000
# use FASTQC to report the difference between S1_1 and S1_2
cd sample.fastq.gz.split
fastqc  --noexsample.fastq.gz.splittract -f fastq sample.part_001.fastq.gz
fastqc  --noexsample.fastq.gz.splittract -f fastq sample.part_002.fastq.gz
# Mild Trimming for SX_1. {unshuffled}
sudo apt install cutadapt
for i in {1..5};do cutadapt -m 10 -q 20 -o 1trimmed$i.fastq.gz sample.part_00$i.fastq.gz; done
# Aggressive Trimming for SX_2. {shuffled}
cd ..
cd shuffled.fastq.gz.split
for i in {1..5};do cutadapt -m 50 -q 100 -o 2trimmed$i.fastq.gz shuffled.part_00$i.fastq.gz; done
# BWA Alignment
cd ..
mkdir -p ~/workdir/assignment/bwa_align/bwaIndex
cd ~/workdir/assignment/bwa_align/bwaIndex
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
ln -s /home/ngs-01/workdir/assignment/bwa_align/bwaIndex/chr22_with_ERCC92.fa .
bwa index -a bwtsw /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa
R1="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_001.fastq.gz"
R2="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_002.fastq.gz"
R3="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_003.fastq.gz"
R4="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_004.fastq.gz"
R5="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_005.fastq.gz"
/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R1  > chr22_with_ERCC92_bwa.sam

/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R2  > chr22_with_ERCC92_bwa.sam
/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R3  > chr22_with_ERCC92_bwa.sam
/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R4  > chr22_with_ERCC92_bwa.sam
/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R5  > chr22_with_ERCC92_bwa.sam
# hisat2 Alignment
conda install -c bioconda hisat2 
mkdir -p ~/workdir/assignment/hisat_align/hisatIndex
cd ~/workdir/assignment/hisat_align/hisatIndex
ln -s ~/home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa .
hisat2_extract_splice_sites.py /home/ngs-01/workdir/assignment/chr22_with_ERCC92.gtf > splicesites.tsv
hisat2_extract_exons.py /home/ngs-01/workdir/assignment/chr22_with_ERCC92.gtf > exons.tsv
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv chr22_with_ERCC92.fa chr22_with_ERCC92
cd ~/workdir/hisat_align
R1="$HOME/workdir/assignment/shuffled.fastq.gz.split/shuffled.part_001.fastq.gz"
R2="$HOME/workdir/assignment/shuffled.fastq.gz.split/shuffled.part_001.fastq.gz"
R3="$HOME/workdir/assignment/shuffled.fastq.gz.split/shuffled.part_001.fastq.gz"
R4="$HOME/workdir/assignment/shuffled.fastq.gz.split/shuffled.part_001.fastq.gz"
R5="$HOME/workdir/assignment/shuffled.fastq.gz.split/shuffled.part_001.fastq.gz"
hisat2 -p 1 -x hisatIndex/chr22_with_ERCC92 --dta --rna-strandness RF -1 $R1 -2 $R2 -S chr22_with_ERCC92_hisat.sam
# Assembly
# Apply reference-based trasncriptome assembly using stringTie.
# Step1. For the 5 samples unshuffled.
# Step2. For the 5 samples shuffled.
# install Samtools
conda install samtools
# convert the SAM file into BAM file 
samtools view -bS chr22_with_ERCC92_hisat.sam > chr22_with_ERCC92_hisat.bam
#convert the BAM file to a sorted BAM file. 
samtools sort chr22_with_ERCC92_hisat.bam -o chr22_with_ERCC92_hisat.sorted.bam
conda install stringtie
stringtie chr22_with_ERCC92_hisat.sorted.bam --rf -l ref_free -o ref_free.gtf
## how many transcript do you have?
cat ref_free.gtf | grep -v "^@" | awk '$3=="transcript"' | wc -l
# Using GTF-Compare to Compare the Generated Annotation Files to a Reference Annotation.
# create virtual evironment with conda
conda create -n ngs-gtf python=3.6 anaconda
source activate ngs-gtf
conda install -c conda-forge pypy3.5
wget https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py
source activate ngs-gtf
conda install gffcompare
mkdir -p ~/workdir/assignment/gtf-compare/method_two && cd ~/workdir/assignment/gtf-compare/method_two
gffcompare -r ../gtfs/ref_sup.gtf ../gtfs/ref_free.gtf
# Apply Differential Expression.
set -euo pipefail ## stop execution on errors, https://explainshell.com/explain?cmd=set+-euxo%20pipefail

# Collect program output here.
RUNLOG=runlog.log

READS=~/workdir/assignment/sample.fastq.gz.split
REF_ERCC=~/workdir/assignment/chr22_with_ERCC92.fa # Reference
INDEX_ERCC=~/workdir/assignment/chr22_with_ERCC92.id

