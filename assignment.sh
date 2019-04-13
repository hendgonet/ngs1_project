# assignment
mkdir ~/workdir/assignment
cd ~/workdir/assignment
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra
mv DRR000033.sra sample.sra
sudo apt install sra-toolkit
# convert sra to one file fastq
fastq-dump --split-spot sample.sra  
# get 5 samples each one of million reads
gzip sample.fastq
seqkit split sample.fastq.gz -s 1000000
# shuffle
seqkit shuffle sample.fastq.gz > shuffled.fastq.gz
# get 5 shuffled samples each one of million reads
seqkit split shuffled.fastq.gz -s 1000000
# use FASTQC to report the difference between S1_1 and S1_2
cd shuffled.fastq.gz.split
fastqc  --noexsample.fastq.gz.splittract -f fastq shuffled.part_001.fastq.gz
cd ..
cd sample.fastq.gz.split
fastqc  --noexsample.fastq.gz.splittract -f fastq sample.part_001.fastq.gz
# Mild Trimming for SX_1. {unshuffled}
sudo apt install cutadapt
cutadapt -m 10 -q 20 -o sample1_1_trimmed.fastq.gz sample.part_001.fastq.gz
cutadapt -m 10 -q 20 -o sample1_1_trimmed.fastq.gz sample.part_002.fastq.gz
cutadapt -m 10 -q 20 -o sample1_1_trimmed.fastq.gz sample.part_003.fastq.gz
cutadapt -m 10 -q 20 -o sample1_1_trimmed.fastq.gz sample.part_004.fastq.gz
cutadapt -m 10 -q 20 -o sample1_1_trimmed.fastq.gz sample.part_005.fastq.gz
# Aggressive Trimming for SX_2. {shuffled}
cd ..
cd shuffled.fastq.gz.split
cutadapt -m 50 -q 100 -o sample1_2_trimmed.fastq.gz shuffled.part_001.fastq.gz
cutadapt -m 50 -q 100 -o sample1_2_trimmed.fastq.gz shuffled.part_002.fastq.gz
cutadapt -m 50 -q 100 -o sample1_2_trimmed.fastq.gz shuffled.part_003.fastq.gz
cutadapt -m 50 -q 100 -o sample1_2_trimmed.fastq.gz shuffled.part_004.fastq.gz
cutadapt -m 50 -q 100 -o sample1_2_trimmed.fastq.gz shuffled.part_005.fastq.gz
# BWA Alignment
cd ..
source activate ngs1
mkdir -p ~/workdir/assignment/bwa_align/bwaIndex
cd ~/workdir/assignment/bwa_align/bwaIndex
ln -s /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa .
bwa index -a bwtsw /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa
R1="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_001.fastq.gz"
R2="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_002.fastq.gz"
R3="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_003.fastq.gz"
R4="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_004.fastq.gz"
R5="$HOME/workdir/assignment/sample.fastq.gz.split/sample.part_005.fastq.gz"
/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R1 > chr22_with_ERCC92_bwa.sam
/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R2  > chr22_with_ERCC92_bwa.sam
/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R3  > chr22_with_ERCC92_bwa.sam
/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R4  > chr22_with_ERCC92_bwa.sam
/usr/bin/time -v bwa mem  /home/ngs-01/workdir/assignment/chr22_with_ERCC92.fa  $R5  > chr22_with_ERCC92_bwa.sam
# hisat2 Alignment
source activate ngs1
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
source activate ngs1
conda install samtools
# convert the SAM file into BAM file 
samtools view -bS chr22_with_ERCC92_hisat.sam > chr22_with_ERCC92_hisat.bam
#convert the BAM file to a sorted BAM file. 
samtools sort chr22_with_ERCC92_hisat.bam -o chr22_with_ERCC92_hisat.sorted.bam
source activate ngs1
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

