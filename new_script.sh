#############################################
##   - NGS2_course                         ##
##   - NGS2_Assignment                     ##
##   - Bash script                         ##
##   - Assignment_Steps:-                  ##
##   - Date: 26 April 2019                 ##
##   - copyright Nada Gamal                ##
##   - Nile University                     ##
#############################################

#############################################

source activate ngs1

#Download Fastq files

#Download the fasta files of ngs2-assignment-data.zip from: 

https://uploadfiles.io/kc0qqbvd

##############################################

#move the downloaded zipped fasta files to it

mkdir -p ~/workdir2/assignment/ngs2_assignment/sample_data 
cd ~/Downloads
mv ngs2-assignment-data.zip ~/workdir2/assignment/ngs2_assignment/sample_data && cd ~/workdir2/assignment/ngs2_assignment/sample_data

###############################################

#unzip the fasta files "ngs2-assignment-data.zip"

unzip ngs2-assignment-data.zip

###############################################

#extract the files compressed with gunzip

cd ~/workdir2/assignment/ngs2_assignment/sample_data/ngs2-assignment-data/

for ((i=1;i<=2;i++)) ; do
gunzip shuffled_SRR8797509_$i.part_001.part_001.fastq.gz

gunzip SRR8797509_$i.part_001.part_001.fastq.gz

done

###############################################

#Download reference file (chr22_with_ERCC92)

cd ~/workdir2/assignment/ngs2_assignment/sample_data
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa

gunzip chr22_with_ERCC92.fa.gz

wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf

################################################

#STAR installation

sudo apt-get update
sudo apt-get install g++
sudo apt-get install make
sudo apt install rna-star

#################################################

#The workflow for RNA variant calling

https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq

#1. Mapping to the reference by using STAR aligner(Specifically, we use the STAR 2-pass method)

#STAR 2-pass alignment steps:

#1) STAR uses genome index files that must be saved in unique directories. The human genome index was built from the FASTA file hg19.fa

#####################################################################################################
#genomeDir=/path/to/hg19                                                                           ##
#mkdir $genomeDir                                                                                  ##
#STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa\  --runThreadN <n>##
######################################################################                             ##
#default parameters                                                                                ##
#STAR --runMode alignReads --genomeDir ./GenomeDir/ --genomeFastaFiles -  --runThreadN 1           ##
#####################################################################################################

#1)STAR uses genome index files that must be saved in unique directories.The human genome index was built from the FASTA file gencode.v29.pc_transcripts.fa

#generate indexed genome for 1st pass alignment:

mkdir ~/workdir2/assignment/ngs2_assignment/genome && cd ~/workdir2/assignment/ngs2_assignment/genome
mkdir idx/
STAR --runThreadN 1 --runMode genomeGenerate --genomeDir idx/ --genomeFastaFiles ~/workdir_2/assignment/ngs2_assignment/sample_data/chr22_with_ERCC92.fa

###########################################################################

#2) Alignment jobs were executed as follows:

#####################################################################################################
#runDir=/path/to/1pass                                                                             ##
#mkdir $runDir                                                                                     ##
#cd $runDir                                                                                        ##
#STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN <n>                      ##
######################################################################                             ##
#default parameters                                                                                ##
#STAR --genomeDir ./GenomeDir/ --readFilesIn Read1 Read2  --runThreadN 1                           ##
#####################################################################################################

#2) Alignment jobs were executed as follows:

#R1=~/workdir2/assignment/ngs2_assignment/sample_data/ngs2-assignment-data/SRR8797509_1.part_001.part_001.fastq.
#R2=~/workdir2/assignment/ngs2_assignment/sample_data/ngs2-assignment-data/SRR8797509_2.part_001.part_001.fastq.

cd ~/workdir2/assignment/ngs2_assignment
mkdir runDir && cd runDir
STAR --runThreadN 1 --genomeDir ~/workdir2/assignment/ngs2_assignment/genome/idx --readFilesIn ~/workdir2/assignment/ngs2_assignment/sample_data/ngs2-assignment-data/SRR8797509_1.part_001.part_001.fastq ~/workdir2/assignment/ngs2_assignment/sample_data/ngs2-assignment-data/SRR8797509_2.part_001.part_001.fastq

#################################################################################

#3) For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

#####################################################################################################
#genomeDir=/path/to/hg19_2pass                                                                     ##
#mkdir $genomeDir                                                                                  ##
#STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa \                 ##
    --sjdbFileChrStartEnd /path/to/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN <n>             ##
######################################################################                             ##
#default parameters                                                                                ##
#STAR --runMode alignReads --genomeDir ./GenomeDir/ --genomeFastaFiles -                           ##
#   --sjdbFileChrStartEnd - --sjdbOverhang 0 --runThreadN 1                                        ##
#####################################################################################################

#3) For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

#generate indexed genome for 2nd pass alignment:

mkdir ~/workdir2/assignment/ngs2_assignment/genome && cd ~/workdir2/assignment/ngs2_assignment/genome
mkdir idy/
STAR --runMode genomeGenerate --genomeDir idy/ --genomeFastaFiles ~/workdir_2/assignment/ngs2_assignment/sample_data/chr22_with_ERCC92.fa --sjdbFileChrStartEnd /path/to/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN 1
