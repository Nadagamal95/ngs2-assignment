#############################################
##   - NGS2_course                         ##
##   - NGS2_Assignment                     ##
##   - Bash script                         ##
##   - Assignment_Steps:-                  ##
##   - Date: 25 April 2019                 ##
##   - copyright Nada Gamal                ##
##   - Nile University                     ##
#############################################

#############################################

source activate ngs1

#Download Fastq files

#Download the fasta files of ngs2-assignment-data.zip from: 

https://uploadfiles.io/kc0qqbvd

##############################################

#make a fqData file & move the downloaded zipped fasta files to it

mkdir -p ~/ngs2_assignment/fqData
cd ~/Downloads
mv ngs2-assignment-data.zip ~/ngs2_assignment/fqData && cd ~/ngs2_assignment/fqData/ 

###############################################

#unzip the fasta files "ngs2-assignment-data.zip"

unzip ngs2-assignment-data.zip

###############################################

#extract the files compressed with gunzip

cd ~/ngs2_assignment/fqData/ngs2-assignment-data/

for ((i=1;i<=2;i++)) ; do
gunzip shuffled_SRR8797509_$i.part_001.part_001.fastq.gz

gunzip SRR8797509_$i.part_001.part_001.fastq.gz

done

###############################################

#Download reference file

mkdir -p ~/ngs2_assignment/sample_data && cd ~/ngs2_assignment/sample_data

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

#uncompress the fa.gz files:

gunzip hg38.fa.gz

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz
gunzip gencode.v29.pc_transcripts.fa.gz

# Download the Transcriptome Annotation File

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz

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

#1) generate indexed genome for 1st pass alignment

mkdir ~/ngs2_assignment/STAR_index2 && cd ~/ngs2_assignment/STAR_index2

GenomeDir=~/ngs2_assignment/sample_data
GenomeFasta=~/ngs2_assignment/sample_data/gencode.v29.pc_transcripts.fa

STAR --genomeDir $GenomeDir --runMode genomeGenerate --genomeFastaFiles $GenomeFasta \
 --sjdbGTFfile ~/ngs2_assignment/sample_data/gencode.v29.annotation.gtfAnnotation --sjdbOverhang 100 --runThreadN 8

#STAR --genomeDir $GenomeDir --runMode genomeGenerate --genomeFastaFiles ~/ngs2_assignment/STAR_index2/gencode.v29.pc_transcripts.fa --sjdbGTFfile ~/ngs2_assignment/STAR_index2/gencode.v29.annotation.gtfAnnotation --sjdbOverhang 75 --runThreadN 1
## '--sjdbOverhang' is 'ReadLength-1'

########################################################

#2) Alignment jobs were executed as follows:

#runDir=~/ngs2_assignment/STAR_index2/1pass
#mkdir $runDir
#cd $runDir
#STAR --genomeDir $genomeDir --readFilesIn ~/ngs2_assignment/STAR_index2/SRR8797509_1.part_001.part_001.fastq ~/ngs2_assignment/STAR_index2/SRR8797509_2.part_001.part_001.fastq --runThreadN 1

#ERROR
######################################################

# run 1st pass

mkdir ~/ngs2_assignment/STAR_index2/1Pass && cd ~/ngs2_assignment/STAR_index2/1Pass

#Read1=~/ngs2_assignment/fqData/ngs2_assignment_data/SRR8797509_1.part_001.part_001.fastq
#Read2=~/ngs2_assignment/fqData/ngs2_assignment_data/SRR8797509_2.part_001.part_001.fastq

STAR --genomeDir ~/ngs2_assignment/sample_data/genome --readFilesIn ~/ngs2_assignment/fqData/ngs2_assignment_data/SRR8797509_1.part_001.part_001.fastq ~/ngs2_assignment/fqData/ngs2_assignment_data/SRR8797509_2.part_001.part_001.fastq --runThreadN 8

#ERROR 
# I TRIED EVERYTHING AND I COULDNOT SOLVE MY ERROR 
################################################################
	
# generate genome with junctions from the 1st pass

genomeDir=~/ngs2_assignment/sample_data/gencode.v29.pc_transcripts.fa
mkdir $genomeDir
STAR --genomeDir $genomeDir --runMode genomeGenerate --genomeFastaFiles $GenomeFasta \
--sjdbFileChrStartEnd - --runThreadN 8

#########################################################

#4) The resulting index is then used to produce the final alignments as follows:

mkdir ~/ngs2_assignment/STAR_index2/2Pass && cd ~/ngs2_assignment/STAR_index2/2Pass
STAR --genomeDir ~/ngs2_assignment/sample_data/genome  --readFilesIn ~/ngs2_assignment/fqData/ngs2_assignment_data/SRR8797509_1.part_001.part_001.fastq ~/ngs2_assignment/fqData/ngs2_assignment_data/SRR8797509_2.part_001.part_001.fastq --runThreadN 8

#ERROR

