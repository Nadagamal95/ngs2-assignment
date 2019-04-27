
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

mkdir ~/ngs2_assignment/STAR_index && cd ~/ngs2_assignment/STAR_index
GenomeDir=~/ngs2_assignment/sample_data/gencode.v29.pc_transcripts.fa
GenomeFasta=~/ngs2_assignment/sample_data/gencode.v29.pc_transcripts.fa
STAR --genomeDir $GenomeDir --runMode genomeGenerate --runThreadN 1 --genomeFastaFiles $GenomeFasta\

# I tried to run it & I got error

#ERROR:
#(ngs1) ngs-01@ngs01-VirtualBox:~/ngs2_assignment/STAR_index$ STAR --genomeDir $GenomeDir --runMode genomeGenerate --genomeFastaFiles $GenomeFasta \ --runThreadN 1
#Apr 27 05:35:33 ..... started STAR run
#Apr 27 05:35:33 ... starting to generate Genome files
#EXITING because of INPUT ERROR: could not open genomeFastaFile:  --runThreadN
#Apr 27 05:36:12 ...... FATAL ERROR, exiting

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

# I tried to run it & I got the same error

################################################################################

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

# generate indexed genome for 1st pass alignment

mkdir ~/ngs2_assignment/STAR_index && cd ~/ngs2_assignment/STAR_index

GenomeDir=~/ngs2_assignment/sample_data/
GenomeFasta=~/ngs2_assignment/sample_data/gencode.v29.pc_transcripts.fa

STAR --genomeDir $GenomeDir --runMode genomeGenerate --genomeFastaFiles $GenomeFasta \
 --sjdbGTFfile ~/ngs2_assignment/sample_data/gencode.v29.annotation.gtfAnnotation --sjdbOverhang 100 --runThreadN 8
## '--sjdbOverhang' is 'ReadLength-1'

# result

#Apr 27 07:18:02 ..... started STAR run
#Apr 27 07:18:02 ... starting to generate Genome files
#terminate called after throwing an instance of 'std::bad_alloc'
#  what():  std::bad_alloc
#Aborted (core dumped)

##################################

#STAR --genomeDir $GenomeDir --runMode genomeGenerate --genomeFastaFiles $GenomeFasta \
# --sjdbGTFfile ~/ngs2_assignment/sample_data/gencode.v29.annotation.gtfAnnotation --sjdbOverhang 75 --runThreadN 4
## '--sjdbOverhang' is 'ReadLength-1'

#result 

#Apr 27 07:18:02 ..... started STAR run
#Apr 27 07:18:02 ... starting to generate Genome files
#terminate called after throwing an instance of 'std::bad_alloc'
#  what():  std::bad_alloc
#Aborted (core dumped)

################################################################################

#4) The resulting index is then used to produce the final alignments as follows:

#####################################################################################################
                                                                                                   ##
#runDir=/path/to/2pass                                                                             ##
mkdir $runDir                                                                                      ##
cd $runDir                                                                                         ##
STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN <n>                       ##
                                                                                                   ##
#####################################################################################################

# run 1st pass

mkdir ~/ngs2_assignment/Pass1 && cd ~/ngs2_assignment/Pass1

GenomeDir=~/ngs2_assignment/sample_data/
GenomeFasta=~/ngs2_assignment/sample_data/gencode.v29.pc_transcripts.fa
Reads=~/ngs2_assignment/fqData/ngs2_assignment_data
CommonPars= --runThreadN 8 --outSAMattributes All

for f in $Reads; do \
STAR $CommonPars --genomeDir $GenomeDir --readFilesIn $f --outFileNamePrefix ./ ; done

###################################################################

# run 1st pass

mkdir ~/ngs2_assignment/STAR_index2/1Pass && cd ~/ngs2_assignment/STAR_index2/1Pass

GenomeDir=~/ngs2_assignment/sample_data/gencode.v29.pc_transcripts.fa
GenomeFasta=~/ngs2_assignment/sample_data/gencode.v29.pc_transcripts.fa
Read1="~/ngs2_assignment/fqData/ngs2_assignment_data/SRR8797509_1.part_001.part_001.fastq"
Read2="~/ngs2_assignment/fqData/ngs2_assignment_data/SRR8797509_2.part_001.part_001.fastq"

STAR --genomeDir $GenomeDir --readFilesIn $Read1 $Read2 --runThreadN 1

#######################################################################

# make splice junctions database file out of SJ.out.tab, 
# filter out non-canonical junctions 
# (if have multiple files, combine them into one and conduct filtering)

mkdir ~/ngs2_assignment/GenomeForPass2
cd ~/ngs2_assignment/GenomeForPass2

awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} \
{if($5 > 0){print $1,$2,$3,strChar[$4]}}' \
../Pass1/SJ.out.tab > SJ.out.tab.Pass1.sjdb	

######################################################################

# generate genome with junctions from the 1st pass
$STAR --genomeDir ./ --runMode genomeGenerate --genomeFastaFiles $GenomeFasta \
--sjdbFileChrStartEnd SJ.out.tab.Pass1.sjdb --sjdbOverhang 100 --runThreadN 1
cd ..

####################################################################

#3) For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

genomeDir=~/ngs2_assignment/STAR_index2/gencode.v29.pc_transcripts.fa
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa \
    --sjdbFileChrStartEnd /path/to/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN <n>

##################################################################

# generate genome with junctions from the 1st pass

genomeDir=~/ngs2_assignment/sample_data/gencode.v29.pc_transcripts.fa
mkdir $genomeDir
STAR --genomeDir $genomeDir --runMode genomeGenerate --genomeFastaFiles $GenomeFasta \
--sjdbFileChrStartEnd - --runThreadN 8

####################################################################

2. Add read groups, sort, mark duplicates, and create index
The above step produces a SAM file, which we then put through the usual Picard processing steps: adding read group information, sorting, marking duplicates and indexing.

java -jar picard.jar AddOrReplaceReadGroups I=star_output.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample 

java -jar picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 

################################################################

#2. Add read groups, sort, mark duplicates, and create index
The above step produces a SAM file, which we then put through the usual Picard processing steps: adding read group information, sorting, marking duplicates and indexing.

java -jar picard.jar AddOrReplaceReadGroups I=star_output.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample 

java -jar picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 

################################################################

#Install GATK
conda install -c bioconda gatk4
 
######################################
 
#3. Split'N'Trim and reassign mapping qualities

java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R ref.fasta -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

################################################################

#4. Indel Realignment (optional)

# GATK Indel Realignment (Optional)
java -jar $GATK -T IndelRealigner -R $GenomeFasta -I split.bam \
-known indels.vcf -targetIntervals intervalListFromRTC.intervals \
-o realigned.bam

##############################################################

#5. Base Recalibration

# GATK Base recalibration (highly recommended, but not works without known SNP data. 
# Skip this step, if can't find dbSNP.vcf file for the organism)
j#ava -jar $GATK -T BaseRecalibrator -R $GenomeFasta -I realigned.bam \
#-knownSites latest_dbsnp.vcf -o recal_data.table; done

for sample in *.dedup.bam;do
  name=${sample%.dedup.bam}

  gatk --java-options "-Xmx2G" BaseRecalibrator \
-R dog_chr5.fa -I $sample --known-sites canis_fam_chr5.vcf \
-O $name.report

  gatk --java-options "-Xmx2G" ApplyBQSR \
-R dog_chr5.fa -I $sample -bqsr $name.report \
-O $name.bqsr.bam --add-output-sam-program-record --emit-original-quals
done


###################################################3

#6. Variant calling

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ref.fasta -I input.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.vcf

####################################################

#7. Variant filtering

java -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg_19.fasta -V input.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o output.vcf 
