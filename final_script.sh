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
STAR --runThreadN 1 --runMode genomeGenerate --genomeDir idx/ --genomeFastaFiles ~/workdir2/assignment/ngs2_assignment/sample_data/chr22_with_ERCC92.fa

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
#    --sjdbFileChrStartEnd /path/to/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN <n>            ##
######################################################################                             ##
#default parameters                                                                                ##
#STAR --runMode alignReads --genomeDir ./GenomeDir/ --genomeFastaFiles -                           ##
#   --sjdbFileChrStartEnd - --sjdbOverhang 0 --runThreadN 1                                        ##
#####################################################################################################

#3) For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

#generate indexed genome for 2nd pass alignment:

mkdir ~/workdir2/assignment/ngs2_assignment/genome2 && cd ~/workdir2/assignment/ngs2_assignment/genome2
mkdir idy/
STAR --runMode genomeGenerate --genomeDir idy/ --genomeFastaFiles ~/workdir2/assignment/ngs2_assignment/sample_data/chr22_with_ERCC92.fa --sjdbFileChrStartEnd ~/workdir2/assignment/ngs2_assignment/runDir/SJ.out.tab --sjdbOverhang 75 --runThreadN 1

####################################################################################

#4) The resulting index is then used to produce the final alignments as follows:

#####################################################################################################
#                                                                                                   ##
#runDir=/path/to/2pass                                                                              ##
#mkdir $runDir                                                                                      ##
#cd $runDir                                                                                         ##
#STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN <n>                       ##
#                                                                                                   ##
#####################################################################################################

#4) The resulting index is then used to produce the final alignments as follows:

#R1=~/workdir2/assignment/ngs2_assignment/sample_data/ngs2-assignment-data/SRR8797509_1.part_001.part_001.fastq.
#R2=~/workdir2/assignment/ngs2_assignment/sample_data/ngs2-assignment-data/SRR8797509_2.part_001.part_001.fastq.

cd ~/workdir2/assignment/ngs2_assignment
mkdir runDir2 && cd runDir2
STAR --runThreadN 1 --genomeDir ~/workdir2/assignment/ngs2_assignment/genome2/idy --readFilesIn ~/workdir2/assignment/ngs2_assignment/sample_data/ngs2-assignment-data/SRR8797509_1.part_001.part_001.fastq ~/workdir2/assignment/ngs2_assignment/sample_data/ngs2-assignment-data/SRR8797509_2.part_001.part_001.fastq

#####################################################################################

#2. Add read groups, sort, mark duplicates, and create index

#####################################################################################################
#java -jar picard.jar AddOrReplaceReadGroups I=star_output.sam O=rg_added_sorted.bam SO=coordinate  ## 
#RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample                                        ##
#                                                                                                   ##
#java -jar picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true        ##
#VALIDATION_STRINGENCY=SILENT M=output.metrics                                                      ##
#                                                                                                   ##
#####################################################################################################

# Install Picard tools

conda install -c bioconda picard 
picard_path=$CONDA_PREFIX/share/picard-2.19.0-0
sudo apt install picard-tools

##########################################################

#2. Add read groups, sort, mark duplicates, and create index

#Only the commands for the 2PASS that will be used as it is the final alignment step

#FOR 1PASS

#cd ~/workdir2/assignment/ngs2_assignment/runDir
#java -jar picard.jar AddOrReplaceReadGroups I=Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20                                        
#Error: Unable to access jarfile picard.jar

#cd ~/workdir2/assignment/ngs2_assignment/runDir
#picard-tools AddOrReplaceReadGroups I=Aligned.out.sam O=rg_added_sorted2.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20                                        
###############################       

#java -jar picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics                                                      
#Error: Unable to access jarfile picard.jar

#cd ~/workdir2/assignment/ngs2_assignment/runDir
#picard-tools MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics                                                      

##############################################################################

#fOR 2PASS

#cd ~/workdir2/assignment/ngs2_assignment/runDir2
#java -jar picard.jar AddOrReplaceReadGroups I=Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20                                        
#Error: Unable to access jarfile picard.jar

cd ~/workdir2/assignment/ngs2_assignment/runDir2
picard-tools AddOrReplaceReadGroups I=Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20     

##############################

#cd ~/workdir2/assignment/ngs2_assignment/runDir2
#java -jar picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics                                                      
#Error: Unable to access jarfile picard.jar

cd ~/workdir2/assignment/ngs2_assignment/runDir2
picard-tools MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics  

##############################################################################

#Install GATK
conda install -c bioconda gatk4

##################################################

#3. Split'N'Trim and reassign mapping qualities

######################################################################################################
# java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R ref.fasta -I dedupped.bam -o split.bam      ##
#  -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS                          ##
######################################################################################################

#3. Split'N'Trim and reassign mapping qualities

cd ~/workdir2/assignment/ngs2_assignment
mkdir split && cd split 
GATK=~/miniconda3/envs/ngs1/share/gatk4-4.1.2.0-0/gatk-package-4.1.2.0-local.jar
GenomeFasta=~/workdir2/assignment/ngs2_assignment/sample_data/chr22_with_ERCC92.fa

java -jar $GATK -T SplitNCigarReads -R $GenomeFasta -I ~/workdir2/assignment/ngs2_assignment/runDir2/dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS                           


###################################################################################

#4. Indel Realignment by GATK (optional)

cd ~/workdir2/assignment/ngs2_assignment
mkdir realignment && cd realignment 
GATK=~/miniconda3/envs/ngs1/share/gatk4-4.1.2.0-0/gatk-package-4.1.2.0-local.jar
GenomeFasta=~/workdir2/assignment/ngs2_assignment/sample_data/chr22_with_ERCC92.fa

java -jar $GATK -T IndelRealigner -R $GenomeFasta -I ~/workdir2/assignment/ngs2_assignment/split/split.bam -known indels.vcf -targetIntervals intervalListFromRTC.intervals -o realigned.bam

####################################################################################

#5. Base Recalibration by GATK  

cd ~/workdir2/assignment/ngs2_assignment
mkdir BaseRecalibrate && cd BaseRecalibrate
GATK=~/miniconda3/envs/ngs1/share/gatk4-4.1.2.0-0/gatk-package-4.1.2.0-local.jar
GenomeFasta=~/workdir2/assignment/ngs2_assignment/sample_data/chr22_with_ERCC92.fa

java -jar $GATK -T BaseRecalibrator -R $GenomeFasta -I ~/workdir2/assignment/ngs2_assignment/realignment/realigned.bam -knownSites latest_dbsnp.vcf -o recal_data.table

###################################################################################

#6. Variant calling

##########################################################################################################
#                                                                                                       ##
#java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ref.fasta -I input.bam -dontUseSoftClippedBases \ ##
#-stand_call_conf 20.0 -o output.vcf                                                                    ##
#                                                                                                       ##
##########################################################################################################

#6. Variant calling

cd ~/workdir2/assignment/ngs2_assignment
mkdir Vcall && Vcall 
GATK=~/miniconda3/envs/ngs1/share/gatk4-4.1.2.0-0/gatk-package-4.1.2.0-local.jar
GenomeFasta=~/workdir2/assignment/ngs2_assignment/sample_data/chr22_with_ERCC92.fa

java -jar $GATK -T HaplotypeCaller -R $GenomeFasta -I input.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.vcf 

########################################################################################

#7. Variant filtering

##########################################################################################################
#                                                                                                       ##
# java -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg_19.fasta -V input.vcf -window 35  \         ##
#  -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o output.vcf        ##
#                                                                                                       ##
##########################################################################################################

#7. Variant filtering

cd ~/workdir2/assignment/ngs2_assignment
mkdir Vcall && Vcall 
GATK=~/miniconda3/envs/ngs1/share/gatk4-4.1.2.0-0/gatk-package-4.1.2.0-local.jar
GenomeFasta=~/workdir2/assignment/ngs2_assignment/sample_data/chr22_with_ERCC92.fa

java -jar $GATK -T VariantFiltration -R $GenomeFasta -V input.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o output.vcf                  


