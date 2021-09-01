#!/usr/bin/env bash

fastqdir=''
maindir=''

mkdir $maindir

#add sample names here same as used in sample sheet
inputfiles='Rep1_siControl_cyto, Rep1_siControl_ER, Rep1_siCNOT1_cyto, Rep1_siCNOT1_ER, Rep2_siControl_cyto, Rep2_siControl_ER, Rep2_siCNOT1_cyto, Rep2_siCNOT1_ER, Rep3_siControl_cyto, Rep3_siControl_ER, Rep3_siCNOT1_cyto, Rep3_siCNOT1_ER'


#initial QC check of data
fastqcdir=$maindir/fastqc_initial
mkdir $fastqcdir

for s in ${inputfiles//,/ }
do	
  fastqc $fastqdir/${s}*.fastq.gz --extract --outdir=$fastqcdir &
done
wait


#adapter and quality trimming
fastpdir=$maindir/fastp
mkdir $fastpdir

for s in ${inputfiles//,/ }
do	
  (fastp -i $fastqdir/${s}*.fastq.gz -w 4 -q 20 -l 30 -o $fastpdir/${s}_fastp.fq.gz) 2> fastp_stats_${s}.txt &
done
wait


#QC check after trimming
fastqcdir2=$maindir/fastqc_fastp
mkdir $fastqcdir2

for s in ${inputfiles//,/ }
do	
  fastqc $fastpdir/${s}_fastp.fq.gz --extract --outdir=$fastqcdir2 &
done
wait

#deduplication based on UMI in CORALL kit
for s in ${inputfiles//,/ }
do	
  gunzip $fastpdir/${s}_fastp.fq.gz
done
wait

dedupdir=$maindir/cdhitdup
mkdir $dedupdir

for s in ${inputfiles//,/ }
do
  /home/sgillen/Programs/cd-hit-v4.6.6-2016-0711/cd-hit-auxtools/cd-hit-dup -i $fastpdir/${s}_fastp.fq -o $dedupdir/${s}_cdhitdup.fq -e 0
done
wait


#check deduplication worked as expected
fastqcdir3=$maindir/fastqc_cdhitdup
mkdir $fastqcdir3

for s in ${inputfiles//,/ }
do	
  fastqc $dedupdir/${s}_cdhitdup.fq --extract --outdir=$fastqcdir3 &
done
wait


#remove the UMI before alignment 
UMIremdir=$maindir/UMIremoved
mkdir $UMIremdir

for s in ${inputfiles//,/ } 
do
  cutadapt --nextseq-trim=20 -u 13 -u -2 -o $UMIremdir/${s}_UMIremoved.fq $dedupdir/${s}_cdhitdup.fq & 
done
wait

#rezip files to save space
for s in ${inputfiles//,/ }
do	
  gzip $dedupdir/${s}_cdhitdup.fq 
done
wait



indexdir=''
#code to generate STAR indexes on human genome fasta file combined with ERCC spike-in sequence fasta file
#STAR --runMode genomeGenerate --genomeDir ./STARindex_withERCC --genomeFastaFiles ./genome_withERCCspikeins.fa --runThreadN 24
aligndir=$maindir/STARaligned
mkdir $aligndir

for s in ${inputfiles//,/ }
do
  STAR --readFilesIn $UMIremdir/${s}_UMIremoved.fq --runThreadN 30 --genomeDir $indexdir/genome_withERCCspikein/STARindex_withERCC --outFilterMultimapNmax 5 --outFilterMismatchNmax 2 --outFilterMatchNmin 18 --outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --outSAMtype BAM Unsorted --outSAMunmapped Within KeepPairs --sjdbGTFfile $indexdir/genome_withERCCspikein/gencode_annotation_withERCCspikeins.gtf --outFileNamePrefix $aligndir/${s}_STAR
done
wait


#sort & index BAM files 
for s in ${inputfiles//,/ };do
  samtools sort $aligndir/${s}_STARAligned.out.bam -o $aligndir/${s}_STARAligned_sorted.bam & #sort bam
done
wait

for s in ${inputfiles//,/ };do
  samtools index $aligndir/${s}_STARAligned_sorted.bam $aligndir/${s}_STARAligned_sorted.bai & #index bam
done
wait


#get read counts
countsdir=$maindir/fcounts
mkdir $countsdir

for s in ${inputfiles//,/ };do
    featureCounts -a $indexdir/genome_withERCCspikein/gencode_annotation_withERCCspikeins.gtf -s 1 -o $countsdir/fcounts_${s}.txt $aligndir/${s}_STARAligned_sorted.bam &
done
wait

for s in ${inputfiles//,/ };do
    cut -f1,7-8 $countsdir/fcounts_${s}.txt > $countsdir/fcounts_${s}_matrix.txt
done
wait





