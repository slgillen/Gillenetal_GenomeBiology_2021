#!/usr/bin/bash
inputdir=''
bcldatadir=''
datadir=''

#countsfiles from featureCounts as inputs for siControl ('con') and siCNOT1 ('not') treated samples
inputfiles='con_0hr, con_05hr, con_1hr, con_2hr, con_4hr, con_8hr, con_16hr, not_0hr, not_05hr, not_1hr, not_2hr, not_4hr, not_8hr, not_16hr' 

#separate samples from sequencing lane
bcl2fastq -R $bcldatadir -o $datadir/libraries -r 5 -w 5 -p 20 --sample-sheet $bcldatadir/SampleSheet.csv --barcode-mismatches 0 --no-lane-splitting

#check data quality
fastqcdir=$datadir/initial_fastqc
mkdir $fastqcdir
for s in ${inputfiles//,/ }
do	
  fastqc $datadir/libraries/${s}*.fastq.gz --extract --outdir=$fastqcdir &
done
wait

#use fastp for data trimming
for s in ${inputfiles//,/ }
do	
  fastp -i $datadir/libraries/${s}*.fastq.gz -w 4 -q 20 -l 20 -o $datadir/${s}_fastp.fq.gz &
done
wait

#unzip the files
for s in ${inputfiles//,/ }
do	
  gunzip $datadir/${s}_fastp.fq.gz 
done
wait

#deduplicate reads based on 12nt UMI
for s in ${inputfiles//,/ }
do
  /home/sgillen/Programs/cd-hit-v4.6.6-2016-0711/cd-hit-auxtools/cd-hit-dup -i $datadir/${s}_fastp.fq -o $datadir/${s}_cdhitdup.fq -e 0
done
wait

#check data after deduplication
fastqcdir=$datadir/fastqc_after_deduplication
mkdir $fastqcdir

for s in ${inputfiles//,/ }
do	
  fastqc $datadir/${s}_cdhitdup.fq --extract --outdir=$fastqcdir &
done
wait


#remove UMIs
for s in ${inputfiles//,/ } 
do
  cutadapt -u 13 -o $datadir/${s}_UMIremoved.fq $datadir/${s}_cdhitdup.fq & 
done
wait


#alignments
aligndir=$datadir/STARalignments
mkdir $aligndir

annotation=$inputdir/gencode.v28.annotation.mostabundantTranscript.gtf
genome=$inputdir/GRCh38.primary_assembly.genome.fa

for s in ${inputfiles//,/ }
do
    STAR --readFilesIn $datadir/${s}_UMIremoved.fq --runThreadN 28 --genomeDir $inputdir --outFilterMultimapNmax 5 --outFilterMismatchNmax 2 --outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --outSAMtype BAM Unsorted --outSAMunmapped Within KeepPairs --sjdbGTFfile ${annotation} --outFileNamePrefix $aligndir/${s} 
done
wait

#sort bam
for s in ${inputfiles//,/ };do
  samtools sort $aligndir/${s}Aligned.out.bam -o $aligndir/${s}_sorted.bam & 
done
wait

#index bam
for s in ${inputfiles//,/ };do
  samtools index $aligndir/${s}_sorted.bam $aligndir/${s}_sorted.bai & 
done
wait


#count reads
countsdir=$datadir/counts
mkdir $countsdir

for s in ${inputfiles//,/ }
do
    featureCounts -a ${annotation} -s 1 -o $countsdir/counts_${s}.txt $aligndir/${s}_sorted.bam &
done
wait

#format counts table
for s in ${inputfiles//,/ }
do
    cut -f1,7-8 $countsdir/counts_${s}.txt > $countsdir/counts_${s}_matrix.txt
done
wait

