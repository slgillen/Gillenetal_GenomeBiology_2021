#!/bin/bash

#library preparation for small RNA samples done with NextFlex small RNA seq kit

#directories
inputdir='/directoryoffinputfiles'
datadir='/directoryforoutputs'

inputfiles='RPF_con1, RPF_con2, RPF_con3'

fastqcdir=$datadir/initial_fastqc
mkdir $fastqcdir

#quality check
for s in ${inputfiles//,/ }
do	
  fastqc $datadir/${s}.fastq --extract --outdir=$fastqcdir &
done
wait

#remove adapters, quality trimming, size selection (25 to 35nt reads + 8nt UMI)
for s in ${inputfiles//,/ }
do
  cutadapt --nextseq-trim=20 -a TGGAATTCTCGGGTGCCAAGG -m 33 -M 43 -o $datadir/${s}_cutadapt.fq $datadir/${s}_cutadapt.fq &
done
wait

#quality check
for s in ${inputfiles//,/ }
do	
  fastqc $datadir/${s}_cutadapt.fq --extract --outdir=$fastqcdir &
done
wait

#read deduplication based on unique molecular indexes (UMI)
for s in ${inputfiles//,/ }
do
  cd-hit-dup -i $datadir/${s}_fastp_cutadapt.fq -o $datadir/${s}_cdhitdup.fq -e 0
done
wait

#check reads remaining
for s in ${inputfiles//,/ }
do	
  fastqc $datadir/${s}_cdhitdup.fq --extract --outdir=$fastqcdir &
done
wait

#remove UMI sequences after deduplication (4nt either side of read)
for s in ${inputfiles//,/ }
do
  cutadapt -u 4 -u -4 -o $datadir/${s}_trimUMI.fq $datadir/${s}_cdhitdup.fq &
done
wait

rRNAdir=$datadir/rRNA_alignments
mkdir $rRNAdir
#Align to rRNA sequences to remove contaminant reads
for s in ${inputfiles//,/ }
do
  (bowtie -S --threads 12 -q -n 1 -l 22 --best --norc --un $rRNAdir/not_rRNA_${s}.fq --al $rRNAdir/rRNA_${s}.fq $inputdir/rRNAhuman $datadir/${s}_trimUMI.fq $rRNAdir/aligned_${s}.sam) 2> $rRNAdir/bowtie_stats_rRNA_${s}.txt
  rm $rRNAdir/aligned_${s}.sam
done
wait

#fastqc of non-rRNA reads
fastqcdir=$datadir/fastqc_after_rRNA_removal
mkdir $fastqcdir
for s in ${inputfiles//,/ }
do	
 fastqc $rRNAdir/not_rRNA_${s}.fq --extract --outdir=$fastqcdir &
done
wait

#Align to protein coding transcriptome - the most abundant transcript per gene in these cells was selected using the matched total RNA-seq data
aligndir=$datadir/transcriptome_alignments
mkdir $aligndir
for s in ${inputfiles//,/ }
do
  bowtie -S --threads 12 -q -n 2 -l 22 --best --norc --un $aligndir/not_aligned_${s}.fq --al $aligndir/RPF_aligned_${s}.fq $inputdir/gencode_v28_proteincoding_mostabundanttranscriptHEK293.fa $rRNAdir/not_rRNA_${s}.fq $aligndir/RPFaligned_${s}.sam) 2> $aligndir/bowtie_alignment_stats_${s}.txt
done
wait

#sam to bam
for s in ${inputfiles//,/ }; do
  samtools view -bS $aligndir/RPFaligned_${s}.sam > $aligndir/RPFaligned_${s}.bam & 
done
wait

#sort bam
for s in ${inputfiles//,/ }; do
  samtools sort $aligndir/RPFaligned_${s}.bam -o $aligndir/sorted_RPFaligned_${s}.bam -@ 10 -m 1G 
done
wait

#index bam
for s in ${inputfiles//,/ }; do
  samtools index $aligndir/sorted_RPFaligned_${s}.bam $aligndir/sorted_RPFaligned_${s}.bai & 
done
wait

#get lengths of aligned reads
for s in ${inputfiles//,/ }
do
    samtools view $aligndir/sorted_RPFaligned_${s}.bam | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c >> $aligndir/lengths_${s}.txt
done
wait

#custom python script (adapted from RiboCount) to get position of read starts along the transcript 
#can be done for all length or specified read lengths 
#offsets can be applied at this stage 
#e.g. -l 28,29,30 -s 12,12,12 would only output data for read lengths 28,29,30 and would apply a P-site offset of 12

for s in ${inputfiles//,/ }; do
  python countingScript.py -b $datadir/sorted_RPFaligned_${s}.bam -f $indir/gencode_v28_proteincoding_mostabundanttranscriptHEK293.fa -o $datadir/CountsOutput -of RPFcountsTable_${s}.txt
done
wait
