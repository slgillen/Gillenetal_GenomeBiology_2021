#directories
inputdir='/directoryoffinputfiles'
datadir='/directoryforoutputs'

#run separately for diff read lengths if want downstream Q.C for each read length
lengths='27,28,29,30,31' #change this depending on the length range you want to look at
for length in ${lengths//,/ }
do
  for s in ${inputfiles//,/ }
  do
    python countingScript.py -b $datadir/sorted_RPFaligned_${s}.bam -f $indir/gencode_v28_proteincoding_mostabundanttranscriptHEK293.fa -l ${length} -o $datadir/CountsOutput -of RPFcountsTable_${s}_${length}.txt
  done
done
wait

