speciesPrefix=$1
for FILE in ${speciesPrefix}-*_R1_*.fastq
do

echo -e "Element:\t${FILE}"

if [[ $FILE =~ (.*-.*)-(.*)_(.*)_(.*)_(.*)_001\.fastq$ ]]; then
  Species=${BASH_REMATCH[1]}
  SampleID=${BASH_REMATCH[2]}
  S=${BASH_REMATCH[3]}
  L=${BASH_REMATCH[4]}
  Run=${BASH_REMATCH[5]}

## Trimming with fastp
#fastp --correction --cut_right --thread 2 --html ${Species}_${SampleID}_fastp.html --json ${Species}_${SampleID}_fastp.json -l 100 -i ${FILE} -I ${Species}-${SampleID}_${S}_${L}_R2_001.fastq.gz -o ${Species}_${SampleID}_R1.trimmed.fastq -O ${Species}_${SampleID}_R2.trimmed.fastq
# minimum length requirement is specified with -l --length_required --> 100bp

## Merge reads with fastp
#fastp --merge -i ${Species}_${SampleID}_R1.trimmed.fastq -I ${Species}_${SampleID}_R2.trimmed.fastq --merged_out ${Species}_${SampleID}_merged.fastq --detect_adapter_for_pe

#Option 2 - discard reads that can't be merged into dummy .fastq files
fastp --merge -i ${SampleID}_R1.trimmed.fastq -I ${SampleID}_R2.trimmed.fastq --merged_out \
    ${SampleID}_merged.fastq --detect_adapter_for_pe --out1 out1.fastq --out2 out2.fastq \
    --unpaired1 unpaired1.fastq --unpaired2 unpaired2.fastq

fi
done
