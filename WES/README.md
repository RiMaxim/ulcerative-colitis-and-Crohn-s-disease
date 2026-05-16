**Конвейр по обработке данных**

1) fastqc -t $2

java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE \
-threads 64 \
-phred33 \
$1_1.fastq.gz $1_2.fastq.gz \
$1_R1_paired.fastq.gz $1_R1_unpaired.fastq.gz \
$1_R2_paired.fastq.gz $1_R2_unpaired.fastq.gz \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/All_adapters.fa:2:30:10 \
LEADING:20 \
TRAILING:20 \
SLIDINGWINDOW:4:20 \
MINLEN:50

SAMPLE_ID=$1

RG=$"@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1"

/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem \
  -t 64 \
  -M -Y \
  -R "$RG" \
  ./reference/Homo_sapiens_assembly38.fasta.gz \
  ${SAMPLE_ID}_R1_paired.fastq.gz \
  ${SAMPLE_ID}_R2_paired.fastq.gz | \
  samtools sort -@ 64 -o ./bam/${SAMPLE_ID}.bam -

samtools indeх ./bam/${SAMPLE_ID}.bam -@ 60

./run2.sh ${SAMPLE_ID}













 1Help            2Save            3Mark            4Replac          5Copy            6Move            7Search          8Delete          9PullDn         10Quit
