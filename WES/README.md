![Logo](Pipeline.png)
**Конвейер по обработке данных**

> $1 - ID входного файла (например Patient1 из файлов Patient1_1.fastq.gz и Patient1_2.fastq.gz)
> $2 - количество потоков (threads)

1) Проверка качества с помощью FastQC (версия v0.11.5)
```
fastqc -t $2 $1_1.fastq.gz $1_2.fastq.gz
```
2) Фильтрация данных с помощью trimmomatic (версия 0.39)
```
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar \
 PE \
 -threads $2 \
 -phred33 \
 $1_1.fastq.gz $1_2.fastq.gz \
 $1_R1_paired.fastq.gz $1_R1_unpaired.fastq.gz \
 $1_R2_paired.fastq.gz $1_R2_unpaired.fastq.gz \
 ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/All_adapters.fa:2:30:10 \
 LEADING:20 \
 TRAILING:20 \
 SLIDINGWINDOW:4:20 \
 MINLEN:50
```
3) Скачивание эталонного генома человека (сборка GATK Homo_sapiens_assembly38.fasta) и словаря последовательности (Homo_sapiens_assembly38.dict) для эталонного генома (требуется 1 раз)
```
wget https://rcs.bu.edu/examples/bioinformatics/gatk/ref/Homo_sapiens_assembly38.dict
wget https://rcs.bu.edu/examples/bioinformatics/gatk/ref/Homo_sapiens_assembly38.fasta
```
4) Индексация генома с помощью bwa-mem2 (версия 2.2.1) (требуется 1 раз)
```
/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index Homo_sapiens_assembly38.fasta
```
5) Выравнивание парных последовательностей на эталонный геном с формированием группы чтения (RG) и последующей сортировкой необходимой для 7 шага.
```
SAMPLE_ID=$1
RG=$"@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1"

/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem \
 -t $2 \
 -M -Y \
 -R "$RG" \
 ./reference/Homo_sapiens_assembly38.fasta.gz \
 ${SAMPLE_ID}_R1_paired.fastq.gz \
 ${SAMPLE_ID}_R2_paired.fastq.gz | \
 samtools sort -@ $2 -o ./bam/${SAMPLE_ID}.bam -
```
6) Загрузка GATK образа  (требуется 1 раз)
```
docker pull broadinstitute/gatk:4.6.2.0
```
7) Запуск инструмента MarkDuplicates из GATK внутри Docker-контейнера для обнаружения и маркировки ПЦР-дупликатов в BAM-файле.
```
docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk MarkDuplicates \
 -I ./bam/$1.bam \
 -O ./bam/$1.marked.bam \
 -M ./bam/$1.metrics.txt \
 --CREATE_INDEX true
```
8) Скачивание панели Vazyme VAHTS Target Capture Core Exome Panel (https://www.vazymeglobal.com/product-center/capture-probe/vahts-target-capture-core-exome-panel). Внутри архива — 4 файла. Для работы используется файл CoreExomePanel.hg38.p12.target.v3(1).bed. Перед использованием его необходимо один раз переименовать в CoreExomePanel.hg38.p12.target.v3.bed (для шаге 9).

9) Вычисление глубины покрытия без дубликатов. 34128352 - суммарная длина панели
```
samtools depth $1.bam -b CoreExomePanel.hg38.p12.target.v3.bed > $1.bed

depth=30
GENOME_SIZE=34128352

awk -v d=$depth -v g=$GENOME_SIZE '
    $3 >= d { count++ }
    END { printf "%.5f%%\n", (count/g)*100 }
' $1.bed

```
10) Загрузка SNP в VCF-формате (версия от 2025-01-15 21:27, 28 ГБ)
```
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
```
11) Замена RefSeq ID в первом столбце на формат chr (например, NC_000001.11 → chr1). Добавление информации о длине хромосом в заголовок файла (например, ##contig=<ID=chr1,length=248956422>). Итоговое название файла dbsnp157.vcf.gz

12) Запуск инструмента IndexFeatureFile из GATK внутри Docker-контейнера для создания индекса у VCF-файла с SNP-датасетом dbsnp157. Выходной файл - dbsnp157.vcf.gz.tbi
```
docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk IndexFeatureFile -I dbsnp157.vcf.gz
```
13) Загрузка InDels в VCF-формате (версия от 2021-05-04 16:17, 20M и 2021-05-04 16:17, 59M) и индексы файлов. Mills_and_1000G_gold_standard.indels.hg38.vcf.gz — это качественный «золотой стандарт». Используется, когда нужно быть максимально уверенным в валидации. Homo_sapiens_assembly38.known_indels.vcf.gz — это максимально полный справочник. Используется для широкого «маскирования» вариабельных участков.
```
wget https://rcs.bu.edu/examples/bioinformatics/gatk/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://rcs.bu.edu/examples/bioinformatics/gatk/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

wget https://rcs.bu.edu/examples/bioinformatics/gatk/ref/Homo_sapiens_assembly38.known_indels.vcf.gz
wget https://rcs.bu.edu/examples/bioinformatics/gatk/ref/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
```
14) Создание индекса fai для генома
```
samtools faidx Homo_sapiens_assembly38.fasta
```
15) Конвертация BED-файла в IntervalList (формат GATK).
```
docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk BedToIntervalList \
 -I CoreExomePanel.hg38.p12.target.v3.bed \
 -O CoreExomePanel.hg38.p12.target.v3.interval_list \
 -R ./reference/Homo_sapiens_assembly38.fasta \
 --SEQUENCE_DICTIONARY ./reference/Homo_sapiens_assembly38.dict \
 --DROP_MISSING_CONTIGS true \
 --SORT true
```
16) Запуск BaseRecalibrator для перекалибровки базовых качеств (Base Quality Score Recalibration, BQSR).
```
docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk BaseRecalibrator \
 -I ./bam/$1.marked.bam \
 -R ./reference/Homo_sapiens_assembly38.fasta \
 --known-sites dbsnp157.vcf.gz \
 --known-sites Homo_sapiens_assembly38.known_indels.vcf.gz \
 --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
 -L CoreExomePanel.hg38.p12.target.v3.interval_list \
 --interval-padding 100 \
 -O ./bam/$1.recal.table
```
17) Запуск ApplyBQSR применяет таблицу перекалибровки качеств (созданную на предыдущем шаге $1.recal.table) к BAM файлу (созданный на шаге 9 на$1.marked.bam), исправляя оценки качества каждого основания.
```
docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk ApplyBQSR \
 -I ./bam/$1.marked.bam \
 -R ./reference/Homo_sapiens_assembly38.fasta \
 --bqsr-recal-file ./bam/$1.recal.table \
 --create-output-bam-index true \
 --add-output-sam-program-record true \
 -O ./bam/$1.recal.bam
```
18) Запуск HaplotypeCaller. Индивидуальный Variant Calling (GATK)
```
docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk HaplotypeCaller \
 -R ./reference/Homo_sapiens_assembly38.fasta \
 -I ./bam/$1.recal.bam \
 -O ./gvcf/$1.g.vcf.gz \
 -L CoreExomePanel.hg38.p12.target.v3.interval_list \
 -ERC GVCF
```
19) Когортный анализ. Создание карты образцов (sample_map.txt). Используется в шаге 20.
```
Patient1 ./gvcf/barcode1.g.vcf.gz
Patient2 ./gvcf/barcode2.g.vcf.gz
...
Patient126 ./gvcf/barcode126.g.vcf.gz
```
20) Когортный анализ. Запуск GenomicsDBImport для создание базы данных. Ниже пример для хромосомы 1 (повторить для всех хромосом chr1-22, X, Y).
```
mkdir -p genomicsdb

CHRS=( {1..22} X Y )

for chr in "${CHRS[@]}"; do
    CHR_NAME="chr${chr}"
    WORKSPACE="./genomicsdb/${CHR_NAME}_db"
    
    echo "========================================"
    echo "Starting import for ${CHR_NAME}..."
    echo "========================================"

    rm -rf "$WORKSPACE"

    docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
     gatk GenomicsDBImport \
     --genomicsdb-workspace-path "$WORKSPACE" \
     -R ./reference/Homo_sapiens_assembly38.fasta \
     --sample-name-map sample_map.txt \
     -L "$CHR_NAME"
done
# docker run --rm -v $(pwd):/data -w /data alpine rm -rf ./genomicsdb
```
21) Когортный анализ. Запуск GenotypeGVCFs для генотипирования. Ниже пример для хромосомы 1 (повторить для всех хромосом chr1-22, X, Y).
```
mkdir -p vcf
CHRS=( {1..22} X Y )

for chr in "${CHRS[@]}"; do
    CHR_NAME="chr${chr}"
    DB_PATH="./genomicsdb/${CHR_NAME}_db"
    OUT_VCF="./vcf/cohort_${CHR_NAME}.vcf.gz"
    
    echo "========================================"
    echo "Starting import for ${CHR_NAME}..."
    echo "========================================"

    docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
     gatk GenotypeGVCFs \
     -R ./reference/Homo_sapiens_assembly38.fasta \
     -V "gendb://$DB_PATH" \
     -O "$OUT_VCF"
done
```
22) Когортный анализ. Объединение хромосом в один файл.
```
docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk MergeVcfs \
 -I ./vcf/cohort_chr1.vcf.gz \
 -I ./vcf/cohort_chr2.vcf.gz \
 -I ./vcf/cohort_chr3.vcf.gz \
 -I ./vcf/cohort_chr4.vcf.gz \
 -I ./vcf/cohort_chr5.vcf.gz \
 -I ./vcf/cohort_chr6.vcf.gz \
 -I ./vcf/cohort_chr7.vcf.gz \
 -I ./vcf/cohort_chr8.vcf.gz \
 -I ./vcf/cohort_chr9.vcf.gz \
 -I ./vcf/cohort_chr10.vcf.gz \
 -I ./vcf/cohort_chr11.vcf.gz \
 -I ./vcf/cohort_chr12.vcf.gz \
 -I ./vcf/cohort_chr13.vcf.gz \
 -I ./vcf/cohort_chr14.vcf.gz \
 -I ./vcf/cohort_chr15.vcf.gz \
 -I ./vcf/cohort_chr16.vcf.gz \
 -I ./vcf/cohort_chr17.vcf.gz \
 -I ./vcf/cohort_chr18.vcf.gz \
 -I ./vcf/cohort_chr19.vcf.gz \
 -I ./vcf/cohort_chr20.vcf.gz \
 -I ./vcf/cohort_chr21.vcf.gz \
 -I ./vcf/cohort_chr22.vcf.gz \
 -I ./vcf/cohort_chrX.vcf.gz \
 -I ./vcf/cohort_chrY.vcf.gz \
 -O ./vcf/cohort_raw.vcf.gz

```
23) Разделение исходного VCF на SNPs и Indels. Фильтрация и последующее объединение.
```
docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk SelectVariants \
 -R ./reference/Homo_sapiens_assembly38.fasta \
 -V ./vcf/cohort_raw.vcf.gz \
 -select-type SNP \
 -O ./vcf/cohort_snps.vcf.gz

docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk SelectVariants \
 -R ./reference/Homo_sapiens_assembly38.fasta \
 -V  ./vcf/cohort_raw.vcf.gz \
 -select-type INDEL \
 -select-type MIXED \
 -O ./vcf/cohort_indels.vcf.gz

docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk VariantFiltration \
 -R ./reference/Homo_sapiens_assembly38.fasta \
 -V ./vcf/cohort_snps.vcf.gz \
 --filter-expression "QUAL < 30.00" --filter-name "HDFLT_QUAL" \
 --filter-expression "QD < 2.00" --filter-name "HDFLT_QD" \
 --filter-expression "FS > 60.00" --filter-name "HDFLT_FS" \
 --filter-expression "SOR > 3.00" --filter-name "HDFLT_SOR" \
 --filter-expression "MQ < 40.00" --filter-name "HDFLT_MQ" \
 --filter-expression "MQRankSum < -12.50" --filter-name "HDFLT_MQRankSum" \
 --filter-expression "ReadPosRankSum < -8.00" --filter-name "HDFLT_ReadPosRankSum" \
 -O ./vcf/cohort_snps_filtered.vcf.gz

docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk VariantFiltration \
 -R ./reference/Homo_sapiens_assembly38.fasta \
 -V ./vcf/cohort_indels.vcf.gz \
 --filter-expression "QUAL < 30.00" --filter-name "HDFLT_INDEL_QUAL" \
 --filter-expression "QD < 2.00" --filter-name "HDFLT_INDEL_QD" \
 --filter-expression "FS > 200.00" --filter-name "HDFLT_INDEL_FS" \
 --filter-expression "SOR > 10.00" --filter-name "HDFLT_INDEL_SOR" \
 --filter-expression "ReadPosRankSum < -20.00" --filter-name "HDFLT_INDEL_ReadPosRankSum" \
 -O ./vcf/cohort_indels_filtered.vcf.gz

docker run --rm -v $(pwd):/data -w /data broadinstitute/gatk:4.6.2.0 \
 gatk MergeVcfs \
 -I ./vcf/cohort_snps_filtered.vcf.gz \
 -I ./vcf/cohort_indels_filtered.vcf.gz \
 -O ./vcf/cohort_filtered.vcf.gz

rm -f ./vcf/cohort_snps.vcf.gz ./vcf/cohort_snps.vcf.gz.tbi
rm -f ./vcf/cohort_indels.vcf.gz ./vcf/cohort_indels.vcf.gz.tbi
rm -f ./vcf/cohort_snps_filtered.vcf.gz ./vcf/cohort_snps_filtered.vcf.gz.tbi
rm -f ./vcf/cohort_indels_filtered.vcf.gz ./vcf/cohort_indels_filtered.vcf.gz.tbi

#bcftools view -H -v snps -f PASS ./vcf/cohort_filtered.vcf.gz | wc -l
```
24) Финальная фильтрация
```
REF="./reference/Homo_sapiens_assembly38.fasta"
INPUT_VCF="./vcf/cohort_filtered.vcf.gz"
PREFIX="./vcf/cohort"

echo "=== 1. Normalization and multiallelic splitting ==="
bcftools norm -f "$REF" -m -both -O v "$INPUT_VCF" > "${PREFIX}.norm.vcf"

echo "=== 2. VCF Sorting ==="
bcftools sort "${PREFIX}.norm.vcf" > "${PREFIX}.sorted.vcf"

echo "=== 3. Filtration: PASS only ==="
bcftools view -i 'FILTER="PASS"' "${PREFIX}.sorted.vcf" > "${PREFIX}.PASS.vcf"

echo "=== 4. Filtration by Call Rate: DP >= 20x in at least 70% of samples ==="
bcftools view -i 'COUNT(FMT/DP>=20)/N_SAMPLES >= 0.7' "${PREFIX}.PASS.vcf" > "${PREFIX}.PASS.DP20.vcf"

echo "=== 5. Final sorting and compression to .gz ==="
bcftools sort "${PREFIX}.PASS.DP20.vcf" -O z > "${PREFIX}.ready_for_vep.vcf.gz"
bcftools index -t "${PREFIX}.ready_for_vep.vcf.gz"

echo "=== 6. Cleaning up temporary files ==="
rm -f "${PREFIX}.norm.vcf" "${PREFIX}.sorted.vcf" "${PREFIX}.PASS.vcf" "${PREFIX}.PASS.DP20.vcf"
```
25) Скачиваем базовый кэш Ensembl v115.2 (Человек, GRCh38) c https://hub.docker.com/r/ensemblorg/ensembl-vep 
```
mkdir -p annotation
mkdir -p annotation/vep_cache annotation/vep_plugins annotation/vep_data

docker run --rm -it \
  -v $(pwd)/annotation/vep_cache:/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_115.2 \
  perl /opt/vep/src/ensembl-vep/INSTALL.pl --AUTO c --SPECIES homo_sapiens --ASSEMBLY GRCh38 --DESTDIR /opt/vep/.vep
```
25) Подготовка клинических данных
26) Конвертация VCF в формат PLINK
27) Запуск статистического теста (Логистическая регрессия)
