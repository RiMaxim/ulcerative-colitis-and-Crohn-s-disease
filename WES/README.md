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

echo "=== 1. Normalize + split multiallelic + compress ==="
bcftools norm -f "$REF" -m -both -Oz -o "${PREFIX}.norm.vcf.gz" "$INPUT_VCF"

bcftools index -t "${PREFIX}.norm.vcf.gz"

echo "=== 2. PASS filter (safe version) ==="
bcftools view -f PASS -Oz -o "${PREFIX}.PASS.vcf.gz" "${PREFIX}.norm.vcf.gz"
bcftools index -t "${PREFIX}.PASS.vcf.gz"

echo "=== 3. DP filter ==="
bcftools view \
-i 'COUNT(FMT/DP>=20)/N_SAMPLES>=0.7' \
-Oz -o "${PREFIX}.PASS.DP20.vcf.gz" \
"${PREFIX}.PASS.vcf.gz"

bcftools index -t "${PREFIX}.PASS.DP20.vcf.gz"

echo "=== 4. Final VEP-ready file ==="
mv "${PREFIX}.PASS.DP20.vcf.gz" "${PREFIX}.ready_for_vep.vcf.gz"
mv "${PREFIX}.PASS.DP20.vcf.gz.tbi" "${PREFIX}.ready_for_vep.vcf.gz.tbi"
```
25) Установка программ для аннотации.
```
echo "=== 1. Deployment of VEP execution environment ==="
docker pull ensemblorg/ensembl-vep:release_115.2

echo "=== 2. Acquisition of pre-indexed VEP cache for GRCh38 ==="
wget https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz

echo "=== 3. ClinVar Database (Medical Classifications) (2026-05-18 11:24, 183M; 2026-05-18 11:24, 595K) ==="
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260517.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260517.vcf.gz.tbi

echo "=== 4. SpliceAI (https://basespace.illumina.com/) ==="
Regisrtaion https://basespace.illumina.com/
wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"
./bs auth
./bs list project
./bs list datasets 66029966
./bs contents dataset -i ds.20a701bc58ab45b59de2576db79ac8d0
./bs download file -i 16534036123 #spliceai_scores.raw.indel.hg38.vcf.gz
./bs download file -i 16534036125 #spliceai_scores.raw.indel.hg38.vcf.gz.tbi
./bs download file -i 16534036127 #spliceai_scores.raw.snv.hg38.vcf.gz
./bs download file -i 16534036128 #spliceai_scores.raw.snv.hg38.vcf.gz.tbi

echo "=== 5. dbNSFP (https://www.dbnsfp.org/download) ==="

wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr1.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr2.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr3.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr4.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr5.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr6.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr7.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr8.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr9.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr10.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr11.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr12.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr13.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr14.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr15.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr16.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr17.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr18.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr19.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr20.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr21.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chr22.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chrX.gz
wget https://zenodo.org/records/14419644/files/dbNSFP4.9a_variant.chrY.gz

zcat dbNSFP4.9a_variant.chr*.gz | bgzip -c > dbNSFP4.9a_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP4.9a_grch38.gz

```
26) Аннотация с помощью VEP.
```
mkdir -p annotation/vep_plugins
git clone https://github.com/Ensembl/VEP_plugins.git annotation/vep_plugins

mkdir -p annotation/vep_data/clinvar
mkdir -p annotation/vep_data/dbNSFP
mkdir -p annotation/vep_data/dbsnp
mkdir -p annotation/vep_data/spliceai

copy files from steps 11-12 and 25.

docker run -it --rm -u $(id -u):$(id -g) -v $(pwd):$(pwd) -w $(pwd) ensemblorg/ensembl-vep:release_115.2 \
vep \
--input_file ./vcf/cohort.ready_for_vep.vcf.gz \
--output_file cohort.annotated.vep.vcf.gz \
--vcf \
--compress_output bgzip \
--assembly GRCh38 \
--species homo_sapiens \
--cache \
--offline \
--dir_cache annotation/vep_cache \
--dir_plugins annotation/vep_plugins \
--fasta reference/Homo_sapiens_assembly38.fasta \
--everything \
--fork 8 \
--force_overwrite \
--custom annotation/vep_data/clinvar/clinvar_20260517.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN \
--custom annotation/vep_data/dbsnp/dbsnp157.vcf.gz,dbSNP,vcf,exact,0,ID \
--plugin dbNSFP,annotation/vep_data/dbNSFP/dbNSFP4.9a_grch38.gz,REVEL_score,CADD_phred,AlphaMissense_score,MPC_score,Polyphen2_HDIV_score,gnomAD_exomes_AF,gnomAD_exomes_NFE_AF \
--plugin SpliceAI,snv=annotation/vep_data/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=annotation/vep_data/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz

#Доступные поля bcftools +split-vep cohort.annotated.vep.vcf.gz -l
#Количество образцов bcftools query -l cohort.annotated.vep.vcf.gz | wc -l
```
27) Индексация cohort.annotated.vep.vcf.gz
```
tabix -s 1 -b 2 -e 2 cohort.annotated.vep.vcf.gz
```
28) Фильтрация по полям с аннотацией.
```
echo "=== 1. Проверить все поля CSQ ==="
bcftools +split-vep cohort.annotated.vep.vcf.gz -l

echo "=== 2. CANONICAL ==="
bcftools +split-vep cohort.annotated.vep.vcf.gz -c CANONICAL:String -i 'CANONICAL="YES"' -Oz -o cohort.canonical.vcf.gz && bcftools index -f cohort.canonical.vcf.gz
bcftools index -n cohort.canonical.vcf.gz

echo "=== 3. MAX_AF ==="
bcftools +split-vep cohort.canonical.vcf.gz -c MAX_AF:Float -i '(MAX_AF < 0.01)' -Oz -o cohort.rare.vcf.gz && bcftools index -f cohort.rare.vcf.gz
bcftools index -n cohort.rare.vcf.gz

echo "=== 4. LoF ==="
bcftools +split-vep cohort.rare.vcf.gz -c IMPACT:String,Consequence:String -i 'IMPACT="HIGH" && Consequence !~ "missense"' -Oz -o cohort.LoF.vcf.gz && bcftools index -f cohort.LoF.vcf.gz && bcftools index -n cohort.LoF.vcf.gz
bcftools filter -e 'INFO/CSQ ~ "missense"' cohort.LoF.vcf.gz -Oz -o cohort.LoF.clean.vcf.gz && mv cohort.LoF.clean.vcf.gz cohort.LoF.vcf.gz && bcftools index -f cohort.LoF.vcf.gz
bcftools index -n cohort.LoF.vcf.gz
bcftools query -f '%INFO/CSQ\n' cohort.LoF.vcf.gz | awk -F'|' '{print $2}' | tr '&' '\n' | sort | uniq -c | sort -nr
bcftools query -f '%INFO/CSQ\n' cohort.LoF.vcf.gz | awk -F'|' '{print $2}' | sort | uniq -c | sort -nr

echo "=== 5. Missense ==="
bcftools +split-vep cohort.rare.vcf.gz -c IMPACT:String,Consequence:String -i 'Consequence ~ "missense_variant" && IMPACT != "HIGH"' -Oz -o cohort.missense.vcf.gz && bcftools index -f cohort.missense.vcf.gz && bcftools index -n cohort.missense.vcf.gz
bcftools filter -e 'INFO/CSQ ~ "stop_gained" || INFO/CSQ ~ "splice_donor" || INFO/CSQ ~ "splice_acceptor" || INFO/CSQ ~ "start_lost"' cohort.missense.vcf.gz -Oz -o cohort.missense.clean.vcf.gz && mv cohort.missense.clean.vcf.gz cohort.missense.vcf.gz && bcftools index -f cohort.missense.vcf.gz
bcftools index -n cohort.missense.vcf.gz
bcftools query -f '%INFO/CSQ\n' cohort.missense.vcf.gz | awk -F'|' '{print $2}' | tr '&' '\n' | sort | uniq -c | sort -nr
bcftools query -f '%INFO/CSQ\n' cohort.missense.vcf.gz | awk -F'|' '{print $2}' | sort | uniq -c | sort -nr

echo "=== 6. Missense Pred ==="
bcftools +split-vep cohort.missense.vcf.gz -s worst -c REVEL_score:Float,AlphaMissense_score:Float,ClinVar_CLNSIG:String -i 'REVEL_score >= 0.5 || AlphaMissense_score >= 0.56 || ClinVar_CLNSIG ~ "Pathogenic|Likely_pathogenic"' -Oz -o cohort.missense.pred.vcf.gz && bcftools index -f cohort.missense.pred.vcf.gz
bcftools index -n cohort.missense.pred.vcf.gz

echo "=== 7. Merge LoF and  Missense Pred ==="
bcftools concat -a cohort.LoF.vcf.gz cohort.missense.pred.vcf.gz -Oz -o cohort.burden.vcf.gz
bcftools index -f cohort.burden.vcf.gz
bcftools index -n cohort.burden.vcf.gz
```
29) Подготовка клинических данных
```
IID	Group	Sex	Age
240125_new_exome_sample1	CD	M	34
240125_new_exome_sample10	UC	M	34
240125_new_exome_sample11	UC	F	70
240125_new_exome_sample12	CD	M	45
240125_new_exome_sample13	Control	F	52
240125_new_exome_sample14	UC	F	27
240125_new_exome_sample15	UC	F	41
240125_new_exome_sample16	CD	M	53
240125_new_exome_sample17	CD	M	45
240125_new_exome_sample18	Control	F	20
240125_new_exome_sample19	UC	F	50
240125_new_exome_sample2	Control	M	22
240125_new_exome_sample20	UC	M	59
240125_new_exome_sample21	Control	F	42
240125_new_exome_sample3	UC	F	27
240125_new_exome_sample4	UC	M	58
240125_new_exome_sample5	CD	M	31
240125_new_exome_sample6	CD	M	45
240125_new_exome_sample7	UC	F	19
240125_new_exome_sample8	UC	F	43
240125_new_exome_sample9	CD	M	25
barcode1	Control	M	47
barcode13	CD	M	30
barcode14	Control	F	55
barcode15	CD	F	44
barcode16	CD	F	32
barcode2	CD	F	26
barcode3	CD	F	33
barcode4	UC	M	27
barcode41	UC	M	26
barcode42	CD	F	56
barcode43	UC	F	24
barcode44	UC	M	36
barcode45	Control	F	73
barcode46	UC	M	50
barcode47	UC	F	20
barcode48	CD	M	20
barcode57	UC	F	70
barcode58	CD	M	43
barcode59	UC	M	38
barcode60	UC	M	55
barcode61	CD	F	24
barcode62	CD	F	29
barcode63	CD	M	45
barcode64	CD	M	46
barcode65	UC	F	20
barcode66	UC	F	30
barcode67	CD	M	60
barcode68	CD	F	28
barcode69	CD	F	30
barcode70	CD	F	23
barcode71	CD	F	37
barcode72	CD	M	22
barcode73	CD	M	42
barcode74	UC	M	35
barcode75	CD	F	19
barcode76	UC	F	51
barcode77	CD	M	38
barcode78	Control	M	47
barcode79	UC	M	60
barcode80	CD	F	31
barcode81	UC	M	58
barcode82	CD	M	38
barcode83	UC	M	36
barcode84	CD	M	40
barcode85	CD	F	21
barcode86	CD	F	56
barcode87	UC	M	42
barcode88	UC	F	31
barcode89	Control	F	22
barcode90	CD	M	21
barcode91	CD	M	57
barcode92	CD	M	44
barcode_1	Control	F	38
barcode_13	UC	M	46
barcode_14	CD	M	24
barcode_15	UC	M	65
barcode_16	CD	M	24
barcode_2	CD	M	61
barcode_3	UC	M	56
barcode_4	UC	F	37
barcode_41	CD	F	42
barcode_42	UC	F	45
barcode_43	CD	M	32
barcode_44	CD	M	20
barcode_45	UC	F	24
barcode_46	UC	F	42
barcode_47	CD	F	38
barcode_48	CD	F	36
barcode_57	CD	F	27
barcode_58	Control	F	21
barcode_59	UC	M	37
barcode_60	UC	F	38
barcode_61	Control	M	39
barcode_62	Control	F	44
barcode_63	Control	M	18
barcode_64	Control	M	18
barcode_65	UC	M	42
barcode_66	CD	F	57
barcode_67	UC	M	30
barcode_68	UC	F	67
barcode_69	CD	F	38
barcode_70	CD	F	57
barcode_71	CD	M	62
barcode_72	CD	F	21
barcode_73	UC	F	32
barcode_74	UC	M	43
barcode_75	CD	M	59
barcode_76	CD	F	51
barcode_77	UC	M	29
barcode_78	UC	M	67
barcode_79	UC	F	43
barcode_89	CD	M	36
barcode_90	CD	F	33
barcode_91	UC	M	42
barcode_92	UC	F	38
barcode_93	CD	F	26
barcode_94	CD	F	49
barcode_95	UC	M	25
```
Проверяем соответствие образцов:
```
cut -f1 phenotype.tsv | tail -n +2 | sort > pheno_ids.txt
bcftools query -l cohort.annotated.vep.vcf.gz | sort > vcf_ids.txt
comm -23 pheno_ids.txt vcf_ids.txt
```
30) Burden analysis. Входные файлы - cohort.burden.vcf.gz (шаг 28) и metadata_WES.txt (шаг 29). В данном коде сравнивается UC и Control.
```
library(SKAT)
library(ggplot2)
library(reshape2)

vcf_in        <- "./cohort.burden.vcf.gz"
vcf_fixed     <- "./cohort.burden.id_fixed.vcf.gz"
plink_prefix  <- ".cohort_burden_plink"
metadata_file <- "./metadata_WES.txt"
bcftools_path <- "./miniforge3/bin/bcftools"
plink_path <- "./miniforge3/bin/plink"

cmd_id <- paste0(bcftools_path, " annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' ", vcf_in, " -Oz -o ", vcf_fixed)
system(cmd_id)
cmd_index <- paste0(bcftools_path, " index -f ", vcf_fixed)
system(cmd_index)
cmd_plink <- paste0(plink_path, " --vcf ", vcf_fixed, " --make-bed --out ", plink_prefix, " --allow-extra-chr --double-id")
system(cmd_plink)

###

cmd_raw_csq <- paste0(bcftools_path, " query -f '%CHROM\\_%POS\\_%REF\\_%ALT\\t%INFO/CSQ\\n' ", vcf_fixed, " > raw_csq.tmp")
system(cmd_raw_csq)

raw_data <- read.table("raw_csq.tmp", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(raw_data) <- c("Variant_ID", "CSQ_String")

extract_gene_name <- function(csq_str) {
  if (is.na(csq_str) || csq_str == "") return(NA)
  first_transcript <- strsplit(csq_str, ",")[[1]][1]
  fields <- strsplit(first_transcript, "\\|")[[1]]
  if (length(fields) >= 4) {
    gene_symbol <- fields[4]
    if (is.na(gene_symbol) || gene_symbol == "" || gene_symbol == "-") {
      if (length(fields) >= 5) gene_symbol <- fields[5] else gene_symbol <- NA
    }
    return(gene_symbol)
  }
  return(NA)
}


genes_list <- sapply(raw_data$CSQ_String, extract_gene_name)
setid_final <- data.frame(Gene = genes_list, Variant_ID = raw_data$Variant_ID, stringsAsFactors = FALSE)
setid_final <- setid_final[!is.na(setid_final$Gene) & setid_final$Gene != "" & setid_final$Gene != "-", ]
cat(sprintf("Найдено %d уникальных генов\n", length(unique(setid_final$Gene))))
write.table(setid_final, file = "./cohort_burden.SetID", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
file.remove("raw_csq.tmp")

###

File.Bed  <- paste0(plink_prefix, ".bed")
File.Bim  <- paste0(plink_prefix, ".bim")
File.Fam  <- paste0(plink_prefix, ".fam")
File.SetID <- "Desktop/WORK/gut/2_stage/R/cohort_burden.SetID"
File.SSD  <- "Desktop/WORK/gut/2_stage/R/cohort_burden.SSD"
File.Info <- "Desktop/WORK/gut/2_stage/R/cohort_burden.SSD.info"

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)

###

metadata <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata <- metadata[metadata$Group %in% c("UC", "Control"), ]
metadata$Phenotype <- ifelse(metadata$Group == "Control", 0, 1)

fam <- read.table(File.Fam, header = FALSE, stringsAsFactors = FALSE)
colnames(fam) <- c("FID", "IID", "Father", "Mother", "Sex_Plink", "Pheno_Plink")

metadata_sorted <- metadata[match(fam$IID, metadata$IID), ]
samples_to_include <- !is.na(metadata_sorted$Phenotype)

n_uc <- sum(metadata_sorted$Phenotype == 1 & samples_to_include, na.rm = TRUE)
n_control <- sum(metadata_sorted$Phenotype == 0 & samples_to_include, na.rm = TRUE)

cat(sprintf("  - UC: %d\n", n_uc))
cat(sprintf("  - Control: %d\n", n_control))

###

SSD.INFO <- Open_SSD(File.SSD, File.Info)

obj_uc <- SKAT_Null_Model(Phenotype ~ 1, 
                          data = metadata_sorted, 
                          out_type = "D", 
                          Adjustment = TRUE)

n_genes <- nrow(SSD.INFO$SetInfo)
gene_names <- SSD.INFO$SetInfo$SetID
n_markers <- SSD.INFO$SetInfo$SetSize

p_burden <- rep(NA, n_genes)
p_skato <- rep(NA, n_genes)
uc_carriers <- rep(0, n_genes)
control_carriers <- rep(0, n_genes)
direction_vec <- rep(NA, n_genes)

for (i in 1:n_genes) {
  if (i %% 50 == 0) cat(sprintf("Processed %d of %d genes (%.1f%%)\n", 
                                i, n_genes, 100 * i / n_genes))
  
  genotypes <- tryCatch({
    Get_Genotypes_SSD(SSD.INFO, i)
  }, error = function(e) NULL)
  
  if (is.null(genotypes)) next
  
  genotypes <- as.matrix(genotypes)
  if (ncol(genotypes) == 0) next
  
  genotypes <- genotypes[samples_to_include, , drop = FALSE]
  phenotypes <- metadata_sorted$Phenotype[samples_to_include]
  
  has_mut <- rowSums(genotypes > 0, na.rm = TRUE) > 0
  if (sum(has_mut, na.rm = TRUE) < 2) next
  
  uc_carriers[i] <- sum(has_mut & phenotypes == 1, na.rm = TRUE)
  control_carriers[i] <- sum(has_mut & phenotypes == 0, na.rm = TRUE)
  
  burden_score <- rowSums(genotypes > 0, na.rm = TRUE)
  if (var(burden_score) > 0) {
    model <- tryCatch({
      glm(phenotypes ~ burden_score, family = binomial())
    }, error = function(e) NULL)
    
    if (!is.null(model)) {
      p_burden[i] <- summary(model)$coefficients[2, 4]
      direction_vec[i] <- ifelse(coef(model)[2] > 0, "Risk", "Protective")
    }
  }
  
  tryCatch({
    obj_skat <- SKAT_Null_Model(phenotypes ~ 1, out_type = "D")
    skat_result <- SKAT(genotypes, obj_skat, method = "SKATO")
    p_skato[i] <- skat_result$p.value
  }, error = function(e) {})
}

Close_SSD()

### 

combined_results <- data.frame(
  Gene = gene_names,
  N_Variants = n_markers,
  P_Burden = p_burden,
  P_SKATO = p_skato,
  UC_Carriers = uc_carriers,
  UC_Total = n_uc,
  Control_Carriers = control_carriers,
  Control_Total = n_control,
  Direction = direction_vec,
  stringsAsFactors = FALSE
)

before <- nrow(combined_results)
combined_results <- combined_results[!is.na(combined_results$P_SKATO), ]
after <- nrow(combined_results)
cat(sprintf("Genes removed due to NA in SKATO: %d\n", before - after))
cat(sprintf("Genes remaining: %d\n", after))

combined_results$Freq_UC <- combined_results$UC_Carriers / combined_results$UC_Total
combined_results$Freq_Control <- combined_results$Control_Carriers / combined_results$Control_Total

combined_results$FDR_Burden <- p.adjust(combined_results$P_Burden, method = "BH")
combined_results$FDR_SKATO <- p.adjust(combined_results$P_SKATO, method = "BH")

combined_results$Interpretation <- "NO_ASSOCIATION"
for (i in 1:nrow(combined_results)) {
  if (!is.na(combined_results$FDR_SKATO[i]) && !is.na(combined_results$FDR_Burden[i])) {
    if (combined_results$FDR_SKATO[i] < 0.05 & combined_results$FDR_Burden[i] < 0.05) {
      if (!is.na(combined_results$Direction[i])) {
        combined_results$Interpretation[i] <- paste0("UNIDIRECTIONAL_", combined_results$Direction[i])
      }
    } else if (combined_results$FDR_SKATO[i] < 0.05 & combined_results$FDR_Burden[i] >= 0.05) {
      combined_results$Interpretation[i] <- "HETEROGENEOUS"
    } else if (combined_results$P_SKATO[i] < 0.05 & combined_results$P_Burden[i] < 0.05) {
      combined_results$Interpretation[i] <- "NOMINAL_ONLY"
    }
  }
}

combined_results <- combined_results[order(combined_results$P_SKATO), ]

write.table(combined_results, "./statistics_UC_vs_Control.txt", row.names = FALSE, sep ="\t")

###
```
