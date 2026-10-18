#!/bin/bash
# Упрощённый запуск hlahd для набора образцов.
#
# Использование:
#   ./run_hlahd.sh <output_dir> <fastq_list.txt>
#
# fastq_list.txt — список fastq-файлов, по два подряд на пару (R1, R2),
# по одному пути на строку. Для каждой пары определяется имя образца
# (по общей части имени файлов) и запускается hlahd. Результат каждого
# образца кладётся в <output_dir>/<sample_name>/.
 
set -euo pipefail
 
THREADS="${THREADS:-28}"
HLAHD_BASE="${HLAHD_BASE:-/opt/hlahd.1.7.0}"
DICT=${HLAHD_BASE}/dictionary
HLA=${HLAHD_BASE}/HLA_gene.split.3.50.0.txt
 
OUT_DIR=$1
LIST=$2
 
mkdir -p "$OUT_DIR"
 
R1=""
while read -r FASTQ; do
    [ -z "$FASTQ" ] && continue
 
    if [ -z "$R1" ]; then
        R1=$FASTQ
        continue
    fi
    R2=$FASTQ
 
    # имя образца — общий префикс файлов до _R1/_R2
    SAMPLE=$(basename "$R1" | sed -E 's/(_R1|\.R1|_1)(\.fastq|\.fq)(\.gz)?$//')
 
    SAMPLE_DIR="$OUT_DIR/$SAMPLE"
    LOG="$SAMPLE_DIR/hlahd-log.txt"
 
    echo "== [$SAMPLE] HLA-typing: $R1 / $R2"
    mkdir -p "$SAMPLE_DIR"
 
    hlahd.sh -t $THREADS -m 50 \
        "$R1" "$R2" \
        "$HLA" "$DICT" "result" "$SAMPLE_DIR" > "$LOG"
 
    cp "$SAMPLE_DIR/result/result/result_final.result.txt" "$SAMPLE_DIR/${SAMPLE}_HLAHD_result.txt"
    rm -rf "$SAMPLE_DIR/result"
 
    echo "== [$SAMPLE] done -> $SAMPLE_DIR/${SAMPLE}_HLAHD_result.txt"
 
    R1=""
done < "$LIST"
 