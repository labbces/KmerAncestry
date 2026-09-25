#!/usr/bin/env bash
set -euo pipefail

PATH=/home/diriano/FASTK/:$PATH
export PATH

K=19
SP=spontaneum
THREADS=10
MEM_GB=15
OUTDIR="FastK_${SP}"

SRR_LIST="$SP.srrs.txt"

mkdir -p "$OUTDIR"
PARTS=()

while IFS= read -r SRR || [[ -n "$SRR" ]]; do
    [[ -z "$SRR" || "$SRR" == \#* ]] && continue

    R1="/DataBig/SRR/NCBI/sra_datasets/${SRR}_1.fastq.gz"
    R2="/DataBig/SRR/NCBI/sra_datasets/${SRR}_2.fastq.gz"
    OUT="${OUTDIR}/${SP}.${SRR}.${K}"

    [[ -s "$R1" && -s "$R2" ]] || {
        echo "Missing or empty FASTQ pair for ${SRR}" >&2
        exit 1
    }

    if [[ -e "${OUT}.done" && -s "${OUT}.ktab" ]]; then
	    echo "Skipping ${SRR}: completed table already exists"
    else
       echo "Processing ${SRR}"
       FastK -v -k"$K" -t1 -N"$OUT" -M"$MEM_GB" -T"$THREADS" \
          "$R1" "$R2" \
            > "${OUT}.out" 2> "${OUT}.err"
       [[ -s "${OUT}.ktab" ]] || {
        echo "FastK did not create a valid table for ${SRR}" >&2
        exit 1
       }

       touch "${OUT}.done"
    fi

    rm -f -- "/tmp/${SRR}_1.fastq" "/tmp/${SRR}_2.fastq"
    PARTS+=("$OUT")
done < "$SRR_LIST"

if [[ -e "${OUTDIR}/${SP}.Fastmerge.done" &&  -s "${OUTDIR}/${SP}.${K}.ktab" ]]; then
 Fastmerge -ht -T"$THREADS" \
    "${OUTDIR}/${SP}.${K}" \
    "${PARTS[@]}" \
    > "${OUTDIR}/${SP}.Fastmerge.out" \
    2> "${OUTDIR}/${SP}.Fastmerge.err"
fi
