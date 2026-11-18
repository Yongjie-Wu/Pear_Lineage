#!/bin/bash

set -euo pipefail

# Usage
if [ $# -ne 1 ]; then
    echo "Usage: bash pear_origin_classifier_genome.sh <sample.genome.fasta>"
    echo "Optional: THREADS=10 ASIAN_WEIGHT=4.428461227 ASIAN_MAX_EURO_SHARE=0.15 EUROPEAN_MIN_EURO_SHARE=0.80"
    exit 1
fi

GENOME="$1"
if [ ! -f "$GENOME" ]; then
    echo "Error: file not found: $GENOME"
    exit 1
fi

for tool in makeblastdb blastn awk; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: required program not found in PATH: $tool" >&2
        exit 1
    fi
done

THREADS=${THREADS:-10}
ASIAN_WEIGHT=${ASIAN_WEIGHT:-4.428461227}
ASIAN_MAX_EURO_SHARE=${ASIAN_MAX_EURO_SHARE:-0.15}
EUROPEAN_MIN_EURO_SHARE=${EUROPEAN_MIN_EURO_SHARE:-0.80}

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: THREADS must be a positive integer" >&2
    exit 1
fi

if ! awk -v w="$ASIAN_WEIGHT" -v a="$ASIAN_MAX_EURO_SHARE" -v e="$EUROPEAN_MIN_EURO_SHARE" '
    BEGIN { exit !(w+0 > 0 && a+0 >= 0 && e+0 <= 1 && a+0 < e+0) }
'; then
    echo "Error: require ASIAN_WEIGHT > 0 and 0 <= ASIAN_MAX_EURO_SHARE < EUROPEAN_MIN_EURO_SHARE <= 1" >&2
    exit 1
fi

base=$(basename "$GENOME")
prefix=${base%.*}

DB_DIR="db"
HIT_DIR="blast_hit"
PERFECT_DIR="blast_perfect"

mkdir -p "$DB_DIR" "$HIT_DIR" "$PERFECT_DIR"

STAT="${prefix}_statistics.tsv"
IDENT="${prefix}_identification.txt"

# Marker sequences
QUERY_TMP=$(mktemp "${prefix}.XXXXXX.query.fa")
CHR_TMP=$(mktemp "${prefix}.XXXXXX.chr.list")
trap 'rm -f "$QUERY_TMP" "$CHR_TMP"' EXIT

cat > "$QUERY_TMP" <<'EOF'
>European_pear
TTGACTCAAAAACGTAAAATAACAAT
>Asian_pear
GGGTAATAAAATGAATCAAAACACTTAGAAACACAAACTAAGACACTTT
EOF

# Include every FASTA sequence; sequences without marker hits are Unknown.
grep '^>' "$GENOME" | sed 's/^>//; s/[[:space:]].*//; s/:.*//' | awk '!seen[$0]++' > "$CHR_TMP"

if [ ! -s "$CHR_TMP" ]; then
    echo "Error: no FASTA headers found in $GENOME"
    rm -f "$QUERY_TMP" "$CHR_TMP"
    exit 1
fi

# Output headers
echo -e "Genome\tAsian\tEuropean\tAsianWeight\tCorrectedEuropean\tCorrectedAsian\tEuropeanShare\tAsianShare\tType\tConfidence\tAsianChr\tEuropeanChr\tHybridChr\tUnknownChr" > "$STAT"
echo -e "Genome\tChr\tType\tEuropean\tAsian\tE/A\tAsianWeight\tCorrectedEuropean\tCorrectedAsian\tEuropeanShare\tAsianShare\tConfidence" > "$IDENT"

echo ">>> ${prefix}"

# Build BLAST db
if [ ! -f "${DB_DIR}/${prefix}_db.nsq" ] && [ ! -f "${DB_DIR}/${prefix}_db.ndb" ]; then
    makeblastdb -in "$GENOME" -dbtype nucl -out "${DB_DIR}/${prefix}_db" >/dev/null 2>&1
fi

# BLAST
blastn \
    -query "$QUERY_TMP" \
    -db "${DB_DIR}/${prefix}_db" \
    -out "${HIT_DIR}/${prefix}.tsv" \
    -evalue 1e-3 \
    -word_size 7 \
    -dust no \
    -soft_masking false \
    -task blastn-short \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp" \
    -num_threads "$THREADS"

# Perfect matches
awk '$3==100 && $13==100 && $4>=25' \
    "${HIT_DIR}/${prefix}.tsv" > "${PERFECT_DIR}/${prefix}.perfect.tsv"

# Genome-level statistics
read asian european < <(
    awk '
    {
        if($1=="Asian_pear"){a++}
        else if($1=="European_pear"){e++}
    }
    END{
        if(a=="") a=0
        if(e=="") e=0
        print a, e
    }' "${PERFECT_DIR}/${prefix}.perfect.tsv"
)

# Chromosome-level identification
# Asian markers are multiplied by ASIAN_WEIGHT because the validated Asian
# motif occurs at a lower density. Default 4.428461227 is the group-median
# marker-density ratio from 8 European and 29 Asian reference genomes:
# 2696.73499 / 537.40328 hits/Gb = 5.018084.
#   No pear-marker hits                         -> Unknown
#   corrected European share <= 0.15            -> Asian
#   corrected European share >= 0.80            -> European
#   otherwise                                   -> Hybrid
awk -v genome="$prefix" -v w="$ASIAN_WEIGHT" \
    -v alo="$ASIAN_MAX_EURO_SHARE" -v ehi="$EUROPEAN_MIN_EURO_SHARE" '
FNR==NR{
    chr=$1
    if(chr!="" && !seen_chr[chr]++){
        chr_order[++n]=chr
        A[chr]=0
        E[chr]=0
    }
    next
}
{
    chr=$2
    gsub(/:.*/,"",chr)

    if(!(chr in seen_chr)){
        seen_chr[chr]=1
        chr_order[++n]=chr
        A[chr]=0
        E[chr]=0
    }

    if($1=="Asian_pear"){A[chr]++}
    else if($1=="European_pear"){E[chr]++}
}
END{
    for(i=1; i<=n; i++){
        c = chr_order[i]
        ce = E[c]
        ca = w*A[c]
        total = ce + ca
        if(total==0){
            type="Unknown"
            es=0
            as=0
            conf=0
        } else {
            es = ce / total
            as = ca / total
            if(es <= alo){
                type="Asian"
            } else if(es >= ehi){
                type="European"
            } else {
                type="Hybrid"
            }
            conf=(es>=as ? es : as)
        }
        printf "%s\t%s\t%s\t%d\t%d\t%d/%d\t%.6g\t%.3f\t%.3f\t%.6f\t%.6f\t%.6f\n", \
            genome,c,type,E[c],A[c],E[c],A[c],w,ce,ca,es,as,conf
    }
}' "$CHR_TMP" "${PERFECT_DIR}/${prefix}.perfect.tsv" \
| sort -k2,2V >> "$IDENT"

# Genome-level call from chromosome composition. A genome is Hybrid whenever
# it contains both Asian and European chromosomes, or any internally mixed
# chromosome. Unknown scaffolds do not override informative chromosomes.
read asian_chr european_chr hybrid_chr unknown_chr < <(
    awk -F '\t' 'NR>1 {n[$3]++} END {
        print n["Asian"]+0, n["European"]+0, n["Hybrid"]+0, n["Unknown"]+0
    }' "$IDENT"
)

awk -v genome="$prefix" -v a="$asian" -v e="$european" -v w="$ASIAN_WEIGHT" \
    -v ac="$asian_chr" -v ec="$european_chr" -v hc="$hybrid_chr" -v uc="$unknown_chr" '
BEGIN {
    ce=e; ca=w*a; total=ce+ca
    if(total==0){ es=0; as=0; conf=0 }
    else { es=ce/total; as=ca/total; conf=(es>=as ? es : as) }
    if(hc>0 || (ac>0 && ec>0)) type="Hybrid"
    else if(ec>0) type="European"
    else if(ac>0) type="Asian"
    else type="Unknown"
    printf "%s\t%d\t%d\t%.6g\t%.3f\t%.3f\t%.6f\t%.6f\t%s\t%.6f\t%d\t%d\t%d\t%d\n", \
        genome,a,e,w,ce,ca,es,as,type,conf,ac,ec,hc,uc
}' >> "$STAT"

# Report chromosome calls.
awk 'NR>1{print $2"\t"$3}' "$IDENT"

echo "DONE"

# Cleanup
rm -f "$QUERY_TMP" "$CHR_TMP"
trap - EXIT
