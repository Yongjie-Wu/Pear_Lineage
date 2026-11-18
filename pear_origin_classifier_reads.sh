#!/bin/bash
set -euo pipefail

# Usage
if [ $# -lt 2 ]; then
    echo "Usage:"
    echo "  bash pear_origin_classifier_reads.sh sample_prefix read1.fq.gz [read2.fq.gz]"
    echo
    echo "Default: scan all reads."
    echo "Optional sampling: SAMPLE_READS=1000000 bash pear_origin_classifier_reads.sh ..."
    echo "Optional decompression threads: THREADS=10 bash pear_origin_classifier_reads.sh ..."
    echo "Optional classification parameters: MIN_TOTAL_HITS=10 BASE_ASIAN_WEIGHT=4.428461227 READ_DETECTABILITY_RATIO=auto"
    echo "  MAX_PURE_MINOR_SHARE=0.029411765 DETECTABILITY_SAMPLE_READS=5000"
    exit 1
fi

PREFIX=$1
R1=$2
R2=${3:-""}

# Parameters
THREADS=${THREADS:-10}
SAMPLE_READS=${SAMPLE_READS:-0}  # 0 = all reads; positive integer = random reads per FASTQ
MIN_TOTAL_HITS=${MIN_TOTAL_HITS:-10}
BASE_ASIAN_WEIGHT=${BASE_ASIAN_WEIGHT:-4.428461227}
READ_DETECTABILITY_RATIO=${READ_DETECTABILITY_RATIO:-auto}
DETECTABILITY_SAMPLE_READS=${DETECTABILITY_SAMPLE_READS:-5000}
MAX_PURE_MINOR_SHARE=${MAX_PURE_MINOR_SHARE:-0.029411765}  # 1/34 chromosome copies

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer" >&2
    exit 1
fi

if ! [[ "$SAMPLE_READS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: SAMPLE_READS must be a non-negative integer" >&2
    exit 1
fi

if ! [[ "$MIN_TOTAL_HITS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MIN_TOTAL_HITS must be a positive integer" >&2
    exit 1
fi

if ! [[ "$DETECTABILITY_SAMPLE_READS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: DETECTABILITY_SAMPLE_READS must be a positive integer" >&2
    exit 1
fi

if ! awk -v w="$BASE_ASIAN_WEIGHT" -v m="$MAX_PURE_MINOR_SHARE" '
    BEGIN { exit !(w+0 > 0 && m+0 > 0 && m+0 < 0.5) }
'; then
    echo "ERROR: require BASE_ASIAN_WEIGHT > 0 and 0 < MAX_PURE_MINOR_SHARE < 0.5" >&2
    exit 1
fi

if [ ! -r "$R1" ]; then
    echo "ERROR: cannot read R1: $R1" >&2
    exit 1
fi

if [ -n "$R2" ] && [ ! -r "$R2" ]; then
    echo "ERROR: cannot read R2: $R2" >&2
    exit 1
fi

# Read detectability
# A 50-bp exact motif has fewer valid start positions than a 26-bp motif and
# is more likely to contain a sequencing error. Estimate that technical bias
# from a small quality/length sample; marker counts below still use all reads.
if [ "$READ_DETECTABILITY_RATIO" = "auto" ]; then
    if ! command -v python3 >/dev/null 2>&1; then
        echo "ERROR: python3 is required for read detectability estimation" >&2
        exit 1
    fi
    READ_DETECTABILITY_RATIO=$(python3 - "$DETECTABILITY_SAMPLE_READS" "$R1" "$R2" <<'PY'
import gzip, math, sys
from collections import Counter

limit = int(sys.argv[1])
paths = [p for p in sys.argv[2:] if p]

def sample(path):
    lengths = Counter()
    psum, pn = [], []
    n = 0
    with gzip.open(path, "rt", errors="replace") as h:
        while n < limit:
            if not h.readline():
                break
            seq = h.readline().rstrip("\r\n")
            h.readline()
            qual = h.readline().rstrip("\r\n")
            if len(seq) != len(qual):
                continue
            lengths[len(seq)] += 1
            while len(psum) < len(qual):
                psum.append(0.0); pn.append(0)
            for i, char in enumerate(qual):
                q = max(0, ord(char) - 33)
                psum[i] += 1.0 - 10.0 ** (-q / 10.0)
                pn[i] += 1
            n += 1
    if n == 0:
        raise SystemExit("no readable FASTQ records for detectability estimate: " + path)
    return lengths, [s/c for s, c in zip(psum, pn)]

def effective_windows(lengths, cycle_p, motif_len):
    total = 0.0
    for read_len, count in lengths.items():
        if read_len < motif_len:
            continue
        logs = [math.log(max(cycle_p[i], 1e-12)) for i in range(read_len)]
        value = sum(logs[:motif_len])
        total += count * math.exp(value)
        for start in range(1, read_len - motif_len + 1):
            value += logs[start + motif_len - 1] - logs[start - 1]
            total += count * math.exp(value)
    return total

e26 = a50 = 0.0
for path in paths:
    lengths, cycle_p = sample(path)
    e26 += effective_windows(lengths, cycle_p, 26)
    a50 += effective_windows(lengths, cycle_p, 50)
if a50 <= 0:
    raise SystemExit("cannot estimate 50-bp marker detectability")
print(f"{e26/a50:.6f}")
PY
    )
fi

if ! awk -v r="$READ_DETECTABILITY_RATIO" 'BEGIN { exit !(r+0 > 0) }'; then
    echo "ERROR: READ_DETECTABILITY_RATIO must be 'auto' or a positive number" >&2
    exit 1
fi
ASIAN_WEIGHT=$(awk -v b="$BASE_ASIAN_WEIGHT" -v r="$READ_DETECTABILITY_RATIO" 'BEGIN { printf "%.6f", b*r }')

# Motifs (forward and reverse complement)
EURO_F="TTGACTCAAAAACGTAAAATAACAAT"
ASIA_F="GGGTAATAAAATGAATCAAAACACTTAGAAACACAAACTAAGACACTTT"
EURO_R="ATTGTTATTTTACGTTTTTGAGTCAA"
ASIA_R="AAAGTGTCTTAGTTTGTGTTTCTAAGTGTTTTGATTCATTTTATTACCC"

RESULT="${PREFIX}_reads_stat.txt"

echo ">>> Processing $PREFIX"
echo "  Base Asian weight:          $BASE_ASIAN_WEIGHT"
echo "  Read detectability ratio:   $READ_DETECTABILITY_RATIO"
echo "  Sample-specific Asian weight: $ASIAN_WEIGHT"
if [ "$SAMPLE_READS" -eq 0 ]; then
    echo "[1/4] Scanning all reads..."
else
    echo "[1/4] Randomly sampling up to $SAMPLE_READS reads from each FASTQ..."
fi

# Count motif-positive sequence lines in one FASTQ.
# The full-read path streams data and does not create a large .seq file.
count_fastq() {
    local fastq=$1

    if [ "$SAMPLE_READS" -eq 0 ]; then
        pigz -dc -p "$THREADS" "$fastq" \
        | LC_ALL=C awk -v ef="$EURO_F" -v er="$EURO_R" -v af="$ASIA_F" -v ar="$ASIA_R" '
            NR % 4 == 2 {
                if (index($0, ef) || index($0, er)) european++
                if (index($0, af) || index($0, ar)) asian++
            }
            END { print european + 0, asian + 0 }
        '
    else
        pigz -dc -p "$THREADS" "$fastq" \
        | awk 'NR % 4 == 2' \
        | shuf -n "$SAMPLE_READS" \
        | LC_ALL=C awk -v ef="$EURO_F" -v er="$EURO_R" -v af="$ASIA_F" -v ar="$ASIA_R" '
            {
                if (index($0, ef) || index($0, er)) european++
                if (index($0, af) || index($0, ar)) asian++
            }
            END { print european + 0, asian + 0 }
        '
    fi
}

echo "  -> R1: $R1"
R1_COUNTS=$(count_fastq "$R1")
read -r EURO_COUNT ASIA_COUNT <<< "$R1_COUNTS"

if [ -n "$R2" ]; then
    echo "  -> R2: $R2"
    R2_COUNTS=$(count_fastq "$R2")
    read -r EURO_R2 ASIA_R2 <<< "$R2_COUNTS"
    EURO_COUNT=$((EURO_COUNT + EURO_R2))
    ASIA_COUNT=$((ASIA_COUNT + ASIA_R2))
fi

echo "[2/4] Motif counting complete"
echo "  European hits: $EURO_COUNT"
echo "  Asian hits:    $ASIA_COUNT"

# Classification
# Asian markers are multiplied by a sample-specific ASIAN_WEIGHT:
# BASE_ASIAN_WEIGHT=4.428461227 corrects the genome marker-density ratio
# (2696.73499/537.40328 hits/Gb; 8 European and 29 Asian references).
# READ_DETECTABILITY_RATIO corrects motif-length and FASTQ-quality bias.
#   raw total hits < MIN_TOTAL_HITS: Unknown
#   corrected minor share < 1/34: corresponding pure origin
#   corrected minor share >= 1/34: Hybrid
echo "[3/4] Classification..."

TOTAL=$((EURO_COUNT + ASIA_COUNT))

if [ "$TOTAL" -lt "$MIN_TOTAL_HITS" ]; then
    TYPE="Unknown"
    CONF="0.000"
    CORRECTED_EUROPEAN="0.000"
    CORRECTED_ASIAN="0.000"
    EUROPEAN_SHARE="0.000000"
    ASIAN_SHARE="0.000000"
else
    read -r TYPE CONF CORRECTED_EUROPEAN CORRECTED_ASIAN EUROPEAN_SHARE ASIAN_SHARE < <(
        awk -v e="$EURO_COUNT" -v a="$ASIA_COUNT" -v w="$ASIAN_WEIGHT" \
            -v minor_cut="$MAX_PURE_MINOR_SHARE" '
        BEGIN {
            ce=e; ca=w*a; es=ce/(ce+ca); as=ca/(ce+ca)
            if(es<minor_cut) type="Asian"
            else if(as<minor_cut) type="European"
            else type="Hybrid"
            conf=(es>=as ? es : as)
            printf "%s %.6f %.3f %.3f %.6f %.6f\n", type,conf,ce,ca,es,as
        }'
    )
fi

# Output
echo "[4/4] Writing result..."

printf 'Sample\tType\tConfidence\tEuropean\tAsian\tTotalHit\tBaseAsianWeight\tReadDetectabilityRatio\tAsianWeight\tCorrectedEuropean\tCorrectedAsian\tEuropeanShare\tAsianShare\n' > "$RESULT"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$PREFIX" "$TYPE" "$CONF" "$EURO_COUNT" "$ASIA_COUNT" "$TOTAL" \
    "$BASE_ASIAN_WEIGHT" "$READ_DETECTABILITY_RATIO" "$ASIAN_WEIGHT" \
    "$CORRECTED_EUROPEAN" "$CORRECTED_ASIAN" "$EUROPEAN_SHARE" "$ASIAN_SHARE" >> "$RESULT"

cat "$RESULT"
echo "DONE"
