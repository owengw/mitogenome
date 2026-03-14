#!/usr/bin/env bash
# =============================================================================
# mito_pipeline.sh — Driver script for mitogenome SLURM array pipeline
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
REFERENCE=""
INPUT_DIR=""
OUTPUT_DIR=""
CONDA_ENV=""
SCRIPT_DIR="$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" && pwd)"
RAW_INPUT=false
ADAPTER="auto"
THREADS=8
MIN_MAPQ=20
MIN_DEPTH=3
MIN_AF=0.7
MIN_BQ=20
KEEP_INTERMEDIATES=false
DRY_RUN=false
FROM_STEP=1

# ── File naming / population ID ───────────────────────────────────────────────
FILE_PATTERN="*_f_paired.fastq.gz"  # glob to find forward read files
REV_PATTERN=""                      # glob for reverse reads (default: auto-derived from FILE_PATTERN)
LANE_PATTERN="_L00[0-9]*"           # sed-stripped suffix to derive sample ID
POP_REGEX='-([A-Z]+)\d'             # Python regex: group(1) = population code
NO_MERGE=false                      # set true for single-lane data (skips cat step)

# ── Taxon / analysis parameters ───────────────────────────────────────────────
SPECIES="Unknown species"
GENETIC_CODE=2
CLADE="Chordata"
CLOCK_RATE_MIN=0.005
CLOCK_RATE_MAX=0.02
MIN_POP_SIZE=5
ANNOT_ENV=""
BEAST_ENV=""
BEAST_BIN=""

# ── Time/memory limits per step ───────────────────────────────────────────────
TIME_TRIM="06:00:00";   MEM_TRIM="16G"
TIME_MAP="12:00:00";    MEM_MAP="32G"
TIME_DEDUP="04:00:00";  MEM_DEDUP="16G"
TIME_CONS="02:00:00";   MEM_CONS="8G"
TIME_QC="01:00:00";     MEM_QC="8G"
TIME_ALIGN="01:00:00";  MEM_ALIGN="16G"
TIME_TREE="04:00:00";   MEM_TREE="16G"

# ── Colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

# ── Help ──────────────────────────────────────────────────────────────────────
usage() {
cat <<EOF
${BOLD}mito_pipeline.sh${RESET} — Mitogenome SLURM array pipeline driver

${BOLD}USAGE${RESET}
  bash mito_pipeline.sh -r <ref.fa> -i <fastq_dir> -o <out_dir> -e <conda_env> [OPTIONS]

${BOLD}REQUIRED${RESET}
  -r  Reference mitogenome FASTA
  -i  Input FASTQ directory
  -o  Output directory (created if absent)
  -e  Full path to main conda environment (mito_pipeline)

${BOLD}TAXON / ANALYSIS PARAMETERS${RESET}
  --species        <str>    Species name for report (default: "Unknown species")
  --genetic-code   <N>      NCBI genetic code table (default: 2, vertebrate mito)
                            Common values: 2=vertebrate, 4=mold/protozoan,
                            5=invertebrate, 9=echinoderm, 13=ascidian mito
  --clade          <str>    MitoZ clade for annotation (default: Chordata)
                            Options: Chordata, Arthropoda, Mollusca, Echinodermata
  --clock-rate-min <float>  BEAST2 clock rate prior lower bound (default: 0.005)
  --clock-rate-max <float>  BEAST2 clock rate prior upper bound (default: 0.02)
                            Suggested ranges:
                              Squamates:  0.005 – 0.02
                              Birds:      0.01  – 0.05
                              Mammals:    0.01  – 0.04
                              Insects:    0.01  – 0.05
  --min-pop-size   <N>      Minimum samples per population for BSP and popgen
                            analyses. Populations below this threshold are
                            excluded from per-population analyses but always
                            included in the combined run (default: 5)
  --annot-env      <path>   Conda env with MitoZ + BioPython + MAFFT (required)
  --beast-env      <path>   Conda env with BEAST2 + PAML (required)
  --beast-bin      <path>   Full path to BEAST2 binary (default: <beast-env>/bin/beast)

${BOLD}PATHS${RESET}
  -S  Directory containing step scripts (default: same dir as this script)

${BOLD}INPUT${RESET}
  -T  Input FASTQs are raw/untrimmed — run fastp step first
  -a  Adapter type if trimming: auto | nextera | truseq | none (default: auto)

${BOLD}PIPELINE PARAMETERS${RESET}
  -t  Threads per job (default: 8)
  -q  Min mapping quality (default: 20)
  -d  Min depth for consensus base call (default: 3)
  -f  Min allele frequency for consensus (default: 0.7)
  -m  Min base quality for mpileup (default: 20)

${BOLD}SLURM RESOURCE OVERRIDES${RESET}
  --time-map   <HH:MM:SS>  Wall time for mapping jobs   (default: 12:00:00)
  --time-trim  <HH:MM:SS>  Wall time for trim jobs      (default: 06:00:00)
  --time-dedup <HH:MM:SS>  Wall time for dedup jobs     (default: 04:00:00)
  --time-cons  <HH:MM:SS>  Wall time for consensus jobs (default: 02:00:00)
  --time-align <HH:MM:SS>  Wall time for alignment      (default: 01:00:00)
  --time-tree  <HH:MM:SS>  Wall time for tree           (default: 04:00:00)
  --mem-map    <Ng>         Memory for mapping jobs      (default: 32G)
  --mem-trim   <Ng>         Memory for trim jobs         (default: 16G)
  --mem-dedup  <Ng>         Memory for dedup jobs        (default: 16G)
  --mem-cons   <Ng>         Memory for consensus jobs    (default: 8G)
  --mem-align  <Ng>         Memory for alignment         (default: 16G)
  --mem-tree   <Ng>         Memory for tree              (default: 16G)

${BOLD}OTHER${RESET}
  --from-step <N>   Start pipeline from step N, skipping earlier steps (default: 1)
  -k                Keep intermediate files
  -n                Dry run — print sbatch commands without submitting
  -h                Show this help

${BOLD}FILE NAMING / POPULATION ID${RESET}
  --file-pattern  <glob>   Glob matching forward read files within input dir
                           (default: "*_f_paired.fastq.gz")
  --rev-pattern   <glob>   Glob matching reverse read files (default: auto-derived
                           by replacing first _f_ with _r_ in --file-pattern)
  --lane-pattern  <str>    sed-compatible pattern stripped from filename to get
                           sample ID — should cover lane suffix + file extension
                           (default: "_L00[0-9]*_f_paired.fastq.gz")
  --pop-regex     <regex>  Python regex with one capture group extracting the
                           population code from a sequence/sample ID
                           (default: "-([A-Z]+)\\d")
                           Examples:
                             "-([A-Z]+)\\d"    matches "1-SI75_..." → SI
                             "_([A-Z]{2,4})_" matches "SampleSI_42" → SI
                             "^([^_]+)"        matches "POP1_ind42" → POP1
  --no-merge               Skip lane-merging in mapping step; use when each
                           sample has only a single pair of FASTQ files

${BOLD}EXAMPLES${RESET}
  # Lizards
  bash mito_pipeline.sh -r ref.fa -i fastq/ -o out/ -e /path/to/mito_pipeline \\
      --species "Plestiodon longirostris" \\
      --annot-env /path/to/mito_annot --beast-env /path/to/mito_beast \\
      --min-pop-size 5

  # Birds
  bash mito_pipeline.sh -r ref.fa -i fastq/ -o out/ -e /path/to/mito_pipeline \\
      --species "Parus major" --clade Chordata --genetic-code 2 \\
      --clock-rate-min 0.01 --clock-rate-max 0.05 \\
      --annot-env /path/to/mito_annot --beast-env /path/to/mito_beast \\
      --min-pop-size 10

  # Insects
  bash mito_pipeline.sh -r ref.fa -i fastq/ -o out/ -e /path/to/mito_pipeline \\
      --species "Drosophila melanogaster" --clade Arthropoda --genetic-code 5 \\
      --clock-rate-min 0.01 --clock-rate-max 0.05 \\
      --annot-env /path/to/mito_annot --beast-env /path/to/mito_beast \\
      --min-pop-size 5

EOF
exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --time-map)        TIME_MAP="$2";        shift 2 ;;
        --time-trim)       TIME_TRIM="$2";       shift 2 ;;
        --time-dedup)      TIME_DEDUP="$2";      shift 2 ;;
        --time-cons)       TIME_CONS="$2";       shift 2 ;;
        --time-align)      TIME_ALIGN="$2";      shift 2 ;;
        --time-tree)       TIME_TREE="$2";       shift 2 ;;
        --mem-map)         MEM_MAP="$2";         shift 2 ;;
        --mem-trim)        MEM_TRIM="$2";        shift 2 ;;
        --mem-dedup)       MEM_DEDUP="$2";       shift 2 ;;
        --mem-cons)        MEM_CONS="$2";        shift 2 ;;
        --mem-align)       MEM_ALIGN="$2";       shift 2 ;;
        --mem-tree)        MEM_TREE="$2";        shift 2 ;;
        --from-step)       FROM_STEP="$2";       shift 2 ;;
        --species)         SPECIES="$2";         shift 2 ;;
        --genetic-code)    GENETIC_CODE="$2";    shift 2 ;;
        --clade)           CLADE="$2";           shift 2 ;;
        --clock-rate-min)  CLOCK_RATE_MIN="$2";  shift 2 ;;
        --clock-rate-max)  CLOCK_RATE_MAX="$2";  shift 2 ;;
        --min-pop-size)    MIN_POP_SIZE="$2";    shift 2 ;;
        --annot-env)       ANNOT_ENV="$2";       shift 2 ;;
        --beast-env)       BEAST_ENV="$2";       shift 2 ;;
        --beast-bin)       BEAST_BIN="$2";       shift 2 ;;
        --file-pattern)    FILE_PATTERN="$2";    shift 2 ;;
        --rev-pattern)     REV_PATTERN="$2";     shift 2 ;;
        --lane-pattern)    LANE_PATTERN="$2";    shift 2 ;;
        --pop-regex)       POP_REGEX="$2";       shift 2 ;;
        --no-merge)        NO_MERGE=true;        shift ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
set -- "${ARGS[@]}"

while getopts ":r:i:o:e:S:a:t:q:d:f:m:Tknh" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
        S) SCRIPT_DIR="$OPTARG" ;;
        T) RAW_INPUT=true ;;
        a) ADAPTER="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        q) MIN_MAPQ="$OPTARG" ;;
        d) MIN_DEPTH="$OPTARG" ;;
        f) MIN_AF="$OPTARG" ;;
        m) MIN_BQ="$OPTARG" ;;
        k) KEEP_INTERMEDIATES=true ;;
        n) DRY_RUN=true ;;
        h) usage ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
[[ -z "$REFERENCE" ]]   && err "Reference FASTA required (-r)"
[[ -z "$INPUT_DIR" ]]   && err "Input directory required (-i)"
[[ -z "$OUTPUT_DIR" ]]  && err "Output directory required (-o)"
[[ -z "$CONDA_ENV" ]]   && err "Conda environment path required (-e)"
[[ ! -f "$REFERENCE" ]] && err "Reference not found: $REFERENCE"
[[ ! -d "$INPUT_DIR" ]] && err "Input directory not found: $INPUT_DIR"
[[ ! -d "$CONDA_ENV" ]] && err "Conda environment not found: $CONDA_ENV"

[[ -z "$ANNOT_ENV" ]]   && err "Annotation conda environment required (--annot-env)"
[[ -z "$BEAST_ENV" ]]   && err "BEAST2 conda environment required (--beast-env)"
[[ ! -d "$ANNOT_ENV" ]] && err "Annotation conda environment not found: $ANNOT_ENV"
[[ ! -d "$BEAST_ENV" ]] && err "BEAST2 conda environment not found: $BEAST_ENV"

if [[ -z "$BEAST_BIN" ]]; then
    BEAST_BIN="${BEAST_ENV}/bin/beast"
fi
[[ ! -f "$BEAST_BIN" ]] && err "BEAST2 binary not found: $BEAST_BIN (use --beast-bin to specify)"

REQUIRED_SCRIPTS=("m01_map.sh" "m02_dedup_filter.sh" "m03_call_consensus.sh"
    "m04_qc_summary.sh" "m05_align.sh" "m06_tree.sh"
    "m07_popgen.sh" "m08_annotate.sh" "m09_dnds.sh" "m10_bsp.sh" "m11_report.sh")
[[ "$RAW_INPUT" == true ]] && REQUIRED_SCRIPTS=("m00_trim.sh" "${REQUIRED_SCRIPTS[@]}")
for s in "${REQUIRED_SCRIPTS[@]}"; do
    [[ ! -f "${SCRIPT_DIR}/${s}" ]] && err "Step script not found: ${SCRIPT_DIR}/${s}"
done

# ── Discover samples ──────────────────────────────────────────────────────────
# Auto-derive reverse pattern if not explicitly set
if [[ -z "$REV_PATTERN" ]]; then
    REV_PATTERN="${FILE_PATTERN/_f_/_r_}"
fi
mapfile -t ALL_F_FILES < <(find "$INPUT_DIR" -maxdepth 1 -name "${FILE_PATTERN}" | sort)
[[ ${#ALL_F_FILES[@]} -eq 0 ]] && err "No files matching '${FILE_PATTERN}' found in ${INPUT_DIR}"

# Strip lane suffix + extension to get unique sample IDs
# LANE_PATTERN is a sed BRE; escape dots for the extension part
mapfile -t SAMPLES < <(
    for f in "${ALL_F_FILES[@]}"; do
      if [[ "$NO_MERGE" == true ]]; then
        basename "$f" | sed 's/\.fastq\.gz$//' | sed 's/\.fq\.gz$//'
      else
        basename "$f" | sed "s/${LANE_PATTERN}//" | sed 's/\.fastq\.gz$//' | sed 's/\.fq\.gz$//'
      fi
    done | sort -u
)

N_SAMPLES=${#SAMPLES[@]}
ARRAY_IDX="1-${N_SAMPLES}"

# Shell-escape species name and pop regex so they survive word-splitting in sbatch EXTRA_ARGS
SPECIES_ESC="$(printf '%q' "${SPECIES}")"
POP_REGEX_ESC="$(printf '%q' "${POP_REGEX}")"

log "Samples discovered: ${BOLD}${N_SAMPLES}${RESET}"
log "Sample → lane mapping:"
for SAMPLE in "${SAMPLES[@]}"; do
    LANES=$(find "$INPUT_DIR" -maxdepth 1 -name "${SAMPLE}${LANE_PATTERN}" \
        | xargs -I{} basename {} | sort | tr '\n' ' ')
    if [[ "$NO_MERGE" == true ]]; then
        echo -e "  ${SAMPLE}  →  ${CYAN}${LANES}${RESET}  ${YELLOW}(no-merge)${RESET}"
    else
        echo -e "  ${SAMPLE}  →  ${CYAN}${LANES}${RESET}"
    fi
done
echo ""

# ── Print analysis parameters ─────────────────────────────────────────────────
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  ANALYSIS PARAMETERS${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "  Species:         ${SPECIES}"
echo -e "  Genetic code:    ${GENETIC_CODE}"
echo -e "  MitoZ clade:     ${CLADE}"
echo -e "  Clock rate:      U(${CLOCK_RATE_MIN}, ${CLOCK_RATE_MAX}) subs/site/Myr"
echo -e "  Min pop size:    ${MIN_POP_SIZE} (popgen + BSP per-population)"
echo -e "  Annot env:       ${ANNOT_ENV}"
echo -e "  BEAST2 env:      ${BEAST_ENV}"
echo -e "  BEAST2 binary:   ${BEAST_BIN}"
echo -e "  File pattern:    ${FILE_PATTERN}"
echo -e "  Rev pattern:     ${REV_PATTERN}"
echo -e "  Lane pattern:    ${LANE_PATTERN}"
echo -e "  Pop regex:       ${POP_REGEX}"
echo -e "  Lane merging:    $( [[ "$NO_MERGE" == true ]] && echo 'disabled (--no-merge)' || echo 'enabled' )"
echo ""

# ── Create output dirs ────────────────────────────────────────────────────────
mkdir -p "${OUTPUT_DIR}"/{trimmed,bam,vcf,consensus,logs,qc,sample_lists,aligned,tree}

# ── Write sample list ─────────────────────────────────────────────────────────
SAMPLE_LIST="${OUTPUT_DIR}/sample_lists/samples.txt"
printf '%s\n' "${SAMPLES[@]}" > "$SAMPLE_LIST"
log "Sample list written → ${SAMPLE_LIST}"

# ── Activate conda for local operations ───────────────────────────────────────
set +u
source ~/.bash_profile
conda activate "$CONDA_ENV"
set -u

# ── Index reference if needed ─────────────────────────────────────────────────
if [[ "$FROM_STEP" -le 1 ]]; then
    if [[ ! -f "${REFERENCE}.bwt" || ! -f "${REFERENCE}.fai" ]]; then
        log "Indexing reference..."
        [[ ! -f "${REFERENCE}.bwt" ]] && bwa index "$REFERENCE"
        [[ ! -f "${REFERENCE}.fai" ]] && samtools faidx "$REFERENCE"
        ok "Reference indexed"
    fi
fi

# ── Common args passed to array step scripts ──────────────────────────────────
COMMON_ARGS="-r ${REFERENCE} -i ${INPUT_DIR} -o ${OUTPUT_DIR} -s ${SAMPLE_LIST} -t ${THREADS} -e ${CONDA_ENV}"
KEEP_FLAG=""; [[ "$KEEP_INTERMEDIATES" == true ]] && KEEP_FLAG="-k"
MERGE_FLAG=""; [[ "$NO_MERGE" == true ]] && MERGE_FLAG="-M"

# ── sbatch helper ─────────────────────────────────────────────────────────────
submit_job() {
    local STEP_NAME="$1"
    local SCRIPT="$2"
    local TIME="$3"
    local MEM="$4"
    local IS_ARRAY="$5"
    local DEPEND_ARG="$6"
    local EXTRA_ARGS="$7"

    local ARRAY_FLAG=""
    [[ "$IS_ARRAY" == true ]] && ARRAY_FLAG="--array=${ARRAY_IDX}%20"
    [[ ! -f "$SCRIPT" ]] && err "Script not found: $SCRIPT"

    local JOB_ID
    if [[ "$DRY_RUN" == true ]]; then
        echo -e "${YELLOW}[DRY-RUN]${RESET} sbatch ${SCRIPT} ..."
        JOB_ID="DRY_RUN_${STEP_NAME}"
    else
        JOB_ID=$(sbatch \
            --job-name=mito_${STEP_NAME} \
            --output=${OUTPUT_DIR}/logs/${STEP_NAME}_%A_%a.out \
            --error=${OUTPUT_DIR}/logs/${STEP_NAME}_%A_%a.err \
            --time=${TIME} \
            --mem=${MEM} \
            --cpus-per-task=${THREADS} \
            --nodes=1 \
            --ntasks-per-node=1 \
            ${ARRAY_FLAG} \
            ${DEPEND_ARG} \
            --wrap="bash ${SCRIPT} ${EXTRA_ARGS}" \
            | awk '{print $NF}')
    fi

    echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  Submitted ${STEP_NAME} → Job ID: ${JOB_ID}${RESET}" >&2
    echo "$JOB_ID"
}

# =============================================================================
# SUBMIT JOBS
# =============================================================================
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  SUBMITTING JOBS (from step ${FROM_STEP})${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

DEPEND=""

# ── Step 00: Trim (optional) ──────────────────────────────────────────────────
if [[ "$RAW_INPUT" == true && "$FROM_STEP" -le 0 ]]; then
    log "Submitting Step 00: Trim (fastp)"
    JOB_TRIM=$(submit_job "m00_trim" "${SCRIPT_DIR}/m00_trim.sh" \
        "$TIME_TRIM" "$MEM_TRIM" true "" \
        "${COMMON_ARGS} -a ${ADAPTER} ${KEEP_FLAG}")
    DEPEND="--dependency=afterok:${JOB_TRIM}"
    COMMON_ARGS="-r ${REFERENCE} -i ${OUTPUT_DIR}/trimmed -o ${OUTPUT_DIR} -s ${SAMPLE_LIST} -t ${THREADS} -e ${CONDA_ENV}"
fi

# ── Step 01: Map ──────────────────────────────────────────────────────────────
if [[ "$FROM_STEP" -le 1 ]]; then
    log "Submitting Step 01: Map (bwa mem)"
    JOB_MAP=$(submit_job "m01_map" "${SCRIPT_DIR}/m01_map.sh" \
        "$TIME_MAP" "$MEM_MAP" true "$DEPEND" \
        "${COMMON_ARGS} ${KEEP_FLAG} ${MERGE_FLAG} -F ${FILE_PATTERN} -V ${REV_PATTERN}")
    DEPEND="--dependency=afterok:${JOB_MAP}"
fi

# ── Step 02: Dedup + filter ───────────────────────────────────────────────────
if [[ "$FROM_STEP" -le 2 ]]; then
    log "Submitting Step 02: Dedup + filter"
    JOB_DEDUP=$(submit_job "m02_dedup" "${SCRIPT_DIR}/m02_dedup_filter.sh" \
        "$TIME_DEDUP" "$MEM_DEDUP" true "$DEPEND" \
        "${COMMON_ARGS} -q ${MIN_MAPQ} ${KEEP_FLAG}")
    DEPEND="--dependency=afterok:${JOB_DEDUP}"
fi

# ── Step 03: Variant call + consensus ─────────────────────────────────────────
if [[ "$FROM_STEP" -le 3 ]]; then
    log "Submitting Step 03: Call + consensus"
    JOB_CONS=$(submit_job "m03_consensus" "${SCRIPT_DIR}/m03_call_consensus.sh" \
        "$TIME_CONS" "$MEM_CONS" true "$DEPEND" \
        "${COMMON_ARGS} -d ${MIN_DEPTH} -f ${MIN_AF} -m ${MIN_BQ} ${KEEP_FLAG}")
    DEPEND="--dependency=afterok:${JOB_CONS}"
fi

# ── Step 04: QC summary ───────────────────────────────────────────────────────
if [[ "$FROM_STEP" -le 4 ]]; then
    log "Submitting Step 04: QC summary"
    JOB_QC=$(submit_job "m04_qc" "${SCRIPT_DIR}/m04_qc_summary.sh" \
        "$TIME_QC" "$MEM_QC" false "$DEPEND" \
        "${COMMON_ARGS}")
    DEPEND="--dependency=afterany:${JOB_QC}"
fi

# ── Step 05: Align ────────────────────────────────────────────────────────────
if [[ "$FROM_STEP" -le 5 ]]; then
    log "Submitting Step 05: MAFFT alignment + consensus"
    JOB_ALIGN=$(submit_job "m05_align" "${SCRIPT_DIR}/m05_align.sh" \
        "$TIME_ALIGN" "$MEM_ALIGN" false "$DEPEND" \
        "-o ${OUTPUT_DIR} -e ${CONDA_ENV} -t ${THREADS}")
    DEPEND="--dependency=afterok:${JOB_ALIGN}"
fi

# ── Step 06: Tree ─────────────────────────────────────────────────────────────
if [[ "$FROM_STEP" -le 6 ]]; then
    log "Submitting Step 06: IQ-TREE phylogeny"
    JOB_TREE=$(submit_job "m06_tree" "${SCRIPT_DIR}/m06_tree.sh" \
        "$TIME_TREE" "$MEM_TREE" false "$DEPEND" \
        "-o ${OUTPUT_DIR} -e ${CONDA_ENV} -t ${THREADS} -n ${SPECIES_ESC}")
    DEPEND="--dependency=afterok:${JOB_TREE}"
fi

# ── Step 07: Population genetics ─────────────────────────────────────────────
if [[ "$FROM_STEP" -le 7 ]]; then
    log "Submitting Step 07: Population genetics (π, Tajima's D, Hd, Fu's Fs)"
    JOB_POPGEN=$(submit_job "m07_popgen" "${SCRIPT_DIR}/m07_popgen.sh" \
        "02:00:00" "16G" false "$DEPEND" \
        "-o ${OUTPUT_DIR} -e ${CONDA_ENV} -p ${MIN_POP_SIZE} -R ${POP_REGEX_ESC}")
    DEPEND="--dependency=afterok:${JOB_POPGEN}"
fi

# ── Step 08: Gene extraction + annotation ────────────────────────────────────
if [[ "$FROM_STEP" -le 8 ]]; then
    log "Submitting Step 08: Gene extraction + annotation"
    JOB_ANNOT=$(submit_job "m08_annot" "${SCRIPT_DIR}/m08_annotate.sh" \
        "02:00:00" "16G" false "$DEPEND" \
        "-o ${OUTPUT_DIR} -e ${ANNOT_ENV} -C ${GENETIC_CODE} -c ${CLADE}")
    DEPEND="--dependency=afterok:${JOB_ANNOT}"
fi

# ── Step 09: dN/dS analysis ───────────────────────────────────────────────────
if [[ "$FROM_STEP" -le 9 ]]; then
    log "Submitting Step 09: dN/dS (yn00)"
    JOB_DNDS=$(submit_job "m09_dnds" "${SCRIPT_DIR}/m09_dnds.sh" \
        "04:00:00" "16G" false "$DEPEND" \
        "-o ${OUTPUT_DIR} -e ${BEAST_ENV} -C ${GENETIC_CODE} -n ${SPECIES_ESC}")
    DEPEND="--dependency=afterok:${JOB_DNDS}"
fi

# ── Step 10: Bayesian Skyline Plot ────────────────────────────────────────────
if [[ "$FROM_STEP" -le 10 ]]; then
    log "Submitting Step 10: BSP (BEAST2)"
    JOB_BSP=$(submit_job "m10_bsp" "${SCRIPT_DIR}/m10_bsp.sh" \
        "48:00:00" "32G" false "$DEPEND" \
        "-o ${OUTPUT_DIR} -e ${BEAST_ENV} -b ${BEAST_BIN} -l ${CLOCK_RATE_MIN} -u ${CLOCK_RATE_MAX} -p ${MIN_POP_SIZE} -n ${SPECIES_ESC} -R ${POP_REGEX_ESC}")
    DEPEND="--dependency=afterok:${JOB_BSP}"
fi

# ── Step 11: Report ───────────────────────────────────────────────────────────
if [[ "$FROM_STEP" -le 11 ]]; then
    log "Submitting Step 11: HTML report"
    JOB_REPORT=$(submit_job "m11_report" "${SCRIPT_DIR}/m11_report.sh" \
        "00:30:00" "8G" false "$DEPEND" \
        "-o ${OUTPUT_DIR} -e ${CONDA_ENV} -s ${SPECIES_ESC} -p ${MIN_POP_SIZE} -S ${SCRIPT_DIR}")
    DEPEND="--dependency=afterok:${JOB_REPORT}"
fi