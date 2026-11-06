# ---- BvAS project configuration (source this from every script) ----
#
# Everything lives on NRDStor (no purge policy; not backed up, so code and metadata are also kept in git).
# Set BULK_ON_NRDSTOR=0 to keep FASTQ and BAM files on $WORK instead (purged after 6 months untouched).
export BULK_ON_NRDSTOR=1

# Which analysis the scripts work on. Override on the command line, for example
#   ANALYSIS=full sbatch workflow/03_align.sbatch      (sbatch passes the variable on to the job)
#   pilot  33 runs, RefBeet-3.0 only, results in results/
#   full   465 runs, RefBeet-3.0 plus organelle genomes, results in results_full/
export ANALYSIS="${ANALYSIS:-pilot}"

if [ -z "${NRDSTOR:-}" ] || [ ! -d "${NRDSTOR}" ]; then
  echo "ERROR: \$NRDSTOR is not set or not mounted on this node." >&2
  return 1 2>/dev/null || exit 1
fi

export PROJ="${NRDSTOR}/bvas"
if [ "$BULK_ON_NRDSTOR" = "1" ]; then export BULK="${PROJ}/bulk"; else export BULK="${WORK}/bvas"; fi

# shared by all analyses
export ENV_DIR="${PROJ}/envs/bvas"
export CONDA_PKGS_DIRS="${PROJ}/conda_pkgs"
export REF_DIR="${PROJ}/ref"
export META_DIR="${PROJ}/metadata"
export LOG_DIR="${PROJ}/logs"
export SCRIPTS_DIR="${PROJ}/scripts"
export FASTQ_DIR="${BULK}/fastq"
export SCRATCH_DIR="${BULK}/tmp"
export GENES_GFF="${REF_DIR}/BeetSet-3.genes.2604.gff"
export GENES_GTF="${REF_DIR}/BeetSet-3.gtf"
export NUCLEAR_FA="${REF_DIR}/RefBeet-3.0.fa"
export ORGANELLE_PREFIX="organelle_"          # sequence name prefix of the added organelle genomes

case "$ANALYSIS" in
  pilot)
    export RES_DIR="${PROJ}/results"
    export QC_DIR="${PROJ}/qc"
    export BAM_DIR="${BULK}/bam"
    export GENOME_FA="$NUCLEAR_FA"
    export STAR_INDEX="${REF_DIR}/star_index"
    export SAMPLES="${META_DIR}/samples_pilot.tsv"
    export COMPARISONS="${META_DIR}/comparisons_pilot.tsv"
    export DELETE_FASTQ_AFTER_QUANT=0
    ;;
  full)
    export RES_DIR="${PROJ}/results_full"
    export QC_DIR="${PROJ}/qc_full"
    export BAM_DIR="${BULK}/bam_full"
    export GENOME_FA="${REF_DIR}/RefBeet-3.0_organelles.fa"     # written by 01c_add_organelles.sh
    export STAR_INDEX="${REF_DIR}/star_index_organelles"
    export SAMPLES="${META_DIR}/samples_all.tsv"
    export COMPARISONS="${META_DIR}/comparisons_all.tsv"
    # STAR and Salmon both read the FASTQ files, so they are deleted at the end of 06_salmon_quant
    export DELETE_FASTQ_AFTER_QUANT=1
    ;;
  *)
    echo "ERROR: unknown ANALYSIS '$ANALYSIS' (use pilot or full)" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac
# A comparison subset (same columns as the full table) can be given to run only new comparisons,
# for example COMPARISONS_FILE=metadata/comparisons_new.tsv sbatch --array=1-3 workflow/08_rmats.sbatch
if [ -n "${COMPARISONS_FILE:-}" ]; then export COMPARISONS="$COMPARISONS_FILE"; fi


export RUNRES_DIR="${RES_DIR}/per_run"        # SJ.out.tab, ReadsPerGene, Log.final.out, flagstat per run
export ST_DIR="${RES_DIR}/stringtie"
export MERGE_DIR="${RES_DIR}/merged"          # merged transcriptome, gffcompare, pooled junctions
export FILT_DIR="${RES_DIR}/filtered"         # structurally filtered transcriptome and its salmon index
export FINAL_GTF="${FILT_DIR}/final.gtf"      # after the expression filter; used by SUPPA2 and rMATS
export SALMON_DIR="${RES_DIR}/salmon"
export SUPPA_DIR="${RES_DIR}/suppa"

# Transcript merging and junction pooling (04b). Runs below MERGE_MIN_UNIQ percent unique mapping are
# left out of the merge (seed and germination libraries, fungus-dominated samples); a junction is pooled
# when it has JUNCTION_MIN_READS reads in at least JUNCTION_MIN_RUN_FRACTION of the merged runs.
export MERGE_MIN_UNIQ=60
export JUNCTION_MIN_READS=3
export JUNCTION_MIN_RUN_FRACTION=0.02
# Expression filter for non-reference transcripts before SUPPA2 and rMATS (filter_transcripts.py):
#   runs   TPM >= MIN_TPM in at least MIN_RUNS runs (used for the pilot set)
#   group  TPM >= MIN_TPM in every run of at least one comparison group
export MIN_TPM=1
export MIN_RUNS=2
case "$ANALYSIS" in
  pilot) export EXPR_FILTER=runs ;;
  full)  export EXPR_FILTER=group ;;
esac
# rMATS event universe: 0 = events defined by the final GTF only (same annotation as SUPPA2),
# 1 = also search unannotated splice sites (--novelSS).
export RMATS_NOVELSS=0
export RMATS_TAG="annot"
export RMATS_DIR="${RES_DIR}/rmats_${RMATS_TAG}"

mkdir -p "$PROJ" "$REF_DIR" "$META_DIR" "$LOG_DIR" "$QC_DIR" "$RUNRES_DIR" "$ST_DIR" \
         "$FASTQ_DIR" "$BAM_DIR" "$SCRATCH_DIR"

activate_env() {
  # module and conda shell functions reference unset variables, so relax "set -u" while they run
  local restore_u=0
  case $- in *u*) restore_u=1; set +u;; esac
  module purge
  module load anaconda
  conda activate "$ENV_DIR"
  [ "$restore_u" = 1 ] && set -u
  return 0
}

# Pick a temp directory for one job: node-local $TMPDIR if it has >= $1 GB free, else $SCRATCH_DIR.
job_tmp() { # usage: TMP=$(job_tmp 60 NAME)
  local need_gb="$1" name="$2" base=""
  if [ -n "${TMPDIR:-}" ] && [ -d "$TMPDIR" ]; then
    local free_gb; free_gb=$(df -Pk "$TMPDIR" | awk 'NR==2{print int($4/1024/1024)}')
    [ "${free_gb:-0}" -ge "$need_gb" ] && base="$TMPDIR"
  fi
  [ -z "$base" ] && base="$SCRATCH_DIR"
  local d="${base}/${name}.${SLURM_JOB_ID:-$$}"; mkdir -p "$d"; echo "$d"
}

# Read a column (by header name) of data line K of the sample sheet
sample_field() { # usage: sample_field LINE COLNAME
  local line="$1" col="$2"
  awk -F'\t' -v L="$line" -v C="$col" 'NR==1{for(i=1;i<=NF;i++) if($i==C) c=i} NR==L+1{print $c}' "$SAMPLES"
}

# Exit quietly when data line K of the sample sheet is marked use_das=excluded (see curation_note).
skip_if_excluded() { # usage: skip_if_excluded LINE
  if [ "$(sample_field "$1" use_das)" = "excluded" ]; then
    echo "run $(sample_field "$1" run) is marked excluded in $SAMPLES; nothing to do"
    exit 0
  fi
}
