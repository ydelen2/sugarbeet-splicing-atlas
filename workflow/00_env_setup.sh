#!/bin/bash
# Run ONCE on a Swan login node:  bash workflow/00_env_setup.sh
# The environment and package cache live on $NRDSTOR (HCC recommends $NRDSTOR for conda envs; $HOME is 20 GiB).
# Package list: environment/bvas.yml. After installation the exact resolved versions are exported to
# environment/bvas.exported.yml so the analysis can be reproduced.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/config.sh"
SPEC="$HERE/../environment/bvas.yml"

# Keep a full log of every run so a failed installation can be diagnosed afterwards
LOGFILE="$LOG_DIR/env_setup.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1
echo "log: $LOGFILE"
echo "host: $(hostname)  start: $(date)"
trap 'echo "FAILED at line $LINENO (exit $?) $(date)"' ERR

set +u   # module and conda shell functions reference unset variables
module purge
module load anaconda

# The user's ~/.condarc is left untouched: the package cache is set through the CONDA_PKGS_DIRS
# environment variable (exported in config.sh) and the channels come from environment/bvas.yml.
mkdir -p "$CONDA_PKGS_DIRS" "$(dirname "$ENV_DIR")"

CREATE="conda"; command -v mamba >/dev/null 2>&1 && CREATE="mamba"

if [ -x "$ENV_DIR/bin/STAR" ] && [ -x "$ENV_DIR/bin/stringtie" ] && [ -x "$ENV_DIR/bin/salmon" ]; then
  echo "Env already exists at $ENV_DIR"
else
  if [ -d "$ENV_DIR" ]; then
    echo "Incomplete environment found at $ENV_DIR, removing and reinstalling"
    rm -rf "$ENV_DIR"
  fi
  if [ -f "$SPEC" ]; then
    $CREATE env create -p "$ENV_DIR" -f "$SPEC"
  else
    echo "Spec file $SPEC not found, using the built-in package list"
    $CREATE create -y -p "$ENV_DIR" -c conda-forge -c bioconda \
      python=3.11 \
      star=2.7.11b hisat2=2.2.1 samtools=1.21 \
      fastp=0.24.0 fastqc multiqc \
      stringtie=2.2.3 gffread gffcompare \
      subread salmon=1.10.3 \
      suppa=2.4 "statsmodels<0.15" rmats=4.3.0 \
      aria2 pigz seqkit bedtools \
      pandas numpy
  fi
fi

conda activate "$ENV_DIR"
set -u
echo "== tool check =="
echo -n "STAR ";      STAR --version
fastp --version 2>&1 | awk 'NR==1'
samtools --version | awk 'NR==1'
echo -n "stringtie "; stringtie --version
echo -n "salmon ";    salmon --version 2>&1 | awk 'NR==1'
rmats.py --version 2>/dev/null || echo "rmats: check manually"
python -c "import statsmodels; print('statsmodels', statsmodels.__version__)"
if suppa.py generateEvents --help >/dev/null 2>&1; then echo "SUPPA2 ok"; else echo "SUPPA2 FAILED:"; suppa.py generateEvents --help 2>&1 | awk '{l[NR]=$0} END{for(i=NR-2;i<=NR;i++) if(i>0) print l[i]}'; fi

mkdir -p "$HERE/../environment"
# The raw export contains the cluster path and HCC channel URLs with an access token; publish only the
# sanitized version (plain name, public channels).
RAW_EXPORT="$ENV_DIR/conda-meta/bvas.export.raw.yml"
conda env export -p "$ENV_DIR" --no-builds > "$RAW_EXPORT"
python "$HERE/../scripts/analysis/sanitize_env_export.py" "$RAW_EXPORT" "$HERE/../environment/bvas.exported.yml"
echo "Exported package versions: environment/bvas.exported.yml"
echo "Environment ready: $ENV_DIR  end: $(date)"
