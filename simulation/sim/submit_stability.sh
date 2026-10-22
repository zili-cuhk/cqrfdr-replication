#!/bin/bash
# Submit the split-stability study (A2) as an LSF job array on seadragon.
#   bash sim/submit_stability.sh                      # SS-MS with the CV screen (+ 2 MDS runs)
#   CQR_STAB_RULE=path bash sim/submit_stability.sh    # SS-MS with the path-ordered screen
# 3 scenarios x 20 data sets = 60 tasks; each runs 50 SS-MS splits + 2 MDS
# runs on one data set (about 1.5-3 h with 4 cores at p = 1000).
# Re-submitting is safe: finished data sets are skipped.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
mkdir -p logs
NTASK=60
bsub <<EOF
#BSUB -J cqr_stab${CQR_STAB_RULE:+_$CQR_STAB_RULE}[1-${NTASK}]
#BSUB -q medium
#BSUB -W 24:00
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -o logs/stab_%J_%I.out
#BSUB -e logs/stab_%J_%I.err
module load R/4.5.1
export CQR_SIM_ROOT=${ROOT} CQR_STAB_RULE=${CQR_STAB_RULE:-cv}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript ${ROOT}/sim/stability.R \$LSB_JOBINDEX
EOF
echo "submitted cqr_stab${CQR_STAB_RULE:+_$CQR_STAB_RULE}: ${NTASK} tasks"
