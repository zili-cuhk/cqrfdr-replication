#!/bin/bash
# Real-data analysis on seadragon (LSF).
#   bash submit_realdata.sh            # all tasks in tasks.csv
#   TASKS=tasks_targeted.csv bash submit_realdata.sh   # the targeted METABRIC analysis
#   ONLY="1-8" bash submit_realdata.sh   # a sub-range (main/mds tasks come first)
#   MEM=48 ONLY="1-4" bash submit_realdata.sh   # more memory (GB) for some tasks
#   CORES=10 ONLY="<mds_dai ids>" bash submit_realdata.sh # more cores (MDS-Dai: 10-fold CV in every split)
# Steps before this (see README.md): 00_install.R, 01_download.R (login node), 02_prepare.R.
# Re-submitting is safe: finished tasks are skipped.
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(pwd)
mkdir -p logs
TASKS=${TASKS:-tasks.csv}
TABLES=${TABLES:-tables}
NTASK=$(( $(wc -l < "$TASKS") - 1 ))
RANGE=${ONLY:-1-$NTASK}
MEM=${MEM:-32}   # GB; hard limit (-M) and reservation
CORES=${CORES:-4}  # cores per task (MDS-Dai: CORES=10, one core per CV fold)
bsub <<EOT
#BSUB -J cqr_real[${RANGE}]%600
#BSUB -q medium
#BSUB -W 24:00
#BSUB -n ${CORES}
#BSUB -R "span[hosts=1] rusage[mem=${MEM}GB]"
#BSUB -M ${MEM}GB
#BSUB -o logs/real_%J_%I.out
#BSUB -e logs/real_%J_%I.err
module load R/4.5.1
export CQR_RD_ROOT=${ROOT} CQR_RD_CORES=${CORES} CQR_RD_TASKS=${TASKS} CQR_RD_TABLES=${TABLES}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript ${ROOT}/R/03_task.R \$LSB_JOBINDEX
EOT
echo "submitted cqr_real: tasks ${RANGE} of ${NTASK}"
