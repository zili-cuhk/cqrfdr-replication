#!/bin/bash
# Plasmode simulation on the real METABRIC design: 4 settings x 50 replicates = 200 tasks.
#   bash plasmode/submit_plasmode.sh            # all tasks
#   ONLY="1-4" bash plasmode/submit_plasmode.sh # a few first, to check timing
#   PM_DATA=BRCA_METABRIC_TARGETED PM_S=5,15 PM_B=0.15,0.3 bash plasmode/submit_plasmode.sh
#   PM_METHODS=msms bash plasmode/submit_plasmode.sh   # add MS-MS to finished tasks (same data)
# Re-submitting is safe: finished tasks are skipped.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
mkdir -p logs
PM_DATA=${PM_DATA:-BRCA_METABRIC}; PM_S=${PM_S:-10,30}; PM_B=${PM_B:-0.1,0.2}
PM_METHODS=${PM_METHODS:-ssms,msms,mds_dai,mds,hdcqr}
RANGE=${ONLY:-1-200}
MEM=${MEM:-32}
bsub <<EOT
#BSUB -J cqr_plasmode[${RANGE}]%400
#BSUB -q medium
#BSUB -W 24:00
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=${MEM}GB]"
#BSUB -M ${MEM}GB
#BSUB -o logs/plasmode_%J_%I.out
#BSUB -e logs/plasmode_%J_%I.err
module load R/4.5.1
export CQR_RD_ROOT=${ROOT} CQR_RD_CORES=4 CQR_PM_DATA=${PM_DATA} CQR_PM_S=${PM_S} CQR_PM_B=${PM_B} CQR_PM_METHODS=${PM_METHODS}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript ${ROOT}/plasmode/plasmode.R \$LSB_JOBINDEX
EOT
echo "submitted cqr_plasmode (${PM_DATA}; s=${PM_S}; b=${PM_B}; ${PM_METHODS}): tasks ${RANGE}"
