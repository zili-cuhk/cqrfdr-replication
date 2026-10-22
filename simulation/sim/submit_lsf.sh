#!/bin/bash
# Submit the simulation as an LSF job array on seadragon.
#   CQR_SIM_SET=df_pilot bash sim/submit_lsf.sh   # step 1: degrees-of-freedom pilot
#   bash sim/submit_lsf.sh                        # step 2: main simulation (CQR_SIM_SET=main)
#   ONLY="1-2401:100" bash sim/submit_lsf.sh      # a sub-range / stride of task ids
#   REPS_PER_TASK=2 bash sim/submit_lsf.sh        # several replications per task
#   CQR_SIM_ONLY=ds bash sim/submit_lsf.sh        # add one method to existing results
# Task id t runs replication ((t-1) %% chunks)*k + 1.. of scenario ceil(t/chunks).
# Re-submitting is safe: finished replications are skipped.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
module load R/4.5.1 2>/dev/null || true
export CQR_SIM_SET=${CQR_SIM_SET:-main}

REPS_PER_TASK=${REPS_PER_TASK:-1}
NSCEN=$(Rscript -e 'source("sim/config.R"); cat(nrow(make_scenarios()))')
REPS=$(Rscript -e 'source("sim/config.R"); cat(SIM_CONFIG$reps)')
CHUNKS=$(( (REPS + REPS_PER_TASK - 1) / REPS_PER_TASK ))
NTASK=$(( NSCEN * CHUNKS ))
RANGE=${ONLY:-1-$NTASK}

mkdir -p logs
# censoring calibration once per scenario set (a few minutes; runs on a compute node)
if [ ! -f "sim/calibration_${CQR_SIM_SET}.csv" ]; then
  echo "calibrating ${CQR_SIM_SET} ..."
  bsub -q medium -W 4:00 -n 1 -K -o "logs/calibrate_${CQR_SIM_SET}.out" \
       "module load R/4.5.1; cd ${ROOT}; CQR_SIM_SET=${CQR_SIM_SET} Rscript sim/calibrate.R" || exit 1
  tail -40 "logs/calibrate_${CQR_SIM_SET}.out"
fi

bsub <<EOF
#BSUB -J cqr_${CQR_SIM_SET}[${RANGE}]%800
#BSUB -q medium
#BSUB -W 24:00
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -o logs/${CQR_SIM_SET}_%J_%I.out
#BSUB -e logs/${CQR_SIM_SET}_%J_%I.err
module load R/4.5.1
export CQR_SIM_ROOT=${ROOT} CQR_SIM_SET=${CQR_SIM_SET} CQR_SIM_ONLY=${CQR_SIM_ONLY:-}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript ${ROOT}/sim/run.R \$LSB_JOBINDEX ${REPS_PER_TASK}
EOF
echo "submitted ${CQR_SIM_SET}: tasks ${RANGE} of ${NTASK} (${NSCEN} scenarios x ${CHUNKS} chunks)"
