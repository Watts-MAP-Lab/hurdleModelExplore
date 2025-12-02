#!/usr/bin/env bash
# Run the hurdle simulation across a range of seeds using up to 4 cores.
# Usage:
#   ./run_sim_4cores.sh <sim_R_script> <start_seed> <end_seed> [cores]
# Example:
#   ./run_sim_4cores.sh ./scripts/rCode/final_simulation.R 1 128 4

set -euo pipefail

# --------- Inputs ----------
SIM_SCRIPT="${1}"
START_SEED="${2:-1}"
END_SEED="${3:-10}"
CORES="${4:-4}"

# --------- Environment hygiene (avoid thread oversubscription) ----------
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# --------- Paths ----------
PROJECT_ROOT="$(pwd)"  # assume you run from project root
OUT_DIR="${PROJECT_ROOT}/data/hurdleCollapse"
LOG_DIR="${PROJECT_ROOT}/logs/hurdleCollapse"
mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# --------- Checks ----------
if ! command -v Rscript >/dev/null 2>&1; then
  echo "[ERROR] Rscript not found in PATH"; exit 1
fi
if [ ! -f "${SIM_SCRIPT}" ]; then
  echo "[ERROR] Simulation script not found: ${SIM_SCRIPT}"; exit 1
fi

echo "[INFO] Starting simulation runs"
echo "[INFO] Script: ${SIM_SCRIPT}"
echo "[INFO] Seeds:  ${START_SEED}..${END_SEED}"
echo "[INFO] Cores:  ${CORES}"
start_epoch=$(date +%s)

# --------- Runner using xargs -P (portable parallelism) ----------
# Each seed writes to logs/hurdleCollapse/seed_<seed>.log
seq "${START_SEED}" "${END_SEED}" | xargs -n 1 -P "${CORES}" -I{} sh -c '
  seed="$1"
  out_file="'"${OUT_DIR}"'/seedVal_${seed}_take3.RDS"
  log_file="'"${LOG_DIR}"'/seed_${seed}.log"

  # Skip if output already exists
  if [ -f "${out_file}" ]; then
    echo "$(date) [SKIP] seed ${seed}: ${out_file} exists" | tee -a "${log_file}"
    exit 0
  fi

  echo "$(date) [START] seed ${seed}" | tee "${log_file}"
  Rscript "'"${SIM_SCRIPT}"'" "${seed}" >> "${log_file}" 2>&1
  status=$?
  if [ "${status}" -eq 0 ]; then
    echo "$(date) [DONE]  seed ${seed} wrote ${out_file}" | tee -a "${log_file}"
  else
    echo "$(date) [FAIL]  seed ${seed} (exit ${status})" | tee -a "${log_file}"
    exit "${status}"
  fi
' sh {}

# --------- Summary ----------
end_epoch=$(date +%s)
elapsed=$((end_epoch - start_epoch))
printf "[INFO] All jobs finished. Elapsed: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

# Quick counts
total_seeds=$((END_SEED - START_SEED + 1))
done_count=$(find "${OUT_DIR}" -maxdepth 1 -name "seedVal_*_take3.RDS" | wc -l | tr -d ' ')
fail_logs=$(grep -l "\[FAIL\]" "${LOG_DIR}"/seed_*.log 2>/dev/null || true)

echo "[INFO] Total requested: ${total_seeds}"
echo "[INFO] Outputs present: ${done_count}"
if [ -n "${fail_logs}" ]; then
  echo "[WARN] Failures logged in:"
  echo "${fail_logs}"
else
  echo "[INFO] No failures detected."
fi