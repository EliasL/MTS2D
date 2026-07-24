#!/usr/bin/env bash
# Recovery run for the expensive L=500 first-noisy edge-flipping minimizations.
# Submit from the MTS2D repository root with:
#   mkdir -p benchmark_jobs
#   sbatch benchmarks/slurm_cauchy_l500_first_minimization_recovery.sh

#SBATCH --job-name=mts2d-l500-firstmin
#SBATCH --partition=LocalQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=benchmark_jobs/%x-%j.out
#SBATCH --error=benchmark_jobs/%x-%j.err

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "This benchmark must be launched by Slurm (sbatch)." >&2
    exit 2
fi

project_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "$project_dir"

result_dir="${BENCHMARK_OUTPUT_DIR:-$project_dir/benchmark_results/cauchy_l500_firstMinRecovery_${SLURM_JOB_ID}}"
mkdir -p "$result_dir"

echo "MTS2D L=500 first-minimization recovery benchmark"
echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "Project: ${project_dir}"
echo "Results: ${result_dir}"
echo "Slurm allocation:"
scontrol show job "$SLURM_JOB_ID" --oneliner
echo "CPU topology:"
lscpu | grep -E '^(CPU\(s\)|On-line CPU|Thread\(s\) per core|Core\(s\) per socket|Socket\(s\)|NUMA node\(s\))'

export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1

srun \
    --ntasks=1 \
    --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --cpu-bind=cores \
    --distribution=block:block \
    python3 tools/run_benchmarks.py \
        --preset full \
        --experiments workload \
        --workload-sizes 500 \
        --workload-loads 0.15 \
        --fixture-seed-stride 0 \
        --no-generate-load-015 \
        --first-minimization-only \
        --first-minimization-histories edgeFlipping \
        --allow-available-fixture-fallback \
        --large-initial-condition-budget-minutes 120 \
        --total-budget-hours 11.75 \
        --build-jobs 16 \
        --output-dir "$result_dir"

report_path=$(find "$result_dir" -name report.md -type f -print -quit)
echo "Finished: $(date --iso-8601=seconds)"
echo "Report: ${report_path:-not generated}"
