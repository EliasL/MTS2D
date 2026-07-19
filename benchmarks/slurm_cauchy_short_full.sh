#!/usr/bin/env bash
# First short integration run of the full benchmark matrix on Cauchy.
# Submit from the MTS2D repository root with:
#   mkdir -p benchmark_jobs
#   sbatch benchmarks/slurm_cauchy_short_full.sh

#SBATCH --job-name=mts2d-bench-short-full
#SBATCH --partition=LocalQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=benchmark_jobs/%x-%j.out
#SBATCH --error=benchmark_jobs/%x-%j.err

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "This benchmark must be launched by Slurm (sbatch)." >&2
    exit 2
fi

project_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "$project_dir"

result_dir="${BENCHMARK_OUTPUT_DIR:-$project_dir/benchmark_results/cauchy_short_full_${SLURM_JOB_ID}}"
mkdir -p "$result_dir"

echo "MTS2D short full-matrix benchmark"
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
        --calls 2 \
        --min-calls 2 \
        --max-calls 2 \
        --warmup-calls 1 \
        --repetitions 2 \
        --workload-sizes 16,32 \
        --case-budget-minutes 5 \
        --total-budget-hours 0.75 \
        --build-jobs 16 \
        --output-dir "$result_dir"

echo "Finished: $(date --iso-8601=seconds)"
echo "Report: ${result_dir}/report.md"
