#!/bin/bash
#SBATCH --time=0-02:00:00   # Maximum time limit
#SBATCH --ntasks=1          # Number of tasks
#SBATCH --cpus-per-task=1   # Number of CPU cores per task
#SBATCH --mem=1G            # Memory per node
#SBATCH --partition=large_336
#SBATCH --array=1-999

echo "Julia is about to run"
julia sim_Harvey_eff.jl
echo "Julia has finished running"

