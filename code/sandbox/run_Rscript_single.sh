#!/bin/bash
#SBATCH --time=0-02:00:00   # Maximum time limit
#SBATCH --ntasks=1          # Number of tasks
#SBATCH --cpus-per-task=1   # Number of CPU cores per task
#SBATCH --mem=10G            # Memory per node
#SBATCH --partition=large_336
#SBATCH --array=1-5
echo "R is about to run"
Rscript single_matrix.R
echo "R has finished running"
# this is a comment at the end of the file
