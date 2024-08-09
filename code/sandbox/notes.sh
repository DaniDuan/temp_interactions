ssh danica@harvey.dept.imperial.ac.uk # Log in

sftp danica@harvey.dept.imperial.ac.uk # file upload

put sim_frame.jl # throw in files
put sim_frame.jl code # throw in files to directory 

get sim_frame.jl # pull out files

sbatch xxx.sh
###### the bash file #######

#!/bin/bash
#SBATCH --time=0-00:30:00   # Maximum time limit (30 minutes)
#SBATCH --ntasks=1          # Number of tasks
#SBATCH --cpus-per-task=1   # Number of CPU cores per task
#SBATCH --mem=10G            # Memory per node
#SBATCH --partition=large_336
#SBATCH --array=1-500 

echo "Julia is about to run"
julia xxx.jl
echo "Julia has finished running"

# at the end of sh file

### example in R
iter <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
### in Julia
index_str = ENV["SLURM_ARRAY_TASK_ID"]
index = parse(Int, index_str)


# https://slurm.schedmd.com/quickstart.html
# submitting a batch
# sbatch is used to submit a job script for later execution. The script will typically contain one or more srun commands to launch parallel tasks.
sbatch temp_div.sh

# scancel is used to cancel a pending or running job or job step. It can also be used to send an arbitrary signal to all processes associated with a running job or job step.
scancel # job number

# checking in the queue
# squeue reports the state of jobs or job steps. It has a wide variety of filtering, sorting, and formatting options. By default, it reports the running jobs in priority order and then the pending jobs in priority order.
squeue

# adding in packages in Julia
# Pkg.add(["Distributions", "LinearAlgebra", "DifferentialEquations", "Plots", "StatsPlots","Sundials","Parameters","CSV","DataFrames","CairoMakie","LsqFit","Logging","TOML"])