#!/bin/bash
#SBATCH --job-name=gp_minibatch         # Job name
#SBATCH --output=gp_minibatch_%j.out    # Output file (%j will be replaced by the job ID)
#SBATCH --error=gp_minibatch_%j.err     # Error file (%j will be replaced by the job ID)
#SBATCH --ntasks=32                     # Number of tasks (MATLAB is single-threaded unless using parallel toolbox)
#SBATCH --mem=128G                      # Memory per node (adjust as necessary)
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --partition=compute             # Partition name (adjust to your cluster partition)
#SBATCH --mail-type=END,FAIL            # Send email on job completion or failure
#SBATCH --mail-user=mfho@umich.edu      # Email address for notifications
#SBATCH --array=0-100%20                # Array for running 100 jobs, 20 jobs at a time

# Load MATLAB module if needed
module load matlab

# Define the number of quasars and offset for each batch
NUM_QUASARS=1000
QSOS_NUM_OFFSET=$(($SLURM_ARRAY_TASK_ID * $NUM_QUASARS))

# Run MATLAB script
matlab -nodisplay -nosplash -r "num_quasars=$NUM_QUASARS; qsos_num_offset=$QSOS_NUM_OFFSET; script_for_gp_minibatch; exit"
