# Sbatch script was developed by Molly Richardson <UP824687__at__myport.ac.uk>
# as part of an MPhys Project with Daniel Thomas

#!/bin/bash

#SBATCH --job-name=fireflypar
#SBATCH --output=/users/mr824687/array_%A_%a.out
#SBATCH --error=/users/mr824687/array_%A_%a.err
#SBATCH --partition=sciama2.q
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-73

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# The lines below run the computations

module add anaconda3

python ./firefly_release-0.1.1/firefly_fits_par.py -v ${SLURM_ARRAY_TASK_ID}
