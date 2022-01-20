#!/bin/sh 
#
########## Begin MOAB/Slurm header ##########
#
# Give job a reasonable name
#SBATCH -J 3D-continuousReference_mu_8
#
# Request number of nodes and CPU cores per node for job
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#
# Request memory
#SBATCH --mem=64000mb
#
# Estimated wallclock time for job
#SBATCH --time=24:00:00
#
#
# Send mail when job begins, aborts and ends
#SBATCH --mail-type=ALL
#
########### End MOAB header ##########

echo "Submit Directory:                     $MOAB_SUBMITDIR"
echo "Working Directory:                    $PWD"
echo "Running on host                       $HOSTNAME"
echo "Job id:                               $MOAB_JOBID"
echo "Job name:                             $MOAB_JOBNAME"
echo "Number of nodes allocated to job:     $MOAB_NODECOUNT"
echo "Number of cores allocated to job:     $MOAB_PROCCOUNT"

# Load module 
module load math/matlab/R2020a

# Start a Matlab program
matlab -nodisplay -r "continuous3DReference_mu_8" >  3D-continuousReference_mu_8.out 2>&1

exit
