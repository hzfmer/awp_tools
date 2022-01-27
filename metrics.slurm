#!/bin/bash
#SBATCH -J metrics       # Job name
#SBATCH -o metrics_%j.out          # Name of stdout output file
#SBATCH -e metrics_%j.err          # Name of stderr error file
#SBATCH -N 30              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 15:00:00        # Run time (hh:mm:ss)
#SBATCH --mem=0
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

echo "Running metrics, Starting at `date`"
cat $0
# Arguments
# -n: number of CPUs
# -f: do filtering, given by `-l` and `-h`
# -l: lower bounds
# -h: upper bounds
# -g: compute SAs; the frequencies ar hardcoded in source code.
# -x: partition along x direction
# -y: partition along y direction 
# -z: partition along z direction
# -p: finer partition; should be divisible by number of nodes along y direction
srun -n 960  ~/scratch/awp_tools/bin_andes/metrics2 -f -l 0.15,0.15,0.15,1,2.5 -h 1,2.5,5,2.5,5 -g -x 48 -y 20 -z 1 -p 3

echo "Ending at `date`"
