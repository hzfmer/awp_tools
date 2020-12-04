#!/bin/bash
#SBATCH -J add_het  # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 54              # Total # of nodes
#SBATCH -t 00:15:00        # Run time (hh:mm:ss)
#SBATCH --mem=0            # Use all available memory 
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

module load intel
module unload Darshan

echo "Running add_het, Starting at `date`"
cat $0
# -n = npx * npy, divisible by nx and ny, respectively
NX=19440
NY=14904
NZ=1250
PX=48
PY=36
NIO=250
NCPU=$((PX * PY))
NZ_SSH=1250
NZ_TAP=938
VPMAX=9500.
VSMAX=6000.

std=10
hurst=5
l=500
s=$(awk -v dividend="${std}" 'BEGIN {printf "%.2f", dividend/100; exit(0)}')
std=$(printf "%02d" $std)
h=$(printf "%03d" $hurst)
FH_HOM=../cvm/la_habra_large_cvmsi_8m.media
FH_SSH=ssh_8m_large_r05h${h}l${l}.out
FH_HET=../la_habra_large_100120/mesh_large_8m_s${std}h${h}l${l}.bin

#srun -n ${NCPU} --cpu-bind=threads --distribution=cyclic  ./add_het --nx ${NX} --ny ${NY} --nz ${NZ} --std ${s} --px ${PX} --py ${PY} --nz_ssh ${NZ_SSH} --nz_tap ${NZ_TAP} --vpmax ${VPMAX} --vsmax ${VSMAX} --fh_ssh ${FH_SSH} --fh_hom ${FH_HOM} --fh_het ${FH_HET}
srun -n ${NCPU} -N 54 ./add_het_nio --nx ${NX} --ny ${NY} --nz ${NZ} --std ${s} --px ${PX} --py ${PY} --nz_ssh ${NZ_SSH} --nz_tap ${NZ_TAP} --nio ${NIO} --vpmax ${VPMAX} --vsmax ${VSMAX} --fh_ssh ${FH_SSH} --fh_hom ${FH_HOM} --fh_het ${FH_HET}
echo -e "\nEnding at `date`"
