#!/bin/bash
#SBATCH -J add_het  # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 78              # Total # of nodes
#SBATCH -t 00:10:00        # Run time (hh:mm:ss)
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
NZ=125
PX=54
PY=46
NCPU=$((PX * PY))
NZ_SSH=125
NZ_TAP=938
VPMAX=9500.
VSMAX=6000.

std=5
hurst=5
l=100
s=$(awk -v dividend="${std}" 'BEGIN {printf "%.2f", dividend/100; exit(0)}')
std=$(printf "%02d" $std)
h=$(printf "%03d" $hurst)
FH_HOM=../cvm/la_habra_large_cvmsi_8m.media
FH_SSH=ssh_8m_large_r05h${h}l${l}.out
#FH_HET=../la_habra_large_100120/mesh_large_8m_s${std}h${h}l${l}.bin
FH_HET=test_fh_het.bin

#srun -n ${NCPU} --cpu-bind=threads --distribution=cyclic  ./add_het --nx ${NX} --ny ${NY} --nz ${NZ} --std ${s} --px ${PX} --py ${PY} --nz_ssh ${NZ_SSH} --nz_tap ${NZ_TAP} --vpmax ${VPMAX} --vsmax ${VSMAX} --fh_ssh ${FH_SSH} --fh_hom ${FH_HOM} --fh_het ${FH_HET}
srun -n ${NCPU} -N 78 --cpu-bind=threads --distribution=cyclic ./add_het --nx ${NX} --ny ${NY} --nz ${NZ} --std ${s} --px ${PX} --py ${PY} --nz_ssh ${NZ_SSH} --nz_tap ${NZ_TAP} --vpmax ${VPMAX} --vsmax ${VSMAX} --fh_ssh ${FH_SSH} --fh_hom ${FH_HOM} --fh_het ${FH_HET}
echo -e "\nEnding at `date`"
