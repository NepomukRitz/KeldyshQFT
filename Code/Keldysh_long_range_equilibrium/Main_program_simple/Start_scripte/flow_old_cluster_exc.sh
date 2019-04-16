#!/bin/bash
#SBATCH -o /home/hpc/uh3o1/ri26yad/fRG/Keldysh_long_range_equilibrium/Aufgeraeumt/Main_program_simple/LOG/myjob.%j.%N.out
#SBATCH -D /home/hpc/uh3o1/ri26yad/fRG/Keldysh_long_range_equilibrium/Aufgeraeumt/Main_program_simple/RUN
#SBATCH -J some_job
#SBATCH --get-user-env
#SBATCH --clusters=lmu_exc
#SBATCH --partition=lmu_exc_uv2k_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=15000
#SBATCH --mail-type=end
#SBATCH --mail-user=lukas.weidinger@physik.lmu.de
#SBATCH --export=NONE
#SBATCH --time=20:00:00

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=1

. /etc/profile

#ulimit -c unlimited

cd /naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple 

/home/hpc/uh3o1/ri26yad/fRG/Keldysh_long_range_equilibrium/Aufgeraeumt/Main_program_simple/Compiled/flow_old 
