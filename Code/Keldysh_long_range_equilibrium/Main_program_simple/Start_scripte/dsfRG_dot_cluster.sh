#!/bin/bash
#SBATCH -o /home/hpc/uh3o1/ri26yad/fRG/Keldysh_long_range_equilibrium/Aufgeraeumt/Main_program_simple/LOG/myjob.%j.%N.out
#SBATCH -D /home/hpc/uh3o1/ri26yad/fRG/Keldysh_long_range_equilibrium/Aufgeraeumt/Main_program_simple/RUN
#SBATCH -J some_job
#SBATCH --get-user-env
#SBATCH --clusters=lmu_asc
#SBATCH --partition=lmu_asc_shm8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16000
#SBATCH --mail-type=end
#SBATCH --mail-user=lukas.weidinger@physik.lmu.de
#SBATCH --export=NONE
#SBATCH --time=240:00:00

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=1

. /etc/profile

#ulimit -c unlimited

cd /naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple 

/home/hpc/uh3o1/ri26yad/fRG/Keldysh_long_range_equilibrium/Aufgeraeumt/Main_program_simple/Compiled/dsfRG_dot 
