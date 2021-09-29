#!/bin/bash
#
#SBATCH --job-name=jobname
#SBATCH --mem=2040
#SBATCH --time=2-20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<username>@physik.uni-muenchen.de
#SBATCH --chdir=/dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/<lrz-ID>/mfrg/Keldysh_mfRG/
#SBATCH --output=/dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/<lrz-ID>/mfrg/Keldysh_mfRG/runs/jobname.%j.%N.out
#SBATCH --error=/dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/<lrz-ID>/mfrg/Keldysh_mfRG/runs/jobname.%j.%N.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=32

echo $HOSTNAME
echo $SLURM_ARRAY_JOB_ID
echo $SLURM_NTASKS

export OMP_NUM_THREADS=32
mpiexec -n $SLURM_NTASKS ./main.o
