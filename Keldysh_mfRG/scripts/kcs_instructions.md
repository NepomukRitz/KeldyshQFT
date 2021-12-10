# Introduction to the KCS-Cluster

Project name for LS vonDelft: **pn34vu**

### Overview

- Login node: kcs-login (`$ ssh <lrz-ID>@kcs-login.cos.lrz.de`)
- Need to have an established VPN connection to the MWN (or be physically connected to it).
- slurm jobs can only request whole machines -> **Always use all 32 cores on each node!**
- 150 nodes, 32 cores/node on 2 sockets
- memory configurations: 180 GB, 370 GB, 760 GB (2x)
- standard queue: runtime: 72h
- long-runner queue: runtime: 30 days, only 1 job/user (Slurm: partition kcs_long, see below)


### File System

 - GPFS file system, efficient for large files 
   - page size 16 MB -> **Do not use many small files!**
   - quota: 100 TB for the chair 
   
<br />

- LRZ DSS (data science storage), organized in containers 
    - LRZ documentation:
      https://doku.lrz.de/display/PUBLIC/Data+Science+Storage
    - Web interface for storage management: https://dssweb.dss.lrz.de (only accessible for chair admins)
    - default container 0000
      - quota (can be changed via web interface -> ask chair admin): \
      1 TB for chair, 64 GB/user, 65000 files/user
      - daily backup
      - path of "home" directory:
        `/dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/\<lrz-ID>`
    - can create new containers, some settings cannot be changed afterwards (see LRZ
      documentation)
        - useful for new projects that need a lot of storage (ask chair admin)
    - share containers via globus (https://doku.lrz.de/display/PUBLIC/DSS+How+Globus+Data+Transfer+and+Globus+Sharing+for+DSS+works)


### Setting up a .bashrc file
It is highly recommended to set up a `.bashrc` file to configure shortcuts 
and load some modules automatically on login.\
In the home-directory entered when logging onto KCS (`/dss/dsshome1/lxc09/<lrz-ID>`), create a `.bashrc`
file using a command-line text editor like vim or nano (nano has to be loaded first using 
the module system, see below).\
A good starting point is the following file (comments start with a `#`):

```
# set up an alias for the "home" directory used for computations
alias kcshome="cd /dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/<lrz-ID>"

# useful to get all information about all files in a given directory in human-readable form
alias ll="ls -lah"

# simplifies navigation through the command-line history
bind '"\e[A": history-search-backward'
bind '"\e[B": history-search-forward'

# load some modules by default
module load nano                # Beginner-friendly command-line text editor
module load gsl                 # TODO: Better to be included in the compile script!
module load mkl                 # TODO: Do we need this?
module load boost/1.61_icc      # TODO: Warning that this module is scheduled for retirement by end of 2019 (!)

module load hdf5/1.8.20-cxx-frt-threadsafe # NOT SUFFICIENT to load in compile script! Needed here as well!
```
After creating or modifying the `.bashrc` file, one has to reload it using `# source .bashrc`. 
Alternatively, one can log off an log on to KCS again, 
as the `.bashrc` file is always automatically loaded upon login.

### Module System

Grants access to pre-installed packages. Useful commands:

| Command | Explanation |
| ------- | ----------- |
| `module avail` | shows available modules |
| `module list`  | shows the currently loaded modules |
| `module load <module_name>`  | loads the module \<module_name> |
| `module unload <module_name>`  | unloads the module \<module_name> |

### Cloning the code (to the right place)

The code should be placed into the "home" directory used for computations, for which the shortcut `kcshome`
was set up in the `.bashrc` file above. It is then simply a matter of navigating to the "home" 
directory and cloning the correct git repository:
```
$ kcshome
$ git clone https://gitlab.physik.uni-muenchen.de/LDAP_ls-vondelft/<project>.git
```
It is easiest to clone the git repo via https and not ssh. 
Otherwise one would have to set up a ssh-key-pair.



### Compiling the code

Depending on the setup used, there should already be a makefile or a compile script in the code base
of the git repository. As of September 2021, the compile script for the Keldysh mfRG code is located at
(starting from `kcshome`) `mfrg/Keldysh_mfRG/scripts/compile_kcs.sh` and reads

```
#!/bin/bash
#  environment variable KELDYSH_MFRG needs to point to the "Keldysh_mfRG" directory of the repository:
#  in ~/.bashrc:
#  export KELDYSH_MFRG="/dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/<lrz-ID>/mfrg/Keldysh_mfRG"

module load gcc
module load hdf5/1.8.20-cxx-frt-threadsafe
module load fftw
module load gsl
module load boost/1.61_icc

export LANG=C
export LC_ALL=C

HDF5="$HDF5_INC $HDF5_CPP_SHLIB $HDF5_SHLIB $SZIP_LIB -lz"
FFTW="$FFTW_INC $FFTW_LIB"
GSL="$GSL_INC $GSL_LIB"
BOOST="$BOOST_INC -L$BOOST_LIBDIR$"

mpiCC -std=c++17 $KELDYSH_MFRG/main.cpp -o $KELDYSH_MFRG/main.o -fopenmp $FFTW $HDF5 $GSL $BOOST
```

As one can already tell from this file, before it can be executed, one has to define another shortcut
to set the correct absolute path to the directory of the code.
To do that, access the `.bashrc` file via `$ nano ~/.bashrc` and add the line 

`export KELDYSH_MFRG="/dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/<lrz-ID>/mfrg/Keldysh_mfRG"`

The compile script can then simply be executed by 

`$ ./mfrg/Keldysh_mfRG/scripts/compile_kcs.sh`.

### Submitting jobs

The job handling of the cluster is organized by SLURM. For basic information on SLURM see https://www.en.it.physik.uni-muenchen.de/dienste/rechencluster/index.html.
To submit a job, one needs to provide 
a corresponding shell script, e.g. by modifying `mfrg/Keldysh_mfRG/scripts/batchfile.sh` (see official Slurm documentation https://slurm.schedmd.com/sbatch.html for details):

```
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

echo $HOSTNAME
echo $SLURM_ARRAY_JOB_ID
echo $SLURM_NTASKS

export OMP_NUM_THREADS=32
mpiexec -n $SLURM_NTASKS ./main.o
```

#### Description:

| Option | Explanation |
| ------ | ----------- |
| `job-name` | Name of job in slurm queue. |
| `mem` | Requested memory (minimum) in MB. |
| `time` | Job wall time: after this time, the job will be killed by slurm. Maximum runtime is 3 days. |
| `mail-type` | Settings for slurm status emails (see below). |
| `chdir` | Directory in which the batchfile will be executed (=path to executable). |
| `output` | Name of log file. `%j` = job-ID, `%N` = node-ID. |
| `error` | Name of error log file. |
| `nodes` | Number of (MPI) nodes on which to execute the job. |
| `ntasks-per-node` | Number of MPI tasks per node. When using OpenMP, set to 1 (only use MPI for inter-node parallelization). |
| `OMP_NUM_THREADS` | Number of OpenMP threads. Should be equal to the number of available cores per node (=32 on KCS). |
| `main.o` | Job executable. |

#### Optional settings:

`#SBATCH --partition=kcs_long` : Submit the job to the KCS long-runner queue (wall time < 30 days).

#### SLURM status emails:

If activated, Slurm informs per email when jobs start/finish/fail etc. Deactivate Slurm emails if you are submitting many jobs, and do NOT forward these emails to another email account (see https://www.en.it.physik.uni-muenchen.de/dienste/rechencluster/index.html).

### Useful SLURM commands

| Command | Explanation |
| ------- | ----------- |
| `sinfo` | view information about SLURM nodes and partitions.  |
| `squeue`  | show all information about pending and running jobs |
| `squeue -u <lrz-ID>`  | show all information about pending and running jobs of \<lrz-ID> |
| `sbatch batchfile.sh`  | submit a job configured in `batchfile.sh`  |
| `scancel <job-ID>`  | cancel the job \<job-ID> |

### Accessing the data

Data can be downloaded from KCS to `<local-path>` on your workstation/laptop using rsync:

`rsync -auvh <lrz-ID>@kcs-login.cos.lrz.de:/dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/<lrz-ID>/mfrg/Keldysh_mfRG/<path-to-result-file> <local-path>`

Hint: if you have defined e.g.

`kcs=<lrz-ID>@kcs-login.cos.lrz.de` and

`kcshome=/dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/<lrz-ID>/mfrg/Keldysh_mfRG`

in your local .bashrc (on the workstation/laptop), the command simplifies:

`rsync -auvh $kcs:$kcshome/<path-to-result-file> <local-path>`.

(If the command is executed from within the destination directory `<local-path>`, then simply replace `<local-path>`=`.` )

### Running unit tests

As long as unit test really are just unit test and do not take much time to execute, 
it should be alright to run them even on the login node. 
