#!/bin/bash

#$ -N job.name  # job name

#$ -V           # use all shell environment variables

# Choose a queue:
#$ -q cdt.7.day     # cdt.7.day, cm.7.day, sopa.1.day 

# Set job runtime
#$ -l h_rt=168:00:00   # time limit = 7days  (in sopa.1.day the time limit = 24h)

# Request memory necessary to run the job
#$ -l h_vmem=1G

# Choose a parallel environment:
#$ -pe mpi 32   # max number of processors = 64

# Path to error and output files
#$ -e /Disk/ds-sopa-personal/s1460633/PhD_project/Drosophila/tmp/
#$ -o /Disk/ds-sopa-personal/s1460633/PhD_project/Drosophila/tmp/

replica=1

if [ ! -d /scratch/Carolina ]; then
mkdir /scratch/Carolina
fi

cp -r /Disk/ds-sopa-personal/s1460633/PhD_project/Drosophila/job.dir.$replica /scratch/Carolina/.   # Copy dir with lammps input and script files

cd /scratch/Carolina/job.dir.$replica

module load mpi/openmpi-x86_64

/usr/lib64/openmpi/bin/mpirun -np 32 /Disk/ds-sopa-personal/s1460633/PhD_project/lammps-14May16-harmlj/src/lmp_mpi -in LAMMPS-script.lam

cd ..
cp -r job.dir.$replica /Disk/ds-sopa-personal/s1460633/PhD_project/Drosophila/job.dir.$replica-finished

rm -r job.dir.$replica

exit 0
