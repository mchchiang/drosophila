#!/bin/bash

#$ -N RUNNAME  # job name

#$ -V           # use all shell environment variables

# Choose a queue:
#$ -q cdt.7.day     # cdt.7.day, cm.7.day, sopa.1.day 

# Set job runtime
#$ -l h_rt=24:00:00   # time limit = 7days  (in sopa.1.day the time limit = 24h)

# Request memory necessary to run the job
#$ -l h_vmem=1G

# Choose a parallel environment:
#$ -pe mpi 64   # max number of processors = 64

# Path to error and output files
#$ -e /Disk/ds-sopa-personal/s1309877/PhD_Project/drosophila/tmp/
#$ -o /Disk/ds-sopa-personal/s1309877/PhD_Project/drosophila/tmp/

if [ ! -d /scratch/s1309877 ]; then
mkdir /scratch/s1309877
fi

rsync -r /Disk/ds-sopa-personal/s1309877/PhD_Project/drosophila/RUNNAME /scratch/s1309877/   # Copy dir with lammps input and script files

cd /scratch/s1309877/RUNNAME

module load mpi/openmpi-x86_64

/usr/lib64/openmpi/bin/mpirun -np 64 /Disk/ds-sopa-personal/s1309877/bin/lammps-16Mar18/src/lmp_mpi -in LAMMPSSCRIPT -log LOGFILE

cd ..
rsync -r RUNNAME /Disk/ds-sopa-personal/s1309877/PhD_Project/drosophila/RUNNAME_finished

rm -r RUNNAME

exit 0
