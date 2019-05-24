#!/bin/bash
#
# Project identification
#
#SBATCH --account=JOBACC_TMP  # no need, this is setup in ~/.bashrc
#
# MAX. JOB LENGTH
#
#SBATCH --time=JOBTIME_TMP       # time (DD-HH:MM:SS)
#
# NODES, CORES/NODE
#
#SBATCH --ntasks=NTASKS_TMP
#
#SBATCH --mem-per-cpu=MEM_CPU_TMP      # memory; default unit is megabytes
#
# ENVIRONMENT VARIABLES
#
#SBATCH --export=ALL
#
# DO NOT RESTART AUTOMATICALLY
#
#SBATCH --no-requeue
#
# JOB NAME
#
#SBATCH --job-name=JOBNAME_TMP
#
# EMAIL JOB RESULTS
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user="USER_EMAIL_TMP"
#
# LOG FILE NAMES
#
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# LOAD NETCDF LIBRARY
module purge
module load hdf5-mpi
module load netcdf-mpi
module load netcdf-fortran-mpi

# SWITCH TO PROJECT DIRECTORY

# tell wrapper script that job started
echo "job started" > jobStart

# SEND JOB
srun oceanM OCEANFILE_TMP > LOGFILE_TMP

# tell wrapper script job completed
echo "job ended" > jobEnd

exit
