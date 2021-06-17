#!/bin/bash

#PBS -l walltime=00:00:20,nodes=7:ppn=4
#PBS -N task
#PBS -q batch

cd $PBS_O_WORKDIR
mpirun --hostfile $PBS_NODEFILE -np 28 ./NumSend
