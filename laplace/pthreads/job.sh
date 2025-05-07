#!/bin/bash
#SBATCH -J laplace_ppd              # Job name
#SBATCH -p fast                     # Job partition
#SBATCH -n 1                        # Number of processes
#SBATCH -t 01:30:00                 # Run time (hh:mm:ss)
#SBATCH --mem=8G                    # Memory requirements
#SBATCH --cpus-per-task=120         # Number of CPUs per process
#SBATCH --output=%x.%j.out          # Name of stdout output file - %j expands to jobId and %x to jobName
#SBATCH --error=%x.%j.err           # Name of stderr output file
#SBATCH --mail-type=END             # Notifies the user that the job finished execution
#SBATCH --mail-user=raphael.santos25@estudante.ufscar.br

size=1000

echo "***************** SEQUENTIAL *****************"

srun singularity run container.sif seq "$size"
srun singularity run container.sif seq_lin "$size"
srun singularity run container.sif seq_ptr "$size"
srun singularity run container.sif seq_linptr "$size"

echo "***************** PARALLEL *****************"

echo "---------------- DEMAND ----------------"

echo "----------- regular -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif d_laplace "$size" "$run"
done

echo "----------- linear -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif d_laplace_lin "$size" "$run"
done

echo "----------- pointer -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif d_laplace_ptr "$size" "$run"
done

echo "----------- linear and pointer -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif d_laplace_linptr "$size" "$run"
done


echo "---------------- REGULAR ----------------"

echo "----------- regular -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif r_laplace "$size" "$run"
done

echo "----------- linear -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif r_laplace_lin "$size" "$run"
done

echo "----------- pointer -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif r_laplace_ptr "$size" "$run"
done

echo "----------- linear and pointer -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif r_laplace_linptr "$size" "$run"
done


echo "---------------- ROUND ROBIN ----------------"

echo "----------- regular -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif rr_laplace "$size" "$run"
done

echo "----------- linear -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif rr_laplace_lin "$size" "$run"
done

echo "----------- pointer -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif rr_laplace_ptr "$size" "$run"
done

echo "----------- linear and pointer -----------"
for run in 1 5 10 20 40 80
do
    srun singularity run container.sif rr_laplace_linptr "$size" "$run"
done

# OBS: if it is an MPI job
# use --mpi=pmi2
# srun --mpi=pmi2 singularity run container.sif mmul_seq 1000
