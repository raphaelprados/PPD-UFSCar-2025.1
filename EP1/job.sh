#!/bin/bash
#SBATCH -J mmul                     # Job name
#SBATCH -p fast                     # Job partition
#SBATCH -n 1                        # Number of processes
#SBATCH -t 01:30:00                 # Run time (hh:mm:ss)
#SBATCH --cpus-per-task=120         # Number of CPUs per process
#SBATCH --output=%x.%j.out          # Name of stdout output file - %j expands to jobId and %x to jobName
#SBATCH --error=%x.%j.err           # Name of stderr output file
#SBATCH --mail-type=END             # Notifies the user that the job finished execution
#SBATCH --mail-user=raphael.santos25@estudante.ufscar.br

echo "*** SEQUENTIAL ***"
srun singularity run container.sif pi_seq 1000000000 >> pi_seq.txt

echo "*** PTHREAD ***"
srun singularity run container.sif pi_pth 1000000000 1 >> pi_pth_1.txt
srun singularity run container.sif pi_pth 1000000000 2 >> pi_pth_2.txt
srun singularity run container.sif pi_pth 1000000000 5 >> pi_pth_5.txt
srun singularity run container.sif pi_pth 1000000000 10 >> pi_pth_10.txt
srun singularity run container.sif pi_pth 1000000000 20 >> pi_pth_20.txt
srun singularity run container.sif pi_pth 1000000000 40 >> pi_pth_40.txt
srun singularity run container.sif pi_pth 1000000000 120 >> pi_pth_120.txt

diff pi_seq.txt pi_pth_1.txt
diff pi_seq.txt pi_pth_2.txt
diff pi_seq.txt pi_pth_5.txt
diff pi_seq.txt pi_pth_10.txt
diff pi_seq.txt pi_pth_20.txt
diff pi_seq.txt pi_pth_40.txt
diff pi_seq.txt pi_pth_120.txt

# OBS: if it is an MPI job
# use --mpi=pmi2
# srun --mpi=pmi2 singularity run container.sif mmul_seq 1000
