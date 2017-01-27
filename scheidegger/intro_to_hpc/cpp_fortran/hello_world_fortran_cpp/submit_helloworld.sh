#!/bin/bash -l

#BATCH --ntasks=1    ## how many cpus used here

#SBATCH --time=01:00:00  ## walltime requested

#SBATCH --output=slurm_test.out ## output file
#SBATCH --error=slurm_test.err    ## error 


### executable
./helloworld.exec

