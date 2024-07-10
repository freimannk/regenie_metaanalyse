#!/bin/bash

#SBATCH --time=96:45:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="regenie_meta"
#SBATCH --partition=amd

module load any/jdk/1.8.0_265
module load nextflow
module load any/singularity/3.5.3
module load squashfs/4.4

NXF_VER=22.04.3 nextflow run main.nf -profile tartu_hpc -resume  