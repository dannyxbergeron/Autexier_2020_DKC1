#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --mem=128000M
#SBATCH --cpus-per-task=10
#SBATCH --get-user-env
#SBATCH --mail-user=danny.bergeron@usherbrooke.ca
#SBATCH --mail-type=ALL


module load nixpkgs/16.09  intel/2016.4
module load bcl2fastq2/2.20


sample_sheet_name=hnRNP_NOP58_NHP2L1_8_samples.csv
bcl2_path=180709_NB502083_0023_AH2JH3BGX7

mkdir -p data/reads


bcl2fastq -r 10 -w 10 -R $bcl2_path -o data/reads --sample-sheet $bcl2_path/$sample_sheet_name --no-lane-splitting
