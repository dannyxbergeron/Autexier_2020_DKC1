#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --mem=128000M
#SBATCH --array=[0-6]

module load samtools/1.5

namelist=('Testis1' \
'BRN_58' \
'Skeletal_muscle_1' \
'Skeletal_muscle_2' \
'Skeletal_muscle_3' \
'BRN_574' \
)


name=${namelist[$SLURM_ARRAY_TASK_ID]}
project=Tissues
project_path=$SCRATCH/sequencing_runs/$project

mkdir -p $project_path/samtools_idxstats/

bamfile=$project_path/STAR/ERCC/${name}/Aligned.sortedByCoord.out.bam
outputfile=$project_path/samtools_idxstats/${name}_idxstat.tsv

echo ${name};

#index inputs
samtools index $bamfile &&
samtools idxstats $bamfile > $outputfile 


echo 'samtools idxstats done'
