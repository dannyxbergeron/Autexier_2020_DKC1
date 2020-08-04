#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --mem=128000M
#SBATCH --array=[0-7]

module load picard/2.18.9

namelist=('decoy' \
'Testis1' \
'BRN_58' \
'Skeletal_muscle_1' \
'Skeletal_muscle_2' \
'Skeletal_muscle_3' \
'BRN_574' \
)

name=${namelist[$SLURM_ARRAY_TASK_ID]}
project=Tissues
project_path=$SCRATCH/sequencing_runs/$project

mkdir -p $project_path/Picard/
mkdir -p $project_path/Picard/$name

BAM_FILE=$project_path/STAR/$name/Aligned.sortedByCoord.out.bam

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
      I=$BAM_FILE \
      O=$project_path/Picard/$name/Picard_insert_size_metrics.txt \
      H=$project_path/Picard/$name/Picard_size_histogram.pdf \

