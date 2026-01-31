#!/bin/bash
#SBATCH -A [proj name]
#SBATCH -p core -n 4
#SBATCH -t 0-24:00:00
#SBATCH -J Est-sfs


module load conda
source conda_init.sh
conda activate [CONDA ENV with estsfs]

mkdir -p output_allsites

## See Est-sfs manual for input and config data
data_input=$1
output_1=$2
output_2=$3

est-sfs config.txt $data_input seedfile.txt $output_1 $output_2


mv $output_1 $output_2 output_allsites
