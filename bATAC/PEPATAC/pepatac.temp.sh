#!/bin/sh
#PBS -q hotel
#PBS -N pepatac
#PBS -V
#PBS -M rlancione@ucsd.edu
#PBS -m abe
#PBS -e /home/rlancione/jeo/pepatac.sh.e
#PBS -o /home/rlancione/jeo/pepatac.sh.o
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=4


# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pepatac 

# Yaml path to replace
path_to_yaml="yaml_path"

# Run pepatac
looper run $path_to_yaml

looper runp $path_to_yaml

looper report $path_to_yaml


