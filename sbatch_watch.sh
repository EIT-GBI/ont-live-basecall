#!/bin/bash
#SBATCH --job-name=test_part2

#SBATCH --partition=gpu

#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --gres=gpu:1

#SBATCH --output=slurm-%j.stdout
#

apptainer exec --bind /mnt:/mnt /mnt/gbi-shared/home/pauline-eitgbi/ont-watcher_v2.sif bash app_exc_watch.sh
