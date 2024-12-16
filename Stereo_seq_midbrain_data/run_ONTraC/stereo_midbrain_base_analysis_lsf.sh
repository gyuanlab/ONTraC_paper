#!/bin/bash
#BSUB -J stereo_midbrain_base_analysis
#BSUB -n 2
#BSUB -P acc_YuanLab
#BSUB -q express
#BSUB -W 2:00
#BSUB -R "rusage[mem=10000] span[hosts=1]"

#BSUB -oo log/job_stereo_midbrain_base_analysis.out
#BSUB -eo log/job_stereo_midbrain_base_analysis.err

JOBID=$1

source /hpc/users/wangw32/.bash_profile

mkdir -p analysis_output

conda activate ONTraC

ONTraC_analysis --NN-dir output/stereo_midbrain_base_NN --GNN-dir output/stereo_midbrain_base_GNN --NT-dir output/stereo_midbrain_base_NT -o analysis_output/stereo_midbrain_base -l log/stereo_midbrain_base.log -s --suppress-cell-type-composition
