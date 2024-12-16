#!/bin/bash
#BSUB -J stereo_midbrain_base
#BSUB -n 12
#BSUB -P acc_YuanLab
#BSUB -q gpuexpress
#BSUB -W 4:00
#BSUB -R "rusage[mem=10000] span[hosts=1]"
#BSUB -gpu num=1
#BSUB -R h100nvl

#BSUB -oo log/job_stereo_midbrain_base.out
#BSUB -eo log/job_stereo_midbrain_base.err

JOBID=$1

source /hpc/users/wangw32/.bash_profile

mkdir -p output log

conda activate ONTraC
ONTraC --meta-input raw_data/stereo_seq_mouse_mid_brain/stereo_input.csv --NN-dir output/stereo_midbrain_base_NN --GNN-dir output/stereo_midbrain_base_GNN --NT-dir output/stereo_midbrain_base_NT --n-cpu 12 --n-neighbors 50 --device cuda --epochs 1000 --batch-size 10 -s 42 --patience 100 --min-delta 0.001 --min-epochs 50 --lr 0.03 --hidden-feats 4 --n-gcn-layers 2 -k 6 --modularity-loss-weight 0.3 --purity-loss-weight 300 --regularization-loss-weight 0.1 --beta 0.03 > log/stereo_midbrain_base.log

bjobs -l -gpu $JOBID >log/job_stereo_midbrain_base_resource_usage.txt
