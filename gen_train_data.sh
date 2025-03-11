#!/bin/bash

# generate the training dataset for YOLOOP

# configurations
# data_dir holds the .mcool file for hi-c contact matrix
data_dir="/Dataset/HiC/hic"
# loop_dir holds the .bedpe file for positive chromatin interaction
loop_dir="/Dataset/HiC/hic/loop_train"
# Take CTCF interactions as example
output_dir="/Dataset/YOLOOP_Train/CTCF"

# set other parameters to default
# window_size=512
box_size=25
# resolution=10000
# balance=0

# Run Data Preparation, take K562 cell line as example
python -u data.py --cm "${data_dir}/k562/k562.mcool" --loop "${loop_dir}/ctcf_k562.bedpe" --output "${output_dir}/k562" --prefix "k562_ctcf_hic" --val_chr 1 9 14 --box $box_size





