"""
@File   : data.py
@Author : Siyuan Chen
@Date   : 2025/3/10
@Desc   :
"""

import json
import os
import time

import argparse
from tqdm import tqdm

import numpy as np
from PIL import Image

from util_data import load_contact_map
from util_data import load_loops


def generate_chr_training_data(cm_path, loop_path, dataset_dir, prefix, chr_name,
                               window_size, box_size, resolution=10000, balance=False, is_val=False):
    """
    create images (contact map windows) and labels (loops) as training data in YOLO format

    Args:
        cm_path (str): path to the file storing the contact matrix, one of .hic/.cool/.mcool formats
        loop_path (str): path to a .bedpe file storing the loop coordinates
        dataset_dir (str): path to the directory of the dataset, to contain the images/ and labels/
        prefix (str): prefix of the dataset to be appended to each sample
        chr_name (str): name of the chromosome of interest, e.g., 'chr1'
        window_size (int): size of the contact map window to be cropped
        box_size (int): size of the bounding box for each loop
        resolution (int): resolution of the contact map to be used, 10000 by default
        balance (bool): whether to use a balanced matrix, False by default
        is_val (bool): whether this chromosome will be used for validation, False by default

    Returns:
        int: number of training samples saved
    """

    # initialize the dataset's file structure
    if is_val:
        img_dir = os.path.join(dataset_dir, 'images', 'val')
        label_dir = os.path.join(dataset_dir, 'labels', 'val')
    else:
        img_dir = os.path.join(dataset_dir, 'images', 'train')
        label_dir = os.path.join(dataset_dir, 'labels', 'train')
    print(f'[info] {prefix}_{chr_name}: image output directory: {img_dir}')
    print(f'[info] {prefix}_{chr_name}: label output directory: {label_dir}')

    os.makedirs(img_dir, exist_ok=True)
    os.makedirs(label_dir, exist_ok=True)

    half_box_size = box_size // 2

    # load the loop coordinates
    loops = load_loops(loop_path, chr_name)

    # load the contact matrix
    start_time = time.perf_counter()
    cm = load_contact_map(cm_path, chr_name, resolution=resolution, balance=balance)  # unbalanced matrix by default
    print(f'[debug] {prefix}_{chr_name}: contact map loaded in {time.perf_counter() - start_time:.6f} seconds')

    x_max = cm.shape[0]
    y_max = cm.shape[1]
    print(f'[info] {prefix}_{chr_name}: contact map has been loaded, shape: {cm.shape}')

    # iterate over all the loops to retrieve the windows centered at each loop and
    # record the loops associated with each window as positive data
    # save the windows as images and associated loops as labels

    data_idx = 0

    skip_num = 0

    for loop in loops:
        # start_time = time.perf_counter()

        # get the contact map window centered at the loop
        # compute the top-left corner coordinates
        # divided by the resolution to get the matrix indices
        # use min and max to handle the out-of-boundary cases
        row_top_left_idx = min(max(loop[1] // resolution - window_size // 2, 0),
                               x_max - 1 - window_size)
        col_top_left_idx = min(max(loop[2] // resolution - window_size // 2, 0),
                               y_max - 1 - window_size)

        # alternatively, load a small chunk of contact map for each loop, which takes <0.15s/iter on Apple M1
        # compared with slicing from the entire contact map, which takes <0.02s/iter
        # x_start = row_top_left_idx * resolution
        # x_end = (row_top_left_idx + window_size) * resolution
        # y_start = col_top_left_idx * resolution
        # y_end = (col_top_left_idx + window_size) * resolution
        # cm_window = load_contact_map(cm_path, chr_name, False, x_start, x_end, y_start, y_end, resolution)

        cm_window = cm[row_top_left_idx:row_top_left_idx + window_size,
                    col_top_left_idx:col_top_left_idx + window_size]

        if cm_window.shape != (window_size, window_size):
            print(f'[warning] loop ({loop[1]:,}, {loop[2]:,}) '
                  f'exceeds the size of the contact map, skipping this annotation')
            skip_num += 1
            continue

        file_name = f'{prefix}_{chr_name}_{data_idx:07}'

        # save the contact map window as a gray-scale image, use .png for lossless compression
        # the axes are interpreted in the reverse order after converting to image objects
        # this will only affect the indexing order when accessing the data
        cm_img = Image.fromarray(cm_window)
        cm_img.convert('L').save(os.path.join(img_dir, f'{file_name}.png'))

        # find all the loops that occurred within this window
        loop_num = len(loops)

        labels = []  # loops within this window

        for i in range(loop_num):
            # convert loop coordinates into indices in the contact matrix
            x_idx = loops[i][1] // resolution
            y_idx = loops[i][2] // resolution

            delta_x = x_idx - row_top_left_idx
            delta_y = y_idx - col_top_left_idx

            if half_box_size <= delta_x <= window_size - half_box_size and \
                    half_box_size <= delta_y <= window_size - half_box_size:
                # record the normalized loop center indices relative to this window
                # need to swap the indices to match with the YOLO index order:
                # column index, followed by row index
                # convert to format string at this step for removing duplicates next
                labels.append(
                    (f'{delta_y / window_size:.10f}', f'{delta_x / window_size:.10f}'))  # normalized by window size

        # save the loops as labels
        unique_labels = set(labels)  # remove duplicated labels due to precision and resolution, etc.

        box_size_normalized = f'{box_size / window_size:.10f}'

        with open(os.path.join(label_dir, f'{file_name}.txt'), 'w') as label_file:
            for label in unique_labels:
                # the width and height of the bounding box always equal the box size
                # YOLO label: column center, row center, width, height (everything should be normalized by image size)
                label_file.write(
                    f'{0:<5}{label[0]:<15}{label[1]:<15}{box_size_normalized:<15}{box_size_normalized:<15}\n')

        data_idx += 1

        if data_idx % 1000 == 0:
            print(f'[info] {prefix}_{chr_name}: {data_idx}/{len(loops)}')

        # print(f'[debug] one data point generated in {time.perf_counter() - start_time:.6f} seconds')

    print(f'[info] {prefix}_{chr_name}: {len(loops)} loops loaded, {skip_num} skipped, {data_idx} saved')
    print(f'[info] {prefix}_{chr_name}: data generation completed for {chr_name}')

    return data_idx


if __name__ == '__main__':

    # config_path = 'config_data.json'

    # load configuration from command line inputs
    parser = argparse.ArgumentParser(description='training data generation for YOLOOP')

    parser.add_argument('--cm', type=str, help='path to the contact matrix file')
    parser.add_argument('--loop', type=str, help='path to the loop annotation file')
    parser.add_argument('--output', type=str, help='path to the output directories for saving the dataset')
    parser.add_argument('--prefix', type=str, help='prefix of the sample names in the dataset')
    parser.add_argument('--window', type=int, help='window size of each sample', default=512)
    parser.add_argument('--box', type=int, help='bounding box size of each loop annotation', default=10)
    parser.add_argument('--resolution', type=int, help='resolution of the contact matrix', default=10000)
    parser.add_argument('--balance', type=bool, help='whether to use a balanced contact matrix, False by default',
                        default=False)
    parser.add_argument('--val_chr', type=str, nargs='*',
                        help='list of chromosomes used for validation, none by default', default=[])

    args = parser.parse_args()
    cm_path = args.cm
    loop_path = args.loop
    dataset_dir = args.output

    prefix = args.prefix

    window_size = args.window
    box_size = args.box
    resolution = args.resolution
    balance = args.balance

    val_chr_list = args.val_chr  # ['1', '9', '14']

    # load the configurations from a .json file
    # with open(config_path, 'r') as config_file:
    #    config = json.load(config_file)

    # cm_path = config['paths']['contact map']
    # loop_path = config['paths']['loop']
    # dataset_dir = config['paths']['output']

    # prefix = config['paths']['prefix']

    # window_size = config['params']['window size']
    # box_size = config['params']['box size']
    # resolution = config['params']['resolution']
    # balance = config['params']['balance']

    print('-------------------- configurations -------------------- ')
    print(f'contact map path: {cm_path}')
    print(f'loop annotation path: {loop_path}')
    print(f'directory for saving data: {dataset_dir}')
    print(f'prefix for saving data: {prefix}')
    print(f'window size: {window_size}')
    print(f'bounding box size: {box_size}')
    print(f'resolution: {resolution}')
    print(f'using balanced contact matrix: {bool(balance)}')
    print(f'chromosomes used for validation:', val_chr_list)
    print('-------------------------------------------------------- ')

    total_data_count = 0

    # generate training data for each chromosome (1 to 22 autosomes)
    if 'mesc' not in cm_path:
        chr_list = ['X', 'Y'] + [f'{i}' for i in range(1, 23)]  # put X, Y in the front
    else:
        chr_list = ['X', 'Y'] + [f'{i}' for i in range(1, 20)]

    print(f'[info] {prefix} generating data from: ')
    print(chr_list)

    print('-' * 50)

    for chr_num in chr_list:
        # perform train-val split
        if chr_num in val_chr_list:
            is_val = True
        else:
            is_val = False

        data_count = generate_chr_training_data(cm_path, loop_path, dataset_dir, prefix, f'chr{chr_num}',
                                                window_size, box_size, resolution, balance, is_val)
        total_data_count += data_count
        print('-' * 50)

    print(f'[info] {prefix}: {total_data_count} data generated and saved in total')
    print('[info] {prefix}: data generation complete')

