#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2022/4/3
# @Author   : Wang Jiuming
# @File     : dataloader.py
# @Software : PyCharm

"""
Usage: train <json_file>
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import cooler
import straw
import json
from tqdm import tqdm
from sklearn.preprocessing import Normalizer

from docopt import docopt


def show(mat):
    """visualize a contact map"""
    figure = plt.figure(figsize=(10, 10))
    axes = figure.add_subplot(111)
    caxes = axes.matshow(np.log10(mat), cmap='Reds')  # OrRd
    plt.show()


def read_loops(loop_path):
    """
    read the loop annotations from a .bedpe file
    :param loop_path: path, path to a loop annotation file
    :return: list of list, list of loop positions
    """
    loops = []
    with open(loop_path, 'r') as loop_file:
        for line in loop_file:
            loops.append(line.split('\t')[:6])

            # loops[i] = ['chr1', 'x_start', 'x_end', 'chr1', 'y_start', 'y_end']
            # Note that all loop[i][j] are strings instead of integers

    return loops


def load_loops(loop_path, chromosome):
    """
    select a list of loop positions from a chromosome
    :param loop_path: path, path to a loop annotation file
    :param chromosome: string, name of the chromosomes of interest, e.g., 'chr1'
    :param start: int, mark the start of the range to find loops, None by default
    :param end: int, mark the end of the range to find loop, None by default
    :return: list of list, list of loop positions
    """
    raw_loops = read_loops(loop_path)
    loops = []
    # select the loops within the range specified by the input

    if 'odc' in loop_path:
        chromosome = 'chr' + chromosome

    for loop in raw_loops:
        if loop[0] != chromosome or loop[1] == 'x1':
            continue

        x = int((int(loop[1]) + int(loop[2])) *0.5)
        y = int((int(loop[4]) + int(loop[5])) *0.5)

        loops.append([loop[0], x, y])
        loops.append([loop[0], y, x])

    print(f'{int(len(loops))} loops loaded from {chromosome}')
    return loops


def load_hic_mat(hic_path, chromosome, resolution=10000):

    # straw API takes a number in the format of a string as chromosome index
    chr_num = chromosome[3:]
    c = cooler.Cooler('{}::/resolutions/{}'.format(hic_path, resolution))
    matrix = c.matrix(balance=False).fetch('{}'.format(chr_num))
    return matrix


if __name__ == '__main__':
    # load arguments from command line input
    args = docopt(__doc__)
    args = json.load(open(args['<json_file>']))

    # define parameters
    protocols = ['hic', 'chiapet', 'hichip', 'dna-sprite']
    cells = ['gm12878', 'k562', 'hap1', 'h1esc']

    resolution = args['resolution']
    window_size = int(args['window_size'])
    box_size = args['box_size']

    hic_path = args['hic_path']
    loop_path = args['loop_path']

    protocol = hic_path.split('/')[-3]
    cell_name = hic_path.split('/')[-2]

    out_dir = os.path.join('data/{}-{}'.format(protocol, cell_name))
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    loop_dic = {}

    # iterate through all chromosomes
    for i in range(1, 23):
        loops = load_loops(loop_path, str(i))
        loops_candidates = load_loops(loop_path, str(i))

        print('=> Loading {}'.format(hic_path))
        chr_matrix = load_hic_mat(hic_path,
                                   f'chr{i}')
        print(f'matrix for chr{i} loaded, shape: {chr_matrix.shape}')

        chr_start = 0
        chr_end = (chr_matrix.shape[0]) * resolution

        # create a numpy array to store the windows
        mat_list = []
        # create a list to store the loops
        loop_dic[i] = []


        flag = 0
        for loop in (loops):
            # find the window centered at the loop
            # x_start and y_start defines the top-left corner
            x_boxleft = loop[1] - int(window_size / 2) * resolution
            y_boxleft = loop[2] - int(window_size / 2) * resolution

            window_mat = chr_matrix[int(x_boxleft / resolution):int(x_boxleft / resolution) + window_size,
                         int(y_boxleft / resolution):int(y_boxleft / resolution) + window_size]
            if window_mat.shape != (window_size, window_size):
                print("window_mat :{}".format(window_mat.shape))
                print("loop centered ({},{}) not suited , continue".format(x_boxleft, y_boxleft))
                continue

            # identify the loops in the window
            window_loops = []
            for a_loop in loops_candidates:
                if flag == 1:
                    if box_size * resolution <= a_loop[1] <= (window_size - box_size) * resolution and \
                            box_size * resolution <= a_loop[2] <= (window_size - box_size) * resolution:
                        # relative loop center coordinate
                        window_loops.append([int(a_loop[1] / resolution), int(a_loop[2] / resolution)])
                    flag = 0
                elif x_boxleft + box_size * resolution <= a_loop[1] <= x_boxleft + (window_size - box_size) * resolution and \
                        y_boxleft + box_size * resolution <= a_loop[2] <= y_boxleft + (window_size - box_size) * resolution:
                    # relative loop center coordinate
                    window_loops.append(
                        [int(a_loop[1] / resolution) - int(x_boxleft / resolution),
                         int(a_loop[2] / resolution) - int(y_boxleft / resolution)])

            loop_dic[i].append(window_loops.copy())

            mat_list.append(np.array(window_mat))

        mat_data = np.zeros((len(mat_list), window_size, window_size))
        for k in range(len(mat_list)):
            mat_data[k,:,:] = mat_list[k]

        # save the loop data and contact map
        # np.save('{}/chr{}-loop.npy'.format(out_dir,i), np.array(loop_data), allow_pickle=True)
        with open('{}/chr{}-loop.json'.format(out_dir,i), 'w') as outfile:
            json.dump(loop_dic, outfile)
        np.save('{}/chr{}-map.npy'.format(out_dir,i), mat_data, allow_pickle=True)

        print('chr{}: contact maps and loops saved'.format(i))

        print('-----------------------------------')

    print('data generation complete')
