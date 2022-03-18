"""
Usage: train <json_file>
"""

import numpy as np
import matplotlib.pyplot as plt
import cooler
import straw
import json

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


def load_loops(loop_path, chromosome, start=None, end=None):
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
    if start is None and end is None:
        for loop in raw_loops:
            if loop[0] != chromosome:
                continue
            loop_start = int((int(loop[1]) + int(loop[2])) / 2)
            loop_end = int((int(loop[4]) + int(loop[5])) / 2)
            loops.append([loop[0], loop_start, loop_end])
            loops.append([loop[0], loop_end, loop_start])
    else:
        for loop in raw_loops:
            if loop[0] != chromosome:
                continue
            # use the middle point to approximate the loop starting/ending point
            # loop_start = x_mid, loop_end = y_mid
            loop_start = int((int(loop[1]) + int(loop[2])) / 2)
            loop_end = int((int(loop[4]) + int(loop[4])) / 2)

            if loop_start < start or loop_end > end:
                continue
            loops.append([loop[0], loop_start, loop_end])
            loops.append([loop[0], loop_end, loop_start])

    print(f'{int(len(loops) / 2)} loops loaded from {chromosome}')
    return loops


def load_a_window(hic_path, chromosome, x_start, x_end, y_start, y_end, resolution=10000):
    """
    load a contact matrix window
    :param hic_path: path, path to a .hic or .mcool file of contact matrix
    :param chromosome: string, the chromosome of interest, e.g., 'chr1' or 'chrX'
    :param x_start: int, x coordinate of the start of the window
    :param x_end: int, x coordinate of the end of the window
    :param y_start: int, y coordinate of the start of the window
    :param y_end: int, y coordinate of the end of the window
    :param resolution: int, resolution of the contact matrix
    :return: numpy array, a block of the contact matrix
    """
    # check for any invalid inputs
    if x_start > x_end or y_start > y_end:
        raise Exception("Invalid arguments: starting position larger than ending position")

    # straw API takes a number in the format of a string as chromosome index
    chr_num = chromosome[3:]
    # determine the file type of the contact matrix
    hic, mcool = False, False
    if '.hic' in hic_path:
        hic = True
    elif '.mcool' in hic_path:
        mcool = True

    m = int((x_end - x_start) / resolution)
    n = int((y_end - y_start) / resolution)
    matrix = np.zeros((m, n))

    # load matrix from a .hic file
    if hic:
        # straw.straw automatically prints 'Hi-C version' to stdout
        # block the message temporarily
        # block_print()
        # data is inclusive of [start, end]
        data = straw.straw('VC', '{}'.format(hic_path), '{}:{}:{}'.format(chr_num, x_start, x_end),
                           '{}:{}:{}'.format(chr_num, y_start, y_end), 'BP', resolution)
        # enable_print()
        # x/ y/ count

        l = len(data[0])
        for i in range(l):
            x = int((data[0][i] - x_start) / resolution)
            y = int((data[1][i] - y_start) / resolution)

            # check for boundary cases
            if x >= m or y >= n:
                continue
            if x < 0 or y < 0:
                continue

            matrix[x, y] = data[2][i]

        # if the window include part of the diagonal
        # then fill in the blank lower left corner with its symmetric part
        # the only case would be x_start <= y_start < x_end
        # since .hic only stores data above the diagonal
        if x_start <= y_start < x_end:
            # sub- marks the starting and ending points of the half-blank lower left window
            sub_x_start = int((y_start - x_start) / resolution)
            sub_x_end = int((x_end - x_start) / resolution)
            # sub_y_start is initialized in the inner loop
            sub_y_end = sub_x_end - sub_x_start
            for j in range(sub_x_start, sub_x_end):
                for k in range(j - sub_x_start + 1, sub_y_end):
                    # find the symmetric point of (j, k) about y = x - sub_x_start
                    matrix[k + sub_x_start, j - sub_x_start] = matrix[j, k]

    # load matrix from a .mcool file
    elif mcool:
        data = cooler.Cooler('{}::/resolutions/{}'.format(hic_path, resolution))
        try:
            matrix = data.matrix(balance=True).fetch('{}:{}-{}'.format(chromosome, x_start, x_end),
                                                     '{}:{}-{}'.format(chromosome, y_start, y_end))
        # handle the balance failed exception
        # write a normalize function to resolve this later
        except:
            print('Normalization failed, proceeding without normalization')
            matrix = data.matrix(balance=False).fetch('{}:{}-{}'.format(chromosome, x_start, x_end),
                                                      '{}:{}-{}'.format(chromosome, y_start, y_end))

    return matrix


if __name__ == '__main__':
    # load arguments from command line input
    args = docopt(__doc__)
    args = json.load(open(args['<json_file>']))

    # define parameters
    protocols = ['hic', 'chiapet', 'hichip', 'dna-sprite']
    cells = ['gm12878', 'k562', 'hap1', 'h1esc']

    resolution = args['resolution']
    window_size = args['window_size']
    box_size = args['box_size']

    hic_path = args['hic_path']
    loop_path = args['loop_path']

    flag = 0  # indicate whether the window is a special one, i.e., top-left corner

    # iterate through all chromosomes
    for i in range(1, 23):
        loops = load_loops(loop_path, f'chr{i}')
        a_loops = load_loops(loop_path, f'chr{i}')
        x_max, y_max = np.max(np.array(loops)[:, 1:].astype(int), axis=0)

        chr_matrix = load_a_window(hic_path, f'chr{i}', 0, x_max + window_size * resolution,
                                   0, y_max + window_size * resolution)
        print(f'matrix for chr{i} loaded, shape: {chr_matrix.shape}')

        # create a numpy array to store the windows
        mat_data = np.zeros((len(loops), window_size, window_size))
        # create a list to store the loops
        loop_data = []

        # iterate through all loops on this chromosome
        for j, loop in enumerate(loops):
            # find the window centered at the loop
            # x_start and y_start defines the top-left corner
            x_start = loop[1] - int(window_size / 2) * resolution
            y_start = loop[2] - int(window_size / 2) * resolution

            # write the contact map window to our numpy array
            try:
                np.copyto(mat_data[j], chr_matrix[int(x_start / resolution):int(x_start / resolution) + window_size,
                                       int(y_start / resolution):int(y_start / resolution) + window_size])
            except ValueError:
                # the loop position is smaller than window_size*resolution
                np.copyto(mat_data[j], chr_matrix[:window_size, :window_size])
                flag = 1

            # identify the loops in the window
            window_loops = []
            for a_loop in a_loops:
                if flag == 1:
                    if box_size * resolution <= a_loop[1] <= (window_size - box_size) * resolution and \
                            box_size * resolution <= a_loop[2] <= (window_size - box_size) * resolution:
                        # relative loop center coordinate
                        window_loops.append([int(a_loop[1] / resolution), int(a_loop[2] / resolution)])
                    flag = 0
                elif x_start + box_size * resolution <= a_loop[1] <= x_start + (window_size - box_size) * resolution and \
                        y_start + box_size * resolution <= a_loop[2] <= y_start + (window_size - box_size) * resolution:
                    # relative loop center coordinate
                    window_loops.append(
                        [int(a_loop[1] / resolution) - int(x_start / resolution),
                         int(a_loop[2] / resolution) - int(y_start / resolution)])

            loop_data.append(np.array(window_loops).copy())

        # save the loop data and contact map
        np.save(f'./data/chr{i}-loop.npy', np.array(loop_data), allow_pickle=True)
        np.save(f'./data/chr{i}-map.npy', mat_data[:len(loop_data)], allow_pickle=True)

        print(f'chr{i}: contact maps and loops saved')
        print(f'map shape: {mat_data[:len(loop_data)].shape}')
        print(f'loop number: {len(loop_data)}')
        print('-----------------------------------')

    print('data generation complete')
