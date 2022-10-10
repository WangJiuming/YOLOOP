import os
import numpy as np
import json
import argparse
from tqdm.auto import tqdm

from utils.data_utils import load_loops, load_contact_map


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser(description='visualization')
    parser.add_argument('--map', default='', help='path to a contact map with .cool/.mcool extension')
    parser.add_argument('--loop', default='', help='ground truth loop labels')
    parser.add_argument('--chr_num', default='1', help='which chromosome to use')
    parser.add_argument('--window', default=512, help='windows size (unit: resolution)')
    parser.add_argument('--box', default=10, help='bounding box size (unit: resolution)')
    parser.add_argument('--resolution', default=10000, help='resolution')
    args = parser.parse_args()

    # define parameters
    protocols = ['hic', 'chiapet', 'hichip', 'dna-sprite']
    cells = ['gm12878', 'k562', 'hap1', 'h1esc']

    map_path = args.map
    loop_path = args.loop

    resolution = int(args.resolution)
    window_size = int(args.window)
    box_size = args['box_size']

    map_name = map_path.split('/')[-1]
    out_dir = os.path.join(f'data/{map_name}')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    loop_dic = {}

    # iterate through all chromosomes
    for i in tqdm(range(1, 23), position=0, desc='chr1-22', leave=False):
        loops = load_loops(loop_path, str(i))
        loops_candidates = load_loops(loop_path, str(i))

        print('=> loading {}'.format(map_path))
        chr_matrix = load_contact_map(map_path, str(i))
        print(f'matrix for chr{i} loaded, shape: {chr_matrix.shape}')

        chr_start = 0
        chr_end = chr_matrix.shape[0] * resolution

        # create a numpy array to store the windows
        mat_list = []
        # create a list to store the loops
        loop_dic[i] = []

        flag = 0
        for loop in tqdm(loops, position=1, desc='loops', leave=True):
            # find the window centered at the loop
            # x_start and y_start defines the top-left corner
            x_boxleft = loop[1] - int(window_size / 2) * resolution
            y_boxleft = loop[2] - int(window_size / 2) * resolution

            window_mat = chr_matrix[int(x_boxleft / resolution):int(x_boxleft / resolution) + window_size,
                                    int(y_boxleft / resolution):int(y_boxleft / resolution) + window_size]
            if window_mat.shape != (window_size, window_size):
                print("window_mat :{}".format(window_mat.shape))
                print("loop centered at ({},{}) not usable, continue".format(x_boxleft, y_boxleft))
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
                    window_loops.append([int(a_loop[1] / resolution) - int(x_boxleft / resolution),
                                         int(a_loop[2] / resolution) - int(y_boxleft / resolution)])

            loop_dic[i].append(window_loops.copy())
            mat_list.append(np.array(window_mat))

        mat_data = np.zeros((len(mat_list), window_size, window_size))
        for k in range(len(mat_list)):
            mat_data[k, :, :] = mat_list[k]

        # save the loop data and contact map
        # np.save('{}/chr{}-loop.npy'.format(out_dir,i), np.array(loop_data), allow_pickle=True)
        with open('{}/chr{}-loop.json'.format(out_dir, i), 'w') as outfile:
            json.dump(loop_dic, outfile)
        np.save('{}/chr{}-map.npy'.format(out_dir, i), mat_data, allow_pickle=True)

        print('chr{}: contact maps and loops saved'.format(i))

    print('data generation complete')
