import numpy as np
from PIL import Image
import json
import argparse
import os
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
    parser.add_argument('--split', default=18, help='train-val-split, how many chromosomes of 22 autosomes '
                                                    'are used to train the model')
    parser.add_argument('--resolution', default=10000, help='resolution')
    args = parser.parse_args()

    map_path = args.map
    loop_path = args.loop

    train_val_split = int(args.split)
    resolution = int(args.resolution)
    window_size = int(args.window)
    box_size = int(args.box)

    map_name = map_path.train_val_split('/')[-1]
    out_dir = os.path.join(f'data/{map_name}')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    loop_dic = {}

    # iterate through all chromosomes
    for i in tqdm(range(1, 23), position=0, desc='chr1-22', leave=False):
        loops = load_loops(loop_path, str(i))
        loops_candidates = load_loops(loop_path, str(i))

        print(f'=> loading {map_path}')
        chr_matrix = load_contact_map(map_path, str(i))
        print(f'matrix for chr{i} loaded, shape: ({chr_matrix.shape[0]}, {chr_matrix.shjape[1]})')

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
        with open(f'{out_dir}/chr{i}-loop.json', 'w') as outfile:
            json.dump(loop_dic, outfile)
        np.save(f'{out_dir}/chr{i}-map.npy', mat_data, allow_pickle=True)

        print(f'chr{i}: contact maps and loops saved')

    print('begin train-val split and data organization')

    # creating directories for data organization
    img_dir = f'{out_dir}/images'
    label_dir = f'{out_dir}/labels'
    if not os.path.exists(img_dir):
        os.mkdir(img_dir)
    if not os.path.exists(label_dir):
        os.mkdir(label_dir)

    if not os.path.exists(img_dir + '/train'):
        os.mkdir(img_dir + '/train')
    if not os.path.exists(img_dir + '/val'):
        os.mkdir(img_dir + '/val')
        if not os.path.exists(label_dir + '/train'):
            os.mkdir(label_dir + '/train')
            if not os.path.exists(label_dir + '/val'):
                os.mkdir(label_dir + '/val')

    sample_count = 0
    for i in range(1, 23):
        mats = np.load(f'{out_dir}/chr{i}-map.npy', allow_pickle=True)
        loops = json.load(open(f'{out_dir}/chr{i}-loop.json'))[str(i)]
        # print(f'contact map of size ({mat.shape[0]}, {mat.shape[1]}, {mat.shape[2]}) loaded from chr{i}')
        print(f'contact map of size {mats.shape} loaded from chr{i}')
        print(f'{len(loops)} loaded from chr{i}')

        sample_num = mats.shape[0]
        # for every cropped contact map window in this chromosome
        for j in range(sample_num):
            # the sample's contact map
            map_img = Image.fromarray(mats[j])
            map_array = np.array(map_img)
            # the sample's loop list (positive labels)
            img_loops = loops[j]
            loop_num = len(img_loops)

            if i < train_val_split:
                # train set
                split_dir = 'train'
            else:
                # validation set
                split_dir = 'val'

            map_img.convert('L').save(f'{out_dir}/images/{split_dir}/{sample_count:08}.png')
            with open(f'{out_dir}/labels/train/{sample_count:08}.txt', 'w') as file:
                # for every loop in this sample's contact map window
                for k in range(0, loop_num):
                    # class, x center, y center, width, height
                    # N.B. reverse the order of x, y in numpy
                    file.write(f'{0:<4}'  # class 0
                               f'{img_loops[k][1] / window_size:<15}'  # x-coordinate
                               f'{img_loops[k][0] / window_size:<15}'  # y-coordinate
                               f'{box_size / window_size:<15}'  # width of the bounding box
                               f'{box_size / window_size:<15}\n')  # height of the bounding box

            sample_count += 1
        print(f'chr{i} done')

    print('data generation all done')
    print('the data are now ready to be used for training the model')

