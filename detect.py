#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2022/4/17
# @Author   : Wang Jiuming
# @File     : detect.py
# @Software : PyCharm
import numpy as np
from PIL import Image
import argparse
from tqdm.auto import tqdm
import json
import cooler
import os
import torch
import time
import resource


os.environ['CUDA_VISIBLE_DEVICES'] = '0'
if torch.cuda.is_available():
    device = 'cuda'
else:
    device = 'cpu'


def sliding_window(grand_matrix, model, window_size=1024, conf_thresh=0.5):
    """
    scan the contact map with a sliding window and return the predicted loop coordinates
    :param grand_matrix: numpy array, the entire contact map
    :param model: object, trained model for detection
    :param window_size: int, the size of the sliding window, preferably a power of 2, equals the stride
    :param conf_thresh: float, the value for thresholding predictions based on the confidence score, default 0.5
    :return: list, absolute coordinates of the predicted loops, boxes[i] = [vertical coord, horizontal coord, score]
    """
    boxes = []

    # N.B. reverse order of width and height from numpy array
    w, h = grand_matrix.shape[1], grand_matrix.shape[0]
    print(f'w={w}, h={h}, window_size={window_size}')
    diag_limit = window_size * 5

    print('=> sliding through the matrix')
    for y_offset in tqdm(np.arange(0, h - window_size + 1, window_size)):
        for x_offset in np.arange(0, w - window_size + 1, window_size):
            if abs(x_offset - y_offset) > diag_limit:
                continue

            mat = grand_matrix[y_offset:(y_offset + window_size), x_offset:(x_offset + window_size)]

            results = model(mat).tolist()[0].xyxy

            for i in range(results.shape[0]):
                coord = results[i]
                # coord = [bottom_left_x, bottom_left_y, upper_right_x, upper_right_y, confidence, class]
                # thresholding by the prediction confidence score
                if coord[4] < conf_thresh:
                    continue
                # absolute loop coordinates
                loop_x = int((coord[0] + coord[2]) / 2) + x_offset
                loop_y = int((coord[1] + coord[3]) / 2) + y_offset

                # reverse to numpy-order for easy processing later
                boxes.append([loop_y, loop_x, coord[4]])

    print('prediction complete')

    return boxes


def save_bedpe_file(loop_list, chr_num, resolution, file_name):
    """
    write the final prediction results to a .bedpe file
    :param loop_list: list of list, loop_list[i] = [vertical coord, horizontal coord, score] for a predicted loop
    :param chr_num: int, the index number of the chromosome
    :param resolution: int, resolution of the contact map used in the experiment
    :param file_name: string, path to the .bedpe file to be written and saved
    :return: None
    """
    with open(file_name, 'w') as f:
        for loop in loop_list:
            score = loop[2]

            # assume that the loop has the size of 1 resolution
            if loop[0] > loop[1]:
                start_1 = int(loop[1]) * resolution
                end_1 = int(loop[1] + 1) * resolution
                start_2 = int(loop[0]) * resolution
                end_2 = int(loop[0] + 1) * resolution
            else:
                start_1 = int(loop[0]) * resolution
                end_1 = int(loop[0] + 1) * resolution
                start_2 = int(loop[1]) * resolution
                end_2 = int(loop[1] + 1) * resolution

            line = "chr{}\t{}\t{}\tchr{}\t{}\t{}\t{}\n".format(chr_num, start_1, end_1, chr_num, start_2, end_2, score)
            f.write(line)


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser(description='yolop_detection')
    parser.add_argument('--map', default='', help='hic file for loop detection')
    parser.add_argument('--weight', default='', help='model weight path')
    parser.add_argument('--output', default='output', help='output path')
    parser.add_argument('--window', default=512, help='window size')
    parser.add_argument('--resolution', default=10000, help='resolution')
    parser.add_argument('--threshold', default=0.6, help='prediction confidence threshold')
    args = parser.parse_args()

    map_path = args.map
    weight_path = args.weight
    window_size = int(args.window)
    resolution = int(args.resolution)
    threshold = float(args.threshold)

    print(f'using window size: {window_size}, resolution: {resolution}, threshold: {threshold}')

    map_name = map_path.split('/')[-1].split('.')[-2]
    output_path = os.path.join(args.output, f'{map_name}-{resolution}')

    if not os.path.exists(os.path.dirname(output_path)):
        os.mkdir(os.path.dirname(output_path))

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # loading the contact map
    print(f'=> loading contact map:{map_path}')
    if '.mcool' in map_path:
        raw_mat = cooler.Cooler(f'{map_path}::/resolutions/{resolution}')
    elif '.cool' in map_path:
        raw_mat = cooler.Cooler(f'{map_path}')
    else:
        raise Exception('File Format Error: cannot recognize the input file format, please use a .cool or .mcool file')

    grand_matrix = raw_mat.matrix(balance=False)

    # load model
    print(f'=> loading model with customized weights from {weight_path}')
    model = torch.hub.load('ultralytics/yolov3', 'custom', path=weight_path, force_reload=True).to(device)

    # Loop through all chromosomes
    for chr_num in range(1, 23):
        time_start = time.perf_counter()
        print(f'=> scanning chr{chr_num}')

        try:
            contact_matrix = grand_matrix.fetch('{}'.format(chr_num))
        except:
            contact_matrix = grand_matrix.fetch('chr{}'.format(chr_num))

        chr_start = 0
        chr_end = (contact_matrix.shape[0]) * resolution
        if contact_matrix.shape[0] % window_size != 0 or contact_matrix.shape[1] % window_size != 0:
            bottom_pad = window_size - contact_matrix.shape[0] % window_size
            right_pad = window_size - contact_matrix.shape[1] % window_size
            # print(bottom_pad, right_pad)
            contact_matrix = np.pad(contact_matrix, ((0, bottom_pad), (0, right_pad)), 'constant')

        # load model
        compute_start = time.perf_counter()
        contact_matrix = np.array(Image.fromarray(contact_matrix).convert('L'))
        boxes = sliding_window(contact_matrix, model, window_size, threshold)

        # print summary information
        compute_time = (time.perf_counter() - compute_start)
        total_time = (time.perf_counter() - time_start)
        print('Summary:')
        print(f'Number of loop predictions:{len(boxes)}')
        print(f'Time usage: {total_time}')
        speed = (total_time / (chr_end - chr_start)) * 1000000
        mem_use = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
        print(f'Speed: {speed:.5f} secs/1Mbp')
        print(f'Memory usage: {mem_use:.1f} MByte')

        # save the predictions
        bedpe_file_name = os.path.join(output_path, f'{map_name}-chr{chr_num}-pred.bedpe')
        save_bedpe_file(boxes, chr_num, resolution, bedpe_file_name)

        speed_file = os.path.join(output_path, "speed_memory.json")
        if os.path.exists(speed_file):
            data = json.load(open(speed_file))
            data[chr_num] = f'chr{chr_num} Total time: {total_time} sec, ' \
                            f'Compute time: {compute_time}, Speed: {speed} secs/1Mbp, Mem: {mem_use} MByte'
        else:
            data = {chr_num: f'chr{chr_num} Total time: {total_time} sec, '
                             f'Compute time: {compute_time}, Speed: {speed} secs/1Mbp, Mem: {mem_use} MByte'}

        with open(speed_file, "w") as f:
            json.dump(data, f, indent=4)
            print("=> {} written".format(speed_file))

    print('genome-wide detection complete')
