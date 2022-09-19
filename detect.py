#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2022/4/17
# @Author   : Wang Jiuming
# @File     : detect.py
# @Software : PyCharm
import json as js
import numpy as np
from PIL import Image
import argparse
from tqdm import tqdm
import cooler
import straw
import os
import struct
import torch
import time
import resource
import multiprocessing
import matplotlib.pyplot as plt
from sklearn.preprocessing import Normalizer
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
# root = os.path.dirname(os.path.abspath(__file__))


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

    print('=> Sliding through the matrix')
    for y_offset in tqdm(np.arange(0, h - window_size + 1, window_size)):
        for x_offset in np.arange(0, w - window_size + 1, window_size):

            hic_mat = grand_matrix[y_offset:(y_offset + window_size), x_offset:(x_offset + window_size)]

            results = model(hic_mat).tolist()[0].xyxy

            for i in range(results.shape[0]):
                coord = results[i]
                # print(coord)
                # print(coord.shape, type(coord))
                # coord = [bottom_left_x, bottom_left_y, upper_right_x, upper_right_y, confidence, class]
                # thresholding by the prediction confidence score
                if coord[4] < conf_thresh:
                    continue
                # absolute loop coordinates
                loop_x = int((coord[0] + coord[2]) / 2) + x_offset
                loop_y = int((coord[1] + coord[3]) / 2) + y_offset

                # reverse to numpy-order for easy processing later
                boxes.append([loop_y, loop_x, coord[4]])

    print('Prediction complete.')

    return boxes


def save_bedpe_file(loop_list, chr_num, resolution, file_name):
    """
    write the final prediciton results to a .bedpe file
    :param loop_list: list of list, loop_list[i] = [vertical coord, horizontal coord, score] for a predicted loop
    :param chr_num: int, the index number of the chromosome
    :param resolution: int, resolution of the contact map used in the experiment
    :param file_name: string, path to the .bedpe file to be written and saved
    :return: None
    """
    with open(file_name, 'w') as f:
        for loop in loop_list:
            score = loop[2]

            if loop[0] > loop[1]:
                start_1 = int(loop[1]) * resolution
                end_1 =  int(loop[1] + 1) * resolution
                start_2 =  int(loop[0]) * resolution
                end_2 = int(loop[0] + 1) * resolution
            else:
                start_1 =  int(loop[0]) * resolution
                end_1 =  int(loop[0] + 1) * resolution
                start_2 =  int(loop[1]) * resolution
                end_2 =  int(loop[1] + 1) * resolution

            line = "chr{}\t{}\t{}\tchr{}\t{}\t{}\t{}\n".format(chr_num, start_1, end_1, chr_num, start_2, end_2, score)
            f.write(line)


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser(description='yolop_detection')
    parser.add_argument('--device', default='cpu', help='default: cpu')
    parser.add_argument('--hic_file', default='', help='hic file for loop detection')
    parser.add_argument('--weight', default='', help='model weight path')
    parser.add_argument('--output', default='results_yolo', help='output path')

    parser.add_argument('--resolution', default=10000, help='resolution')

    parser.add_argument('--window', default=512, help='window size')
    parser.add_argument('--threshold', default=0.6, help='prediction confidence threshold')

    args = parser.parse_args()
    device = args.device
    if torch.cuda.is_available():
        device = 'cuda'
    else:
        device = 'cpu'
    hic_path = args.hic_file

    weight_path = args.weight
    window_size = int(args.window)
    resolution = int(args.resolution)
    threshold = float(args.threshold)

    print(f'using window size: {window_size}, {type(window_size)}, resolution: {resolution}, threshold: {threshold}')

    info_list = hic_path.split('/')
    info = info_list[-4] + '_' + info_list[-3] + '_' + info_list[-2]
    output_path = os.path.join(args.output, info + "_" + str(window_size))

    if not os.path.exists(os.path.dirname(output_path)):
        os.mkdir(os.path.dirname(output_path))

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # loading mcool data once
    print('=> Loading Hi-C file :{}'.format(hic_path))
    raw_cool_data = cooler.Cooler('{}::/resolutions/{}'.format(hic_path,resolution))
    grand_matrix = raw_cool_data.matrix(balance=False)
    # Load model
    print("=> Loading module from {}".format(weight_path))
    model = torch.hub.load('ultralytics/yolov3', 'custom', path=weight_path, force_reload=True).to(device)

    # Loop through all chromosomes
    for chrom in range(1,23):
        time_start = time.perf_counter()
        print('=> Doing chromosome {}'.format(chrom))

        try:
            contact_matrix = grand_matrix.fetch('{}'.format(chrom))
        except:
            contact_matrix = grand_matrix.fetch('chr{}'.format(chrom))

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

        compute_time = (time.perf_counter() - compute_start)
        total_time = (time.perf_counter() - time_start)
        # visualize(contact_matrix, boxes)
        print('Summary:')
        print('Loop predictions:{}'.format(len(boxes)))
        print('Time usage: {}'.format(total_time))
        memMb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
        print("Time usage: %5.3f secs/1Mbp %5.1f MByte" % ((total_time / (chr_end - chr_start)) * 1000000, memMb))
        bedpe_file_name = os.path.join(output_path, info + '_{}_chr{}_pred.bedpe'.format('yolov3',
                                                                                        str(chrom)))
        save_bedpe_file(boxes, chrom, resolution, bedpe_file_name)

        speed_file = os.path.join(output_path, "speed_memory.json")
        if os.path.exists(speed_file):
            data = js.load(open(speed_file))
            data[chrom] = "Chr {} Total Time: {} sec, Compute Time:{}, Speed:{} secs/1Mbp Mem:{} MByte".format(chrom, total_time,compute_time,(total_time / (chr_end - chr_start)) * 1000000, memMb)
        else:
            data = {}
            data[chrom] = "Chr {} Total Time: {} sec, Compute Time:{}, Speed:{} secs/1Mbp Mem:{} MByte".format(chrom, total_time,compute_time,(total_time / (chr_end - chr_start)) * 1000000, memMb)

        with open(speed_file, "w") as f:
            js.dump(data, f,indent=4)
            print("=> {} wirtten".format(speed_file))
