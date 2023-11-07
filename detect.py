import os
import sys
from pathlib import Path
import argparse

import torch.cuda

import numpy as np
import cooler

# add our local modified YOLO project directory to the top of Python paths
sys.path.insert(0, os.path.abspath('./ultralytics/'))

from ultralytics import YOLO


def write_bedpe(loop_list, bedpe_path, resolution):
    """
    write the loop prediction results to a .bedpe file

    Args:
        loop_list (list): list of loop coordinates and confidence scores
        bedpe_path (str): path to the output .bedpe file
        resolution (int): resolution of the contact map
    """
    loop_list = sorted(loop_list, key=lambda x: x[0])
    with open(bedpe_path, 'w') as bedpe_file:
        for loop in loop_list:
            x_start = loop[1] * resolution
            x_end = (loop[1] + 1) * resolution
            y_start = loop[2] * resolution
            y_end = (loop[2] + 1) * resolution
            bedpe_file.write(f'{loop[0]}\t{x_start}\t{x_end}\t{y_start}\t{y_end}\t{loop[3]}\n')


def sliding_window(model, cm, window_size, conf_thresh, restrict=True):
    """
    use a sliding window to detect loops from a chromosome's contact matrix

    Args:
        model (pytorch model object): pretrained YOLO model
        cm (numpy array): contact matrix for a chromosome
        window_size (int): size of the sliding window
        conf_thresh (float): threshold for the detection confidence score
        restrict (bool): whether to detect only around the main diagonal, True by default

    Returns:
        list: list of detected loop coordinates
    """
    preds = []

    # reverse order of width and height from numpy array
    w, h = cm.shape[1], cm.shape[0]
    print(f'[info] width: {w}, height: {h}, window size: {window_size}')
    diag_limit = window_size * 5

    print('[info] scanning the contact matrix ...')

    for y_offset in np.arange(0, h - window_size + 1, window_size):
        for x_offset in np.arange(0, w - window_size + 1, window_size):
            # skip the edges if restricting the detection to the diagonal parts only
            if restrict and abs(x_offset - y_offset) > diag_limit:
                continue

            mat = cm[y_offset:y_offset + window_size, x_offset:x_offset + window_size]

            results = model(mat)[0]

            for result in results:
                coords = result.boxes.xyxy
                confs = result.boxes.conf

                for coord, conf in zip(coords, confs):
                    # coord = [bottom_left_x, bottom_left_y, upper_right_x, upper_right_y, confidence, class]
                    # threshold by the prediction confidence score
                    if conf < conf_thresh:
                        continue

                    # absolute loop coordinates
                    loop_x = int((coord[0] + coord[2]) / 2) + x_offset
                    loop_y = int((coord[1] + coord[3]) / 2) + y_offset

                    # reverse to numpy-order for easy processing later
                    preds.append([loop_y, loop_x, conf.item()])

        return preds


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='YOLOOP for efficient chromatin loop detection')

    parser.add_argument('--cm', type=str, help='path to the input contact matrix with .mcool/.cool extension')
    parser.add_argument('-r', '--resolution', type=int, help='resolution of the contact matrix')
    parser.add_argument('-b', '--balance', action='store_true', help='whether to use a balanced matrix or not')
    parser.add_argument('-m', '--model', type=str, help='YOLOOP model checkpoint to be loaded')
    parser.add_argument('--out', type=str, help='output directory for saving the prediction results')
    parser.add_argument('-t', '--thresh', type=float, default=0.0, help='threshold for the confidence score')
    parser.add_argument('--device', type=str, default='none', help='device to be used, e.g., cuda, cuda:0, cpu')

    args = parser.parse_args()

    # get the device
    if args.device != 'none':
        device = args.device
    elif torch.cuda.is_available():
        device = 'cuda'
    else:
        device = 'cpu'

    # load the pre-trained model
    model_path = args.model
    model = YOLO(model_path).to(device)

    # load the contact map
    cm_path = args.cm
    resolution = args.resolution
    print(f'[info] input path: {cm_path}, resolution: {resolution}')

    if '.mcool' in cm_path:
        cm_obj = cooler.Cooler(f'{cm_path}::/resolutions/{resolution}')
    elif '.cool' in cm_path:
        cm_obj = cooler.Cooler(cm_path)
    else:
        raise Exception(f'file format {Path(cm_path).suffix} not supported')

    # get the chromosomes to make predictions on
    chr_list = cm_obj.chromnames
    if 'MT' in chr_list:
        chr_list.remove('MT')
    print('[info] chromosome list: ', chr_list)

    genome_cm = cm_obj.matrix(balance=args.balance)

    print('[info] begin prediction')

    all_results = []

    for chr_name in chr_list:
        chr_cm = genome_cm.fetch(chr_name)  # get the contact matrix for this chromosome

        results = sliding_window(model, chr_cm, 512, args.thresh)

        print(f'[info] number of predictions for chromosome {chr_name.replace("chr", "")}: {len(results)}')

        results = [[chr_name] + result for result in results]

        all_results += results

    # save the prediction results
    output_dir = args.out
    os.makedirs(output_dir, exist_ok=True)

    write_bedpe(all_results, os.path.join(output_dir, f'yoloop_pred_{Path(cm_path).stem}.bedpe'), resolution)

    print('[info] detection: done')
