import os
import sys
import time
import argparse
from tqdm import tqdm
import json
# add our local modified YOLO project directory to the top of Python paths
sys.path.insert(0, os.path.abspath('./ultralytics/'))
import resource
from ultralytics import YOLO
# from PIL import Image
import numpy as np
import cooler

from util_data import write_bedpe


def sliding_window(model, cm, window_size, conf_thresh, restrict=True):
    """
    use a sliding window to detect loops from a chromosome's contact matrix

    Args:
        model (object): pretrained YOLO model
        cm (array): contact matrix for a chromosome
        window_size (int): size of the sliding window
        conf_thresh (float): threshold for the detection confidence score
        restrict (bool): whether to detect only around the main diagonal, True by default

    Returns:
        list: list of detected loop coordinates
    """
    preds = []

    # reverse order of width and height from numpy array
    w, h = cm.shape[1], cm.shape[0]
    print(f'[debug] width: {w}, height: {h}, window size: {window_size}')
    diag_limit = window_size * 5

    print('[info] sliding through the contact matrix')
    for y_offset in tqdm(np.arange(0, h - window_size + 1, window_size)):
        for x_offset in np.arange(0, w - window_size + 1, window_size):
            # skip the edges if restricting the detection to the diagonal parts only
            if restrict and abs(x_offset - y_offset) > diag_limit:
                continue

            mat = cm[y_offset:y_offset + window_size, x_offset:x_offset + window_size]

            results = model(mat)[0]

            #print(results.boxes.data)
            #input('proceed? ')
            for result in results:
                coords = result.boxes.xyxy
                confs = result.boxes.conf

                for coord, conf in zip(coords, confs):
                    # coord = [bottom_left_x, bottom_left_y, upper_right_x, upper_right_y, confidence, class]
                    # threshold by the prediction confidence score
            #        if conf < conf_thresh:
            #            continue
                    # absolute loop coordinates
                    loop_x = int((coord[0] + coord[2]) / 2) + x_offset
                    loop_y = int((coord[1] + coord[3]) / 2) + y_offset

                    # reverse to numpy-order for easy processing later
                    preds.append([loop_y, loop_x, conf.item()])

    return preds


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='YOLOOP for chromatin loop detection')

    parser.add_argument('--cm', type=str, help='path to the input contact matrix')
    parser.add_argument('--r', type=float, help='resolution of the contact matrix')
    parser.add_argument('--model', type=str, help='YOLOOP model checkpoint to be loaded')
    parser.add_argument('--out', type=str, help='output directory for saving the prediction results')

    args = parser.parse_args()

    # load the pre-trained model
    model_path = args.model
    model = YOLO(model_path)

    # load the contact map
    cm_path = args.cm
    print('[debug] cm path: ', cm_path)

    cell_type = cm_path.split('/')[4]
    cm_obj = cooler.Cooler(cm_path)

    resolution = args.r
    output_dir = args.out

    chr_list = cm_obj.chromnames
    if 'MT' in chr_list:
        chr_list.remove('MT')
    print('[debug] chromosome list: ', chr_list)

    genome_cm = cm_obj.matrix(balance=False)

    print('[info] begin prediction')

    all_results = []

    # chr_list = ['chr1','chr9','chr14']
    print(chr_list)
    time_json_file = os.path.join(output_dir, 'time_yoloop_{}.json'.format(cell_type))
    time_dic = {}

    for chr_name in chr_list:
        start_time = time.time()

        chr_cm = genome_cm.fetch(chr_name)  # get the contact matrix for this chromosome

        # time the sliding window computation

        results = sliding_window(model, chr_cm, 512, 0.2)

        print(len(results))

        print(f'[info] time for {chr_name}: {time.time() - start_time} seconds')
        print(f'[info] number of predictions: {len(results)}')
        memMb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
        time_dic[chr_name] = [time.time() - start_time,memMb]
        # next, write the prediction results to the output file
        results =[[chr_name] + result for result in results]

        all_results += results


    write_bedpe(all_results,
                os.path.join(output_dir, 'yoloop_pred_{}.bedpe'.format(cell_type)),
                int(resolution))

    print('[info] detection: done')
    with open(time_json_file, "w") as f:
        json.dump(time_dic, f, indent=4)


