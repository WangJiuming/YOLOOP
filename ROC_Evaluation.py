import numpy as np
from tqdm import tqdm
import os
import struct
import argparse
import csv
import json as js


gt_path = {}

# gt_path['gm12878'] read_gt_loops= '/Dataset/HiC/hic/gm12878/hic-gm12878-loop.txt'
gt_path['gm12878'] = '/Dataset/HiC/hichip/gm12878/hichip-gm12878-loop.bedpe'
gt_path['k562'] = '/Dataset/HiC/hic/k562/hic-k562-loop.txt'
gt_path['imr90'] = '/Dataset/HiC/hic/imr90/hic-imr90-loop.txt'
gt_path['kbm7'] = '/Dataset/HiC/hic/kbm7/hic-kbm7-loop.txt'
gt_path['huvec'] = '/Dataset/HiC/hic/huvec/hic-huvec-loop.txt'
gt_path['odc'] = '/Dataset/HiC/hic/odc/hic-odc-loop.bedpe'


def read_gt_loops(loop_path):
    loops = []
    with open(loop_path, 'r') as loop_file:
        for line in loop_file:
            loops.append(line.split('\t')[:6])
    return loops


def read_pred_loops(loop_path):
    loops = []
    with open(loop_path, 'r') as loop_file:
        for line in loop_file:
            loops.append(line.strip('\n').split('\t'))

    return loops


def get_recall_precision(gt_centers, pred_centers,resolution):
    Total_Positive = len(gt_centers)
    TP = 0

    for gt_center in tqdm(gt_centers):

        d_min = 10 ** 10

        prediction = None

        for i in range(len(pred_centers)):
            d = np.linalg.norm(np.array(gt_center) - np.array(pred_centers[i]))

            if abs(d) <= d_min:
                d_min = d
                prediction = pred_centers[i]

        if d_min <= resolution:
            TP += 1

    recall = TP / Total_Positive
    if len(pred_centers) == 0:
        precision = 0
    else:
        precision = TP / len(pred_centers)

    if precision > 1:
        precision = 1

    print("TP:{} Total Loop :{} Prediction:{}".format(TP, Total_Positive, len(pred_centers)))

    return recall, precision



def load_grount_truth(cell_type,args):


    gt_loops = read_gt_loops(gt_path[cell_type])
    print(f'number of all loops = {len(gt_loops)}')
    gt_centers_dict = {}



    for chrom in range(int(args.chr_start), int(args.chr_end)):
        print('chr' + str(chrom))
        gt_centers_dict[str(chrom)] = []


        for gt in gt_loops:

            if gt[1] == 'x1':
                continue

            if gt[0] != 'chr' + str(chrom) :
                continue

            x = int((int(gt[1]) + int(gt[2])) * 0.5)
            y = int((int(gt[4]) + int(gt[5])) * 0.5)

            gt_centers_dict[str(chrom)].append([x, y])
            # gt_centers_dict[str(chrom)].append([y, x])

        print(len(gt_centers_dict[str(chrom)]))
    return gt_centers_dict


def load_HiC_LDNet(args):
    HiC_LDNet_center_dic = {}
    HiC_LDNet_scores_dic = {}

    for chrom in tqdm(range(args.chr_start, args.chr_end)):

        pred_file = 'results_yolo/HiC_' + args.experiment + '_' + args.cell_type + "_" + str(args.window_size)+ \
                    '/' + 'HiC_' + args.experiment + '_' + args.cell_type + '_yolov3_chr{}_pred.bedpe'.format(chrom)


        # the chromosome information already incorporated in the file path
        pred_loops = read_pred_loops(pred_file)
        print(f'total number of loops predicted by our model = {len(pred_loops)}')

        HiC_LDNet_center_dic[str(chrom)] = []
        HiC_LDNet_scores_dic[str(chrom)] = []

        for pred in pred_loops:
            # threshold the predicted loops by the prediction confidence
            if float(pred[6]) < float(args.threshold):
                continue

            HiC_LDNet_center_dic[str(chrom)].append(
                [0.5 * (int(pred[1]) + int(pred[2])), 0.5 * (int(pred[4]) + int(pred[5]))])
            HiC_LDNet_scores_dic[str(chrom)].append(float(pred[6]))

    for key in HiC_LDNet_center_dic.keys():
        print('chr_' + key, "prediction loops", len(HiC_LDNet_center_dic[key]))

    return HiC_LDNet_center_dic

def load_Peakachu(args):
    peakachu_pred_scores = {}
    peakachu_centers = {}

    for chrom in tqdm(range(int(args.chr_start), int(args.chr_end))):

        peakachu_path = '/home/chens0a/HiC-Yolo/baseline/' + '{}-{}/'.format(args.experiment,
                                                                                             args.cell_type)\
                        + args.experiment + '-' + args.cell_type + f'-peakachu/chr{chrom}.bed'
        print("Reading {}".format(peakachu_path))
        peakachu_loops = read_pred_loops(peakachu_path)
        print(f'total number of loops predicted by peakachu = {len(peakachu_loops)}')


        peakachu_pred_scores[str(chrom)] = []
        peakachu_centers[str(chrom)] = []

        for pred in peakachu_loops:

            # pred[0] is in format of 'chr1', chrom is 1
            if "chr" + str(chrom) != pred[0]:
                continue

            # threshold the predictions by the prediction confidence
            if float(pred[6]) < float(args.threshold):
                continue

            x_center = int((int(pred[1]) + int(pred[2])) / 2)
            y_center = int((int(pred[4]) + int(pred[5])) / 2)

            peakachu_centers[str(chrom)].append([x_center, y_center])
            peakachu_pred_scores[str(chrom)].append(float(pred[6]))

    for key in peakachu_centers.keys():
        print('chr_' + key, "prediction loops", len(peakachu_centers[key]))

    return peakachu_centers

def load_chromosight(args):
    chromosight_path = '/Dataset/HiC/' + '{}/{}/'.format(args.experiment,args.cell_type) \
                       + 'default_8_thread.tsv'

    chromosight_loops = []

    with open(chromosight_path) as f:
        lines = csv.reader(f, delimiter='\t')
        for line in lines:
            chromosight_loops.append(line)

    print(f'total number of loops predicted by chromosight = {len(chromosight_loops)}')

    chromosight_centers = {}
    chromesight_pred_scores = {}

    for chrom in tqdm(range(args.chr_start, args.chr_end)):

        chromosight_centers[str(chrom)] = []
        chromesight_pred_scores[str(chrom)] = []

        for pred in chromosight_loops:

            # pred[0] is in format of '1', chrom is 1
            if pred[0] != str(chrom):
                continue

            if float(pred[10]) < float(args.threshold):
                continue

            x_center = int((int(pred[1]) + int(pred[2])) / 2)
            y_center = int((int(pred[4]) + int(pred[5])) / 2)

            chromosight_centers[str(chrom)].append([x_center, y_center])
            chromesight_pred_scores[str(chrom)].append(float(pred[10]))

    for key in chromosight_centers.keys():
        print('chr_' + key, "prediction loops", len(chromosight_centers[key]))
    return chromosight_centers

def load_HiCExplorer(args):
    explorer_path = '/home/chens0a/HiC-Yolo/baseline/' + '{}-{}/'.format(args.experiment, args.cell_type) \
                       + '{}-{}-hicexplorer.bedgraph'.format(args.experiment, args.cell_type)

    explorer_loops = read_pred_loops(explorer_path)
    print(f'total number of loops predicted by HiCExplorer = {len(explorer_loops)}')


    explorer_centers = {}


    for chrom in tqdm(range(args.chr_start, args.chr_end)):

        explorer_centers[str(chrom)] = []


        for pred in explorer_loops:
            print()
            # pred[0] is in format of '1', chrom is 1
            if pred[0] != str(chrom):
                continue
            x_center = int((int(pred[1]) + int(pred[2])) / 2)
            y_center = int((int(pred[4]) + int(pred[5])) / 2)

            explorer_centers[str(chrom)].append([x_center, y_center])


    for key in explorer_centers.keys():
        print('chr_' + key, "prediction loops", len(explorer_centers[key]))
    return explorer_centers

def load_HICCUPS(args):
    hiccups_path = '/home/chens0a/HiC-Yolo/baseline/' + '{}-{}/'.format(args.experiment, args.cell_type) \
                    + '{}-{}-hiccups/'.format(args.experiment, args.cell_type)+'merged_loops.bedpe'


    hiccups_loops = read_pred_loops(hiccups_path)
    print(f'total number of loops predicted by HiCExplorer = {len(hiccups_loops)}')

    hiccups_centers = {}

    for chrom in tqdm(range(int(args.chr_start), int(args.chr_end))):

        hiccups_centers[str(chrom)] = []

        for pred in hiccups_loops:
            # pred[0] is in format of '1', chrom is 1
            if pred[0] != str(chrom) or pred[1]=='x1':
                continue


            x_center = int((int(pred[1]) + int(pred[2])) / 2)
            y_center = int((int(pred[4]) + int(pred[5])) / 2)

            hiccups_centers[str(chrom)].append([x_center, y_center])

    for key in hiccups_centers.keys():
        print('chr_' + key, "prediction loops", len(hiccups_centers[key]))
    return hiccups_centers





def get_Final_result(eval_dic,gt_centers_dict):
    final_result = {
        "Recall": 0,
        "Precision": 0,
        "F1 Score": 0,
    }

    total_TP = 0
    total_Pred = 0
    total_loops = 0

    for chrom in eval_dic.keys():
        chrom_rec = eval_dic[chrom]["Recall"]
        chrom_pre = eval_dic[chrom]["Precision"]
        chrom_f1 = eval_dic[chrom]["F1"]

        num_chrom_gt_loops = len(gt_centers_dict[chrom])
        TP = chrom_rec * num_chrom_gt_loops

        if chrom_pre == 0:
            Num_Pred = 0
        else:
            Num_Pred = TP / chrom_pre


        total_TP += TP
        total_Pred += Num_Pred
        total_loops += num_chrom_gt_loops




    final_result["Recall"] = total_TP / total_loops

    if total_Pred !=0:
        final_result["Precision"] = total_TP / total_Pred
    else:
        final_result["Precision"] = 0


    if final_result["Precision"] == 0 and final_result["Recall"] == 0:
        final_result["F1 Score"] = 0
    else:
        final_result["F1 Score"] = 2 * final_result["Recall"] * final_result["Precision"] / (
                    final_result["Precision"] + final_result["Recall"])

    return final_result, total_TP, total_loops, total_Pred



if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser(description='dfsa')
    parser.add_argument('--method', default='HiC-LDNet', help='method')
    parser.add_argument('--experiment', default='hic', help='experiment')
    parser.add_argument('--cell_type', default='gm12878', help='cell_type')
    parser.add_argument('--threshold', default=0.6, help='threshold')
    parser.add_argument('--chr_start', default=1, help='chr_start')
    parser.add_argument('--chr_end', default=23, help='chr_end')
    parser.add_argument('--window_size', default=512, help='Only use for HiC-LDNet')
    parser.add_argument('--resolution', default=100000, help='resolution')

    args = parser.parse_args()


    # Load Ground Truth Loops
    gt_center_dic = load_grount_truth(args.cell_type, args)

    pred_center_dic = None
    if args.method == "HiC-LDNet":
        pred_center_dic = load_HiC_LDNet(args)
    elif args.method == "peakachu":
        pred_center_dic = load_Peakachu(args)
    elif args.method == "chromosight":
        pred_center_dic = load_chromosight(args)
    elif args.method == "hicexplorer":
        pred_center_dic = load_HiCExplorer(args)
    elif args.method == "hiccups":
        pred_center_dic = load_HICCUPS(args)

    Eval_Dic = {}

    for chrom in range(int(args.chr_start), int(args.chr_end)):
        Eval_Dic[str(chrom)] = {}

        pred_centers = pred_center_dic[str(chrom)]
        print('=> Running chromosom {}'.format(chrom))
        recall, precision = get_recall_precision(gt_center_dic[str(chrom)],
                                                 pred_centers,
                                                 args.resolution)

        if precision == 0 and recall  == 0:
            f1 = 0
        else:
            f1 = 2 * recall * precision / (precision + recall)
        Eval_Dic[str(chrom)]["Recall"] = recall
        Eval_Dic[str(chrom)]["Precision"] = precision
        Eval_Dic[str(chrom)]["F1"] = f1


    Final_Dic, total_TP, total_loops, total_Pred = get_Final_result(Eval_Dic,gt_center_dic)

    file_path = os.path.join("/home/chens0a/HiC-Yolo/results",
                             args.experiment + '_' + args.cell_type)
    if not os.path.exists(file_path):
        os.mkdir(file_path)


    if args.method != "HiC-LDNet":
        file = os.path.join(file_path,args.method + ".json")
    else:
        file = os.path.join(file_path, args.method+"_"+ str(args.window_size) + "_.json")

    if os.path.exists(file):
        data = js.load(open(file))
        data[args.threshold] = [Final_Dic["Recall"],
                                Final_Dic["Precision"],
                                Final_Dic["F1 Score"],
                                total_TP,
                                total_loops,
                                total_Pred
                                ]

    else:
        data = {}
        data[args.threshold] = [Final_Dic["Recall"],
                                Final_Dic["Precision"],
                                Final_Dic["F1 Score"],
                                total_TP,
                                total_loops,
                                total_Pred
                                ]

    with open(file, "w") as f:
        js.dump(data, f,indent=4)
        print("=> {} wirtten".format(file))