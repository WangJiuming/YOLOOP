import numpy as np
import matplotlib.pyplot as plt
import cooler
import argparse

from data import load_loops, load_contact_map
from utils.visualize_utils import show, show_loop, show_compare


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser(description='visualization')
    parser.add_argument('--map', default='', help='path to a contact map with .cool/.mcool extension')
    parser.add_argument('--labels', default='', help='optional, ground truth loop labels')
    parser.add_argument('--preds', default='', help='loop detection results')
    parser.add_argument('--threshold', default=0.0, help='threshold for prediction confidence, default: 0')
    parser.add_argument('--chr_num', default='1', help='which chromosome to plot')
    parser.add_argument('--start', default=None, help='plot the region from where')
    parser.add_argument('--end', default=None, help='plot the region to where')
    parser.add_argument('--resolution', default=10000, help='resolution')
    parser.add_argument('--compare', default=False, help='to show both labels and predictions or not')
    args = parser.parse_args()

    contact_map = load_contact_map(args.map, args.chr_num, int(args.start), int(args.end), int(args.resolution))
    print(f'The specified contact map region has size: {contact_map.shape} (in {args.resolution}b)')

    if int(args.compare):
        if args.labels == '' or args.preds == '':
            raise Exception('Missing Arguments: both labels and predictions need to be provided for comparison')
        labels = load_loops(args.labels, args.chr_num, int(args.start), int(args.end))
        preds = load_loops(args.preds, args.chr_num, int(args.start), int(args.end))
        show_compare(contact_map, preds, labels, int(args.start), int(args.end), args.chr_num, int(args.resolution))
    else:
        if args.preds == '':
            raise Exception('Missing Arguments: predictions need to be provided for plotting')
        loops = load_loops(args.preds, args.chr_num, int(args.start), int(args.end))
        show_loop(contact_map, loops, int(args.start), int(args.end), args.chr_num, int(args.resolution))
