import numpy as np
import matplotlib.pyplot as plt
import cooler
import argparse

from data import load_loops, load_contact_map


def show(mat):
    """visualize a contact map"""
    figure = plt.figure(figsize=(10, 10))
    axes = figure.add_subplot(111)
    caxes = axes.matshow(np.log10(mat), cmap='Reds')
    plt.show()


def show_loop(mat, loops, start, end, chr_num, resolution=10000):
    """visualize a contact map together with loop predictions"""
    # plot the contact map
    figure = plt.figure(figsize=(10, 10))
    axes = figure.add_subplot(111)
    caxes = axes.matshow(np.log10(mat), cmap='Reds')

    # add the loops
    set_label = False
    for loop in loops:
        circle = plt.Circle(((loop[1]-start)/resolution, (loop[2]-start)/resolution), 1.5,
                            color='black', fill=False, lw=1)
        axes.add_patch(circle)
        if not set_label:
            circle.set_label('predictions')
            set_label = True

    # plt.show()
    plt.legend()
    plt.title(f'Chr {chr_num}: {start} to {end}\nresolution: {resolution}')
    plt.axis('off')
    plt.savefig('pred.png')


def show_validate(mat, preds, labels, start, end, chr_num, resolution=10000):
    # plot the contact map
    figure = plt.figure(figsize=(10, 10))
    axes = figure.add_subplot(111)
    caxes = axes.matshow(np.log10(mat), cmap='Reds')

    # remove symmetric loops in the loop lists
    preds = [loop for loop in preds if loop[1] > loop[2]]  # upper-right
    labels = [loop for loop in labels if loop[1] < loop[2]]  # lower-left

    # add the predicted loops
    set_label = False
    for loop in preds:
        pred_circle = plt.Circle(((loop[1]-start)/resolution, (loop[2]-start)/resolution), 1.5,
                                 color='black', fill=False, lw=1.5)
        axes.add_patch(pred_circle)
        if not set_label:
            pred_circle.set_label('predictions')
            set_label = True

    # add the ground truth loops
    set_label = False
    for loop in labels:
        label_circle = plt.Circle(((loop[1]-start)/resolution, (loop[2]-start)/resolution), 1.5,
                                  color='yellow', fill=False, lw=1.5)
        axes.add_patch(label_circle)
        if not set_label:
            label_circle.set_label('labels')
            set_label = True

    # plt.show()
    plt.legend()
    plt.title(f'Chr {chr_num}: {start} to {end}\nresolution: {resolution}')
    plt.axis('off')
    plt.savefig('compare.png')


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
            raise Exception('MissingArguments: both labels and predictions need to be provided for comparison')
        labels = load_loops(args.labels, args.chr_num, int(args.start), int(args.end))
        preds = load_loops(args.preds, args.chr_num, int(args.start), int(args.end))
        show_validate(contact_map, preds, labels, int(args.start), int(args.end), args.chr_num, int(args.resolution))
    else:
        if args.preds == '':
            raise Exception('MissingArguments: predictions need to be provided for plotting')
        loops = load_loops(args.preds, args.chr_num, int(args.start), int(args.end))
        show_loop(contact_map, loops, int(args.start), int(args.end), args.chr_num, int(args.resolution))
