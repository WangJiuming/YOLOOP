import numpy as np
import matplotlib.pyplot as plt


def show(mat, save=False):
    """visualize a contact map"""
    figure = plt.figure(figsize=(10, 10))
    axes = figure.add_subplot(111)
    caxes = axes.matshow(np.log10(mat), cmap='Reds')
    plt.axis('off')
    if save:
        plt.savefig('contact_map.png')
    else:
        plt.show()


def show_loop(mat, loops, start, end, chr_num, resolution=10000, save=False):
    """
    visualize a contact map together with loop calls
    :param mat: array, contact map to be plotted against
    :param loops: list of list, the loop coordinates in the region to be plotted
    :param start: int, starting coordinate of the range to be plotted
    :param end: int, ending coordinate of the range to be plotted
    :param chr_num: int, number of the chromosome, e.g., 12
    :param resolution: int, resolution of the contact map, 10kb by default
    :param save: bool, whether to save the figure or not, if set False, then plot to the screen, False by default
    :return: None
    """
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
    if save:
        plt.savefig('pred.png')
    else:
        plt.show()


def show_compare(mat, preds, labels, start, end, chr_num, resolution=10000, save=False):
    """
    visualize a contact map together with two sets of loop calls (preds and labels)
    :param mat: array, contact map to be plotted against
    :param preds: list of list, the first set of loop coordinates in the region to be plotted
    :param labels:list of list, the second set of loop coordinates in the region to be plotted
    :param start: int, starting coordinate of the range to be plotted
    :param end: int, ending coordinate of the range to be plotted
    :param chr_num: int, number of the chromosome, e.g., 12
    :param resolution: int, resolution of the contact map, 10kb by default
    :param save: bool, whether to save the figure or not, if set False, then plot to the screen, False by default
    :return: None
    """
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
    if save:
        plt.savefig('compare.png')
    else:
        plt.show()
