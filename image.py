#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2022/4/12
# @Author   : Wang Jiuming
# @File     : image.py
# @Software : PyCharm

import numpy as np
from PIL import Image
import os
import json

window_size = 512
box_size = 10
cell = 'odc'
# cell = 'h1esc'
protocol = 'hic'

root = None

count = 0
train_val_split = 18  # 4-1 train-test split for 22 chromosomes


if not os.path.exists("{}-{}".format(protocol,cell)):
    os.mkdir("{}-{}".format(protocol,cell))



img_dir = os.path.join("{}-{}/images".format(protocol,cell))
label_dir = os.path.join("{}-{}/labels".format(protocol,cell))
if not os.path.exists(img_dir):
    os.mkdir(img_dir)
if not os.path.exists(label_dir):
    os.mkdir(label_dir)

if not os.path.exists(img_dir + "/train"):
    os.mkdir(img_dir + "/train")
if not os.path.exists(img_dir + "/val"):
    os.mkdir(img_dir + "/val")
    if not os.path.exists(label_dir + "/train"):
        os.mkdir(label_dir + "/train")
        if not os.path.exists(label_dir + "/val"):
            os.mkdir(label_dir + "/val")



for i in range(1, 23):
    hic_maps = np.load('data/{}-{}/chr{}-map.npy'.format(protocol,cell,i),
                       allow_pickle=True)
    chr_dic = json.load(open('data/{}-{}/chr{}-loop.json'.format(protocol,cell,i)))
    chr_loop_list = chr_dic[str(i)]



    map_num = hic_maps.shape[0]

    # for every cropped contact map window in this chromosome
    for j in range(map_num):
        hic_img = Image.fromarray(hic_maps[j])
        np_array = np.array(hic_img)

        image_loop = chr_loop_list[j]

        loop_num = len(image_loop)

        if i < train_val_split:
            hic_img.convert('L').save(f'{protocol}-{cell}/images/train/{count:08}.png')

            with open(f'{protocol}-{cell}/labels/train/{count:08}.txt', 'w') as file:

                # for every loop in this cropped window
                for k in range(0, loop_num):
                    # class, x center, y center, width, height
                    # N.B. reverse the order of x, y in numpy
                    file.write('{0:<4}{1:<15}{2:<15}{3:<15}{3:<15}\n'.format(0,
                                                                             image_loop[k][1] / window_size,
                                                                             image_loop[k][0] / window_size,
                                                                             box_size / window_size))
        else:
            hic_img.convert('L').save(f'{protocol}-{cell}/images/val/{count:08}.png')

            with open(f'{protocol}-{cell}/labels/val/{count:08}.txt', 'w') as file:

                for k in range(0, loop_num):
                    # start from 1, skip the first loop (center loop is repeated)
                    # class, x center, y center, width, height
                    # N.B. reverse the order of x, y in numpy
                    file.write('{0:<4}{1:<15}{2:<15}{3:<15}{3:<15}\n'.format(0,
                                                                             image_loop[k][1]/ window_size,
                                                                             image_loop[k][0]/ window_size,
                                                                             box_size / window_size))
        count += 1
    print(f'chromosome {i} done')
