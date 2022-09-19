#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2022/4/3
# @Author   : Wang Jiuming
# @File     : dataloader.py
# @Software : PyCharm

import numpy as np

import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torchvision import transforms

import matplotlib.pyplot as plt
import random
import cooler
import straw
import json

from docopt import docopt


# hyperparameters
# batch_size = 8


class ContactMapDataset(Dataset):
    def __init__(self, path_prefix, transform):
        self.data = np.load(f'{path_prefix}-map.npy', allow_pickle=True)
        self.boxes = np.load(f'{path_prefix}-loop.npy', allow_pickle=True)
        print(f'self.boxes.shape = {self.boxes[0].shape}')
        self.transform = transform
        self.size = self.data.shape[1]

        self.aux = np.array([5, 5, 1, 1]).reshape((1, 4))
        # print(np.broadcast_to(self.meta, (4, 4)))

    def __len__(self):
        return self.data.shape[0]

    def __getitem__(self, item):
        # print(f'here {self.boxes[item].shape}, {np.broadcast_to(self.meta, (self.boxes[item].shape[0], 4)).shape}')
        return self.transform(self.data[item]).reshape((1, self.size, self.size)), \
               self.transform(np.hstack((self.boxes[item],
                                         np.broadcast_to(self.aux, (self.boxes[item].shape[0], 4)))))
        # format the label, normalize the xywh


# if __name__ == '__main__':
#     path = '/Users/jmwang/Desktop/AIST4010 Foundation of Applied Deep Learning/project/data/chr1'
#
#     train_ds = ContactMapDataset(path, transforms.Compose([transforms.ToTensor()]))
#
#     train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
#
#     amap, box = next(iter(train_loader))
#
#     print(amap.shape, box.shape)
