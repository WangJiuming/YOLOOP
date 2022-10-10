import numpy as np

import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torchvision import transforms


class ContactMapDataset(Dataset):
    def __init__(self, path_prefix, transform):
        self.data = np.load(f'{path_prefix}-map.npy', allow_pickle=True)
        self.boxes = np.load(f'{path_prefix}-loop.npy', allow_pickle=True)
        print(f'self.boxes.shape = {self.boxes[0].shape}')
        self.transform = transform
        self.size = self.data.shape[1]

        self.aux = np.array([5, 5, 1, 1]).reshape((1, 4))

    def __len__(self):
        return self.data.shape[0]

    def __getitem__(self, item):
        return self.transform(self.data[item]).reshape((1, self.size, self.size)), \
               self.transform(np.hstack((self.boxes[item],
                                         np.broadcast_to(self.aux, (self.boxes[item].shape[0], 4)))))
        # format the label, normalize the xywh
