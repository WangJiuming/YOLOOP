"""
@File   : train.py
@Author : Siyuan Chen
@Date   : 2025/3/10
@Desc   : 
"""

import os
import sys
import json
import argparse

os.environ["OMP_NUM_THREADS"] = '16'

import comet_ml
from comet_ml import Experiment

# add our local modified YOLO project directory to the top of Python paths
sys.path.insert(0, os.path.abspath('./ultralytics/'))

from ultralytics import YOLO

if __name__ == '__main__':
    # 传递参数
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',
                        '--file',
                        help='Config file')
    args = parser.parse_args()
    config_path = args.file

    with open(config_path, 'r') as config_file:
        config = json.load(config_file)

    dataset_cfg_path = config['paths']['dataset']  # path to the .yaml configuration file for the training dataset
    model_cfg_path = config['paths']['model']  # path to the .yaml configuration file for the YOLO model

    dataset_name = os.path.splitext(os.path.basename(dataset_cfg_path))[0]

    print(f'[info] dataset configuration: {dataset_name}')
    print(f'[info] model configuration: {model_cfg_path}')

    epochs = config['params']['epochs']
    save_period = config['params']['save period']
    device = config['params']['device']
    batch_size = config['params']['batch size']
    print('[info] device:', device)
    print(f'[info] total epochs: {epochs}')
    print(f'[info] save checkpoint every {save_period} epochs')
    print(f'[info] batch size: {batch_size}')

    window_size = config['params']['window size']
    print(f'[info] window size: {window_size}')

    # initialize logger
    # create an experiment with your api key
    experiment = Experiment(
        api_key="dzIeausd4DXyJ5uXl3JXbNOx4",
        project_name=config['params']['project name'],
    )

    comet_ml.init()

    # initialize the model from scratch
    model = YOLO(model_cfg_path)  # build a new model from .yaml
    print(model.info(detailed=True))  # print the model architecture, check if the input layer has 1 channel

    # train and validate the model (train-val split is specified in the dataset configuration file)
    print(f'[info] start training')
    print(f'[info] training results to be saved at: yoloop/{dataset_name}/')

    results = model.train(data=dataset_cfg_path,
                          epochs=epochs,
                          save_period=save_period,
                          imgsz=window_size,
                          batch=batch_size,
                          device=device,
                          verbose=True,
                          project='yoloop178',
                          name=dataset_name,
                          pretrained=False,
                          patience=10,
                          )

    experiment.end()

    print(f'[info] training: done')


