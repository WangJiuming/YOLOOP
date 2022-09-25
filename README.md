# YOLOOP: YOLO for chromatin LOOP detection - an efficient and robust approach


## Installation

First, install the prerequisite packages from the yoloop-spec.yml, which specifies the required packages for using YOLOOP. We suggest use anaconda for installing the dependencies into a virtual environment.

```bash
conda create -f env-spec.yml
```

Then, activate the environment as usual. The name of the virtual environment is "yoloop" by default.

```bash
conda activate yoloop
```

After installing all the packages successfully, you may proceed to download the code directly from GitHub followed by moving to the working directory.

```bash
git clone https://github.com/WangJiuming/YOLOOP.git
```
```bach
cd YOLOOP
```

Now, you are ready to use YOLOOP!

## Use YOLOOP for chromatin loop detection
### download sample data
YOLOOP performs extremely efficient loop detection across contact maps obtained with various sequencing protocols and various cell lines. All of the datasets (ie contact maps and loop annotations) are in the public domain. Their source and access number are listed in the Supplementary of our paper. In the meanwhile, please also feel free to use any of your own datasets!

YOLOOP supports the two currently most commonly used two file formats of contact maps, namely cooler and hic. For a better performance, we suggest use cooler format over hic for less memory IO overhead. If only hic format is available, you may also check out this very convenient tool to convert it from hic to cool.

In the following tutorial, we will use the [HUVEC dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) downloaded from GEO with access number GSE63525.

After downloading the dataset, you may modify the path to the .hic and .bedpe files in dataloader.json.

### download the pretrained model checkpoint

To help users conveniently use our model, we offer a variety of pretrained model checkpoints optimized on different datasets. You may download the checkpoint of interest directly [here](https://drive.google.com/drive/folders/1yyqtltWRwDi-YRTHjii7hD1W08XiUevf?usp=sharing), which will save the user tremendous amount of time by avoiding training the model from scratch. 

### genome-wide chromatin loop detection with YOLOOP

After obtaining the model checkpoint, you are ready to perform chromatin loop detection efficiently with YOLOOP.

(to be finished later)

### analyze the prediction results

The output format of the prediction results is as follows.

## Train YOLOOP from scratch

### pre-process the training data

TBC

procedure:
1. modify the path to .hic and .bedpe files in dataloader.json
2. run "$ python data.py data.json" to generate chr{i}-map.npy and chr{i}-loop.npy, output saved at ./data/
3. run "$ python image.py" to convert .npy to .png and .txt, output will be saved at ./hic-gm12878/
   (no need to change the code, the program will read data from ./data/ directly, empty output directories are provided)
4. use data.ipynb to check if the converted images and labels are correct
5. download yolov3 from GitHub via "git clone https://github.com/ultralytics/yolov3"
6. move hic-gm12878.yaml to ./yolov3/
7. "$ cd yolov3"
8. "$ pip install -r requirements.txt"
9. "$ python train.py --img 256 --batch 16 --epoch 100 --data hic-gm12878.yaml --weights '' --cfg yolov3.yaml"
   (we can modify the "weights" argument to a path to a .pt file to resume training)
10. training results will be saved at ./yolov3/runs/train/exp{i}





###########################################################################









