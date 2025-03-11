# A Fast and Adaptive Detection Framework for Genome-wide Chromatin Loop Mapping from Hi-C data

![header](header.png)

This is the open source code for YOLOOP. 

Note: The current version of our code is for a demo to facilitate the review process. All source materials to develop our model will be released upon publication, including all the codes, data and models for pre-training, evaluation and adaptive fine-tuning.


## 1. Installation

First, download the code from GitHub and move to the working directory.

```bash
git clone https://github.com/WangJiuming/YOLOOP.git
```
```bash
cd YOLOOP
```

Then, install the prerequisite packages from the ```environment.yml```, which specifies the required packages for using YOLOOP. We recommend use anaconda for installing the dependencies into a virtual environment.

```bash
conda env create --name yoloop --file environment.yml
```

By default, we installed the latest PyTorch with CUDA version 12.1. If your local environment configurations is different, please install based on the [PyTorch Installation Guide](https://pytorch.org/get-started/locally/).

After installing all the packages successfully, you may proceed to activate the environment as usual. The name of the virtual environment is "yoloop".

```bash
conda activate yoloop
```

Now, you are ready to use YOLOOP!

## 2. Use YOLOOP for chromatin loop detection
### 2.1 Download sample data
YOLOOP performs extremely efficient loop detection across contact maps obtained with various sequencing protocols and from various cell lines. All the datasets (i.e., contact maps and loop annotations) are in the public domain. Their sources and access numbers are listed in the Supplementary Information of our paper. In the meanwhile, please also feel free to use any of your own datasets!

YOLOOP supports one of the currently most commonly used file formats of contact maps, cooler. You can find more about it at its official [documentation](https://cooler.readthedocs.io/en/latest/index.html). For a better performance, we highly recommend use it for less memory IO overhead. If only hic format is available, you may also check out [this](https://github.com/4dn-dcic/hic2cool) very convenient tool to convert it from hic to cool.

In the following tutorial, we will use the [GM12878 dataset](https://data.4dnucleome.org/files-processed/4DNFIXP4QG5B/) downloaded from 4DN Portal with access number 4DNFIXP4QG5B.
```bash
mkdir -p data
wget -P ./data https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/d6abea45-b0bb-4154-9854-1d3075b98097/4DNFIXP4QG5B.mcool
```

### 2.2 Download the pre-trained model checkpoint

[//]: # (To help users conveniently use our model, we offer a variety of pretrained model checkpoints optimized on different datasets. You may download the checkpoint of interest directly [here]&#40;https://drive.google.com/drive/folders/1yyqtltWRwDi-YRTHjii7hD1W08XiUevf?usp=sharing&#41;, which will save the user tremendous amount of time by avoiding training the model from scratch. The model was trained with a window size of 512 for 100 epochs.)
We have included one model checkpoint at ```./models/gm12878_hic_10kb.pt```. This model was pre-trained on the GM12878 Hi-C contact map with CTCF ChIA-PET interactions at 10kb resolution.

### 2.3 Genome-wide chromatin loop detection with YOLOOP

After obtaining the model checkpoint, you are ready to perform chromatin loop detection efficiently with YOLOOP by calling the detection procedure. A standard calling would be as the following.
```bash
python detect.py --cm ./data/4DNFIXP4QG5B.mcool --r 10000 --model ./models/gm12878_hic_10kb.pt --out ./results 
```
The program will detect cuda devices automatically, and we strongly suggest use cuda for a much better performance.
Besides setting the paths, here are also several hyperparameters that we may tune. A complete configuration of the procedure would be as follows.
```
usage: detect.py [-h] [--cm CM] [-r RESOLUTION] [-b] [-m MODEL] [--out OUT]
                 [-t THRESH] [--device DEVICE]

YOLOOP for efficient chromatin loop detection

optional arguments:
  -h, --help                show this help message and exit
  --cm CM                   path to the input contact matrix with .mcool/.cool extension
  -r RESOLUTION, --resolution RESOLUTION    resolution of the contact matrix
  -b, --balance             whether to use a balanced matrix or not
  -m MODEL, --model MODEL   YOLOOP model checkpoint to be loaded
  --out OUT                 output directory for saving the prediction results
  -t THRESH, --thresh THRESH    threshold for the confidence score
  --device DEVICE           device to be used, e.g., cuda, cuda:0, cpu
```

### 2.4 Analyze the prediction results

After the detection is complete, the results will be saved in a .bedpe file in the specified directory. An example of the prediction results is as follows.

```
chr1	610000	620000	chr1	37880000	37890000	0.760546875
```

The above example consists of seven columns. The first three columns indicate the x-coordiante of the loop and the following three columns indicate the y-coordinate. The last column shows the confidence level of the prediction.

### 2.5 Motif logo analysis

Reproducible run for custom scripts can be found under ```reproducibility/```.


## 3. Training your own model with YOLOOP

We provide an example for you to train your own detection model with YOLOOP. We use the K562 cell line with CTCF-supported loops. 

### 3.1 Training Data Preparation

run
```bash
gen_train_data.sh
```
in which you need to provide the path for the Hi-c contact file in .mcool format in  ```data_dir```, and the .bedpe file for positive chromatin interaction in ```loop_dir```. For reference, the ground truth positive interactions we use is in [reproducibility/ground_truth](https://github.com/WangJiuming/YOLOOP/tree/main/reproducibility/ground_truth)
The 
It will take several minutes depending on the size of the .mcool file.
After that, you will get a data folder containing the training samples with certain binding factors.

```
└── CTCF
    └── k562
        ├── images
        └── labels
```


### 3.2 Train

Take K562 cell line as an example, we now train a YOLOOP model on the K562 cell line with orthogonal ChIA-PET CTCF binding factors.

```bash
CUDA_VISIBLE_DEVICES=0 python train.py -f config_train_k562.json
```

### Complete Data and Training

You are able to train YOLOOP on multiple cell lines with .mcool files on [GM12878](https://data.4dnucleome.org/files-processed/4DNFIXP4QG5B/), [IMR-90](https://data.4dnucleome.org/files-processed/4DNFIJTOIGOI/), [HCT116](https://data.4dnucleome.org/files-processed/4DNFILP99QJS/), and [K562](https://data.4dnucleome.org/files-processed/4DNFITUOMFUQ/). For those who wish to train on your own ```.hic``` file as raw contact matrix, you are able to convert them into ```.mcool``` file using [hic2cool](https://pypi.org/project/hic2cool/0.4.0/). The detection framework should be trained with multiple binding factor annotations as [ground_truth](https://github.com/WangJiuming/YOLOOP/tree/main/reproducibility/ground_truth)
The data folder should appear as follows after conversion.
```angular2html
.├── CTCF
│   ├── gm12878
│   ├── hct-116
│   ├── imr90
│   └── k562
├── RAD21
│   ├── gm12878
│   └── mesc
└── SMC1
    └── sc_hic

```
