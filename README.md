# YOLOOP:  Fast and Adaptive Detection Framework for Genome-wide Chromatin Loop Mapping from Hi-C data


## Installation

First, download the code from GitHub and move to the working directory.

```bash
git clone https://github.com/WangJiuming/YOLOOP.git
```
```bash
cd YOLOOP
```

Then, install the prerequisite packages from the ```environment.yml```, which specifies the required packages for using YOLOOP. We recommend use anaconda for installing the dependencies into a virtual environment.

```bash
conda create --name yoloop --file environment.yml
```

After installing all the packages successfully, you may proceed to activate the environment as usual. The name of the virtual environment is "yoloop".

```bash
conda activate yoloop
```

Now, you are ready to use YOLOOP!

## Use YOLOOP for chromatin loop detection
### download sample data
YOLOOP performs extremely efficient loop detection across contact maps obtained with various sequencing protocols and from various cell lines. All the datasets (i.e., contact maps and loop annotations) are in the public domain. Their sources and access numbers are listed in the Supplementary Information of our paper. In the meanwhile, please also feel free to use any of your own datasets!

YOLOOP supports one of the currently most commonly used file formats of contact maps, cooler. You can find more about it at its official [documentation](https://cooler.readthedocs.io/en/latest/index.html). For a better performance, we highly recommend use it for less memory IO overhead. If only hic format is available, you may also check out [this](https://github.com/4dn-dcic/hic2cool) very convenient tool to convert it from hic to cool.

In the following tutorial, we will use the [GM12878 dataset](https://data.4dnucleome.org/files-processed/4DNFIXP4QG5B/) downloaded from 4DN Portal with access number 4DNFIXP4QG5B.

After downloading the dataset, you may modify the path to the .hic and .bedpe files in ```config.json```.

### download the pretrained model checkpoint

To help users conveniently use our model, we offer a variety of pretrained model checkpoints optimized on different datasets. You may download the checkpoint of interest directly [here](https://drive.google.com/drive/folders/1yyqtltWRwDi-YRTHjii7hD1W08XiUevf?usp=sharing), which will save the user tremendous amount of time by avoiding training the model from scratch. The model was trained with a window size of 512 for 100 epochs.

### genome-wide chromatin loop detection with YOLOOP

After obtaining the model checkpoint, you are ready to perform chromatin loop detection efficiently with YOLOOP by calling the detection procedure. A standard calling would be as the following.
```
python detect.py --map <path_to_contact_map> --weight <path_to_model_checkpoint> --ouput <path_to_output_dir> 
```
The program will detect cuda devices automatically, and we strongly suggest use cuda for a much better performance.
Besides setting the paths, here are also several hyperparameters that we may tune. A complete configuration of the procedure would be as follows.
```
python detect.py --map <path_to_contact_map> --weight <path_to_model_checkpoint> --ouput <path_to_output_dir> --window 256 --threshold 0.5
```

### analyze the prediction results

After the detection is complete, the results will be saved in a .bedpe file in the specified directory. An example of the prediction results is as follows.

```
chr1	610000	620000	chr1	37880000	37890000	0.760546875
```

The above example consists of seven columns. The first three columns indicate the x-coordiante of the loop and the following three columns indicate the y-coordinate. The last column shows the confidence level of the prediction.

To further understand the prediction results, we offer a convenient visualization tool to visualize the prediction results by mapping the loops back to the contact map. To use the program, simply use the following command.

```
python visualize.py --map <path_to_contact_map> --loop <path_to_chr_loop> --threshold 0.5 --chr 1 --start 210000000 --end 212000000 --resolution 10000
```

The visualize.py program requires several command line arguments. Besides paths to the contact map file (.cool or .mcool) and the loop prediction file (.bedpe), users may also threshold the loops based on the prediction confidence level. For example, using the command line above, the program will plot the region on chromosome 1 from 210,000,000 to 212,000,000 under resolution 10kb. The output is shown below.

<img src="https://github.com/WangJiuming/YOLOOP/blob/main/images/pred.png" width="500">

For comparison between ground truth labels and predictions, simply set the compare flag in the visualize.py program to be True (which, by default, is False).

```
python visualize.py --compare 1 --map <path_to_contact_map> --labels <path_to_labels> --preds <path_to_predictions> --threshold 0.5 --chr 1 --start 210000000 --end 212000000 --resolution 10000
```
The ground truth labels and the model's predictions will be shown on each side of the diagnol for easy comparison. The result of the above command is shown here.

<img src="https://github.com/WangJiuming/YOLOOP/blob/main/images/compare.png" width="500">

The black circles in upper-right part of the image represents predictions, while the yellow circles in the lower-left part represents the labels. As is shown, YOLOOP can accurately recover the ground truth chromatin loops.

