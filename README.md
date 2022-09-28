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
```bash
cd YOLOOP
```

Now, you are ready to use YOLOOP!

## Use YOLOOP for chromatin loop detection
### download sample data
YOLOOP performs extremely efficient loop detection across contact maps obtained with various sequencing protocols and from various cell lines. All of the datasets (ie contact maps and loop annotations) are in the public domain. Their source and access number are listed in the Supplementary of our paper. In the meanwhile, please also feel free to use any of your own datasets!

YOLOOP supports the two currently most commonly used two file formats of contact maps, namely cooler and hic. For a better performance, we suggest use cooler format over hic for less memory IO overhead. If only hic format is available, you may also check out [this](https://github.com/4dn-dcic/hic2cool) very convenient tool to convert it from hic to cool.

In the following tutorial, we will use the [HUVEC dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) downloaded from GEO with access number GSE63525.

After downloading the dataset, you may modify the path to the .hic and .bedpe files in dataloader.json.

### download the pretrained model checkpoint

To help users conveniently use our model, we offer a variety of pretrained model checkpoints optimized on different datasets. You may download the checkpoint of interest directly [here](https://drive.google.com/drive/folders/1yyqtltWRwDi-YRTHjii7hD1W08XiUevf?usp=sharing), which will save the user tremendous amount of time by avoiding training the model from scratch. The model was trained with a window size of 512 for 100 epochs.

### genome-wide chromatin loop detection with YOLOOP

After obtaining the model checkpoint, you are ready to perform chromatin loop detection efficiently with YOLOOP.

First, you need to configure the hyperparameters for the model in data.json. Specifically, you need to customize the path to the dataset.
```python
python config.py -detect <path_to_contact_map>
```
For detection (with the -detect flag), the path to loop annotations is not required.

Then, you are ready to call the detection procedure. The standard calling would be as the following.
```python
python detect.py --device cuda --map <path_to_contact_map> --weight <path_to_model_checkpoint> --ouput <path_to_output_dir> 
```
We strongly suggest use cuda for a much better performance.
Besides setting the paths, here are also a few hyperparameters that we may tune. A complete configuration of the procedure would be as follows.
```python
python detect.py --device cuda --map <path_to_contact_map> --weight <path_to_model_checkpoint> --ouput <path_to_output_dir> --window 256 --threshold 0.5
```

### analyze the prediction results

After the detection is complete, the results will be saved in a .bedpe file in the specified directory. An example of the prediction results is as follows.

```
chr1	610000	620000	chr1	37880000	37890000	0.760546875
```

The above example consists of seven columns. The first three columns indicate the x-coordiante of the loop and the following three columns indicate the y-coordinate. The last column shows the confidence level of the prediction.

To further understand the prediction results, we offer a convenient visualization tool to visualize the prediction results by mapping the loops back to the contact map. To use the program, simply use the following command.

```python
python visualize.py --map <path_to_contact_map> --loop <path_to_chr_loop> --threshold 0.5 --chr 1 --start 210000000 --end 212000000 --resolution 10000
```

The visualize.py program requires several command line arguments. Besides paths to the contact map file (.cool or .mcool) and the loop prediction file (.bedpe), users may also threshold the loops based on the prediction confidence level. For example, using the command line above, the program will plot the region on chromosome 1 from 210,000,000 to 212,000,000 under resolution 10kb. The output is shown below.



## Cite us
