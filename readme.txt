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









