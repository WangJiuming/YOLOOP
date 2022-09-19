# !/bin/bash -v

#GM12878
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/hic/gm12878/hic-gm12878.mcool --weight=/home/chens0a/yolov3-master/hic-gm12878-runs/train/exp/weights/best.pt


#K562
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/hic/k562/hic-k562.mcool --weight=/home/chens0a/yolov3-master/hic-k562-runs/train/exp/weights/best.pt

#IMR-90
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/hic/imr90/hic-imr90.mcool --weight=/home/chens0a/yolov3-master/hic-imr90-runs/train/exp/weights/best.pt

#KBM7
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/hic/kbm7/hic-kbm7.mcool --weight=/home/chens0a/yolov3-master/hic-kbm7-runs/train/exp/weights/best.pt

#HUVEC
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/hic/huvec/hic-huvec.mcool --weight=/home/chens0a/yolov3-master/hic-huvec-runs/train/exp/weights/best.pt




# CHIA-PET DNA-SPRITE HICHIP
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/chia-pet/gm12878/chia-pet-gm12878.mcool --weight=/home/chens0a/yolov3-master/hic-gm12878-runs/train/exp/weights/best.pt
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/dna-sprite/gm12878/dna-sprite-gm12878.mcool --weight=/home/chens0a/yolov3-master/hic-gm12878-runs/train/exp/weights/best.pt
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/hichip/gm12878/hichip-gm12878.mcool --weight=/home/chens0a/yolov3-master/hic-gm12878-runs/train/exp/weights/best.pt

# SC-HiC ODC
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/hic/odc/hic-odc.mcool --weight=/home/chens0a/yolov3-master/hic-gm12878-runs/train/exp/weights/best.pt
time python detect.py --window=512 --threshold=0.1 --device='cuda' --hic_file=/Dataset/HiC/hic/odc/hic-odc.mcool --weight=/home/chens0a/yolov3-master/hic-odc-runs/train/exp/weights/best.pt




