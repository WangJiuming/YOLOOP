# !/bin/bash -v

# Ours GM12878
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --window_size=512 --method='HiC-LDNet' --experiment='hic' --cell_type='gm12878' --threshold=$threshold
#done

# Ours K562
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --window_size=512 --method='HiC-LDNet' --experiment='hic' --cell_type='k562' --threshold=$threshold
#done


# Ours IMR-90
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --window_size=512 --method='HiC-LDNet' --experiment='hic' --cell_type='imr90' --threshold=$threshold
#done

# Ours KBM7
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --window_size=512 --method='HiC-LDNet' --experiment='hic' --cell_type='kbm7' --threshold=$threshold
#done

# Ours HUVEC
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --window_size=512 --method='HiC-LDNet' --experiment='hic' --cell_type='huvec' --threshold=$threshold
#done

#HICHIP GM12878
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --window_size=512 --method='HiC-LDNet' --experiment='hichip' --cell_type='gm12878' --threshold=$threshold
#done

#ODC SingleCell
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --window_size=512 --method='HiC-LDNet' --experiment='hic' --cell_type='odc' --threshold=$threshold
#done

# Chromosight
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --method='chromosight' --experiment='hic' --cell_type='odc' --threshold=$threshold
#done

#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --method='chromosight' --experiment='hic' --cell_type='k562' --threshold=$threshold
#done
##
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --method='chromosight' --experiment='hic' --cell_type='imr90' --threshold=$threshold
#done
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --method='chromosight' --experiment='hic' --cell_type='huvec' --threshold=$threshold
#done

#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --method='chromosight' --experiment='hic' --cell_type='kbm7' --threshold=$threshold
#done



#python ROC_Evaluation.py --method='hicexplorer' --experiment='hic' --cell_type='gm12878' --threshold=0.1
#python ROC_Evaluation.py --method='hicexplorer' --experiment='hic' --cell_type='k562' --threshold=0.1
#python ROC_Evaluation.py --method='hicexplorer' --experiment='hic' --cell_type='imr90' --threshold=0.1
#python ROC_Evaluation.py --method='hicexplorer' --experiment='hic' --cell_type='huvec' --threshold=0.1


#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --method='peakachu' --experiment='hic' --cell_type='gm12878' --threshold=$threshold  --chr_start=1 --chr_end=23
#done
#
#
#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --method='peakachu' --experiment='hic' --cell_type='imr90' --threshold=$threshold --chr_start=1 --chr_end=23
#done

#for threshold in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
#do
#python ROC_Evaluation.py --method='peakachu' --experiment='hic' --cell_type='huvec' --threshold=$threshold --chr_start=1 --chr_end=23
#done


for threshold in 0.3 0.4 0.5 0.6
do
python ROC_Evaluation.py --method='peakachu' --experiment='hic' --cell_type='odc' --threshold=$threshold --chr_start=1 --chr_end=23
done

# HICCUPS OCD
#python ROC_Evaluation.py --method='hiccups' --experiment='hic' --cell_type='odc' --threshold=$threshold --chr_start=1 --chr_end=23


