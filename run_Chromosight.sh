# !/bin/bash -v

time chromosight detect --threads 1 hic-k562.mcool::/resolutions/10000 default_1_thread
time chromosight detect --threads 1 hic-gm12878.mcool::/resolutions/10000 default_1_thread
time chromosight detect --threads 8 hic-imr90.mcool::/resolutions/10000 default_8_thread

time chromosight detect --threads 8 hic-kbm7.mcool::/resolutions/10000 default_8_thread
time chromosight detect --threads 8 hic-huvec.mcool::/resolutions/10000 default_8_thread

time chromosight detect --threads 1 --min-dist 20000 --max-dist 200000 human-hic-gm12878.mcool::/resolutions/10000 20k-200k
