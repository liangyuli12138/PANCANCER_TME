les sample_list |while read i;do echo /ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/cellbin_122/bin/python3 /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/05.cell_correct/cell_correct.cpu.py -m /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/04.filter_bright/02.filter/tif_file/$i""_mask.bright.tif -g /jdfssz1/ST_TSCBI/P22Z10200N0433/USER/wubin2/project/pancancer/09.spatial/00.data_0517/gem.file/$i.gem.gz -o  /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/$i/ -d 15;done > run.cpu.sh

les sample_list |while read i;do echo cp result/$i/cell_mask_corect.compression.tif filter.brigth.tif/$i.cell_mask_corect.tif ;done|sh

#les sample_list |while read i;do echo -e "cut -f 1-5,8 result/$i/cell_mask_corect.gem > result/$i/cell_mask_corect.gem.tmp\nmv result/$i/cell_mask_corect.gem.tmp result/$i/cell_mask_corect.gem";done|sh &

