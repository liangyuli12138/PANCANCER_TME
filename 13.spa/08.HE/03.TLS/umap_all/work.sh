cp ../at/*input at/

for i in ../at/*at;do echo perl get.at.pl $i;done|sh

perl -e 'while(<>){chomp;@a=split;$o=`cat CellBin_RCTD_plot.hg.R`;$o=~s/aaaa/$a[1].$a[0]/g;$o=~s/bbbb/$a[0]/g;open OUT,">shell/plot/plot.$a[0].R";print OUT "$o";open OUT,">shell/plot/plot.$a[0].sh";print OUT "/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/umap_all/shell/plot/plot.$a[0].R\n";}' sn.list

perl -e 'while(<>){chomp;@a=split;print "perl get.cell.pl merge.list /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/region/$a[0]_immue_region.csv /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$a[0]_cellbin.final.celltype.obs.csv at_cell/$a[0].cell at_cell/$a[0].cell.merge\n"}' sn.list |sh

perl -e 'while(<>){chomp;@a=split;$o=`cat CellBin_RCTD_plot.hg.R`;$o=~s/aaaa/$a[1].$a[0]/g;$o=~s/bbbb/$a[0]/g;open OUT,">shell/plot_cell/plot.$a[0].R";print OUT "$o";open OUT,">shell/plot_cell/plot.$a[0].sh";print OUT "/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/umap_all/shell/plot_cell/plot.$a[0].R\n";}' sn.list

perl -e 'while(<>){chomp;@a=split;$o=`cat CellBin_RCTD_plot.hg.R`;$o=~s/aaaa/$a[1].$a[0]/g;$o=~s/bbbb/$a[0]/g;open OUT,">shell/plot_merge/plot.$a[0].R";print OUT "$o";open OUT,">shell/plot_merge/plot.$a[0].sh";print OUT "/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/umap_all/shell/plot_merge/plot.$a[0].R\n";}' sn.list

