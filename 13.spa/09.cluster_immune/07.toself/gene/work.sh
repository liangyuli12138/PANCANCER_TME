perl -e 'while(<>){chomp;@a=split(/\t/);$t="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/at.pla/$a[0].ori.at";open OUT,">at/$a[0].input";open OUT1,">at/$a[0].at";print OUT "cell\n";print OUT1 "cell,tlsclu\n";open IN,$t;<IN>;while(<IN>){chomp;if(/Cluster/){s/\"//g;@b=split(/,/);print OUT "$b[1]\n";print OUT1 "$b[1],$b[3]\n";}}}' sn.list


perl -e 'while(<>){chomp;@a=split;$o=`cat get.mean.py`;$o=~s/aaaa/$a[0]/g;open OUT,">shell/$a[0].py";print OUT "$o";open OUT,">shell/$a[0].sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].py\n"}' sn.list
(base) wubin2 /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/gene


perl -e 'print ",CD19,MS4A1,CD3D,CXCL9,CXCL10,CXCL11,CXCL12,CXCL13,CXCR5,ICOS,CTLA4,PDCD1,LTA,IL21,IL6,IL17A,FCER2\n";while(<>){chomp;@a=split(/\t/);$t="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/gene/out/$a[0].geneAmean.csv";open IN,$t;<IN>;while(<IN>){chomp;print "$a[0]","$_\n"}} ' sn.list > merge.gene.csv

#perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$_};open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);print "$a[0],$ha{$a[0]}\n"}' merge.gene.csv ../colocation/immune.cluster.r0.5.obs > merge.gene.csv.at
perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$_};open IN1,$ARGV[1];<IN1>;print "cell,CD19,MS4A1,CD3D,CXCL9,CXCL10,CXCL11,CXCL12,CXCL13,CXCR5,ICOS,CTLA4,PDCD1,LTA,IL21,IL6,IL17A,FCER2\n";while(<IN1>){chomp;@a=split(/,/);print "$ha{$a[0]}\n"}' merge.gene.csv ../colocation/immune.cluster.r0.5.obs > merge.gene.csv.at

perl -e 'while(<>){chomp;print "sc.pl.umap(adata_concat,color=\"$_\",frameon=False, na_color=\"grey\",save=\".$_.gene.png\")\n"}' gene.list >> gene.py


