perl -e 'while(<>){chomp;@a=split;$o=`cat get.mean.py`;$o=~s/aaaa/$a[0]/g;open OUT,">shell/$a[0].py";print OUT "$o";open OUT,">shell/$a[0].sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].py\n"}' sn.list

perl -e 'print ",CD19,MS4A1,CD3D,CXCL9,CXCL10,CXCL11,CXCL12,CXCL13,CXCR5,ICOS,CTLA4,PDCD1,LTA,IL21,IL6,IL17A,FCER2\n";while(<>){chomp;@a=split(/\t/);$t="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/gene.normalizing/out/$a[0].geneAmean.csv";open IN,$t;<IN>;while(<IN>){chomp;print "$a[0]","$_\n"}} ' sn.list > merge.gene.csv

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$_};open IN1,$ARGV[1];<IN1>;print "cell,CD19,MS4A1,CD3D,CXCL9,CXCL10,CXCL11,CXCL12,CXCL13,CXCR5,ICOS,CTLA4,PDCD1,LTA,IL21,IL6,IL17A,FCER2\n";while(<IN1>){chomp;@a=split(/,/);print "$ha{$a[0]}\n"}' merge.gene.csv ../colocation/immune.cluster.r0.5.obs > merge.gene.csv.at

