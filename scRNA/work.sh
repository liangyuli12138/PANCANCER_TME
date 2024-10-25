
perl -e 'while(<>){chomp;@a=split;print "mkdir -p $a[0]/$a[1]\n"}' donor.list.c |sh

perl -e 'while(<>){chomp;@a=split;print "ln -s /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancancer/08.filter_by_gene/07.cluster_sub/01.extract/donor/plot_clustering.$a[0].$a[1].remarker.h5ad $a[0]/$a[1]\n"}' donor.list.c |sh


perl -e 'while(<>){chomp;@t=split;$f=`cat run_infercnvpy.py`;$f=~s/aaaa/$t[0]/g;$f=~s/bbbb/$t[1]/g;open OUT,">$t[0]/$t[1]/run_infercnvpy.py";print OUT "$f";open OUT,">$t[0]/$t[1]/run_infercnvpy.$t[1].sh";print OUT "python run_infercnvpy.py\n"}' donor.list.c

perl -e 'while(<>){chomp;@a=split;print "ln -s /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancancer/05.merge_ori_maxtir/05.cancer_cnv/01.find_cancer/genes.ref.gtf $a[0]/$a[1]\n"}' donor.list.c |sh

perl -e 'while(<>){chomp;@a=split;print "cd $a[0]\/$a[1]\nqsub -cwd -l vf=30G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q run_infercnvpy.$a[1].sh\ncd -\n"}' donor.list.c |sh

perl -e 'while(<>){chomp;@a=split;print "cp -r $a[0]/$a[1]/figures find_cancer_single_donor_0506/figures.$a[0].$a[1]\n"}' donor.list.c |sh

