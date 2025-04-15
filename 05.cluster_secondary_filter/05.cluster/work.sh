perl -e '@x=(2000);@y=(50);@z=(0.3,0.5,0.6,0.8,1,1.2,1.5,1.8,2,2.5,3);while(<>){chomp;$t=$_;$o=`cat cluster.py`;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){`mkdir -p $t/$x[$i]_$y[$j]_$z[$k]`;`ln -s /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/04.data.filter_round3/$t/pancancer.ref.0804.final.$t.umap.h5ad $t/$x[$i]_$y[$j]_$z[$k]`}}}}' celltype.list1

perl -e '@x=(2000);@y=(50);@z=(0.3,0.5,0.6,0.8,1,1.2,1.5,1.8,2,2.5,3);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){$o=`cat cluster.py`;$p="$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].py";open OUT,">$p";$o=~s/aaaa/$t/g;$o=~s/bbbb/$x[$i]/g;$o=~s/cccc/$y[$j]/g;$o=~s/dddd/$z[$k]/g;print OUT "$o";$p="$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].sh";open OUT,">$p";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/05.cluster/$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].py\n"}}}}' celltype.list1

perl -e '@x=(2000);@y=(50);@z=(0.3,0.5,0.6,0.8,1,1.2,1.5,1.8,2,2.5,3);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){print "cd $t/$x[$i]_$y[$j]_$z[$k]\nqsub -cwd -l vf=50G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q cluster.$x[$i]_$y[$j]_$z[$k].sh\ncd -\n"}}}}' celltype.list1 |sh

ls */*/core*|perl -e 'while(<>){chomp;@a=split(/\//);print "cd $a[0]/$a[1]\nrm -r figures \*sh.* core*\nqsub -cwd -l vf=100G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q  cluster.$a[1].sh\ncd -\n"}'|sh

perl -e '@x=(2000);@y=(50);@z=(0.3,0.5,0.6,0.8,1,1.2,1.5,1.8,2,2.5,3);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){print "cp $t/$x[$i]_$y[$j]_$z[$k]/figures/umap.on.cluster.png sub_CD4CD8NK_cluster/$t.$z[$k].umap.on.cluster.png\ncp $t/$x[$i]_$y[$j]_$z[$k]/figures/umap.Phenotype.cluster.png sub_CD4CD8NK_cluster/$t.$z[$k].umap.Phenotype.cluster.png\ncp $t/$x[$i]_$y[$j]_$z[$k]/figures/umap.Tissue.cluster.png sub_CD4CD8NK_cluster/$t.$z[$k].umap.Tissue.cluster.png\n"}}}}' celltype.list1 |sh &

for i in *pca/*pca;do for j in $i/*;do echo -e "mkdir -p sub_CD4CD8_pcatest_cluster/`basename $i`/`basename $j`\ncp $j/figures/*png sub_CD4CD8_pcatest_cluster/`basename $i`/`basename $j`";done;done|sh &

for i in NK_pca/NK_pca/;do for j in $i/*;do echo -e "mkdir -p sub_NK_pcatest_cluster/`basename $i`/`basename $j`\ncp $j/figures/*png sub_NK_pcatest_cluster/`basename $i`/`basename $j`";done;done|sh &

