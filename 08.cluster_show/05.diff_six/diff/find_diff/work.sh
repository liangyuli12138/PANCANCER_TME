ls diff/*csv|perl -e 'print "ori_suq\tgene\tscore\tavg_log2FC\tp_val\tp_val_adj\tpct.1\tpct.2\tCluster\n";while(<>){chomp;@a=split(/\//);$id=$1 if(/\/(\S+)\.csv/);$t=`head -n 101 $_|tail -n 100 `;$t=~s/,/\t/g;@t=split(/\n/,$t);for($i=0;$i<@t;$i++){print "$t[$i]\t$id\n"}}' > merge.all.diffgene.xls.before.filter

cp merge.all.diffgene.xls.before.filter merge.all.diffgene.xls

cut -f 9 merge.all.diffgene.xls|sort -u |awk -F "." '!/Cluster/ {print $2}'|while read i;do echo mkdir -p find_plot/$i;done|sh

les cell.sort.list |awk '!/Mul/' |while read i;do echo perl get.py.pl $i;done|sh

les cell.sort.list |awk '!/Mul/'|cut -d "," -f 1|while read i;do echo -e "cd find_plot/$i\nqsub -cwd -l vf=50G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q find.$i.sh\ncd -";done|sh

