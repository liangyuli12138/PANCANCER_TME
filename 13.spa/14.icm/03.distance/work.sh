cp ../01.plot/at/* at &

perl -e 'while(<>){chomp;@a=split;print "perl merge.pl merge.list /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$a[0]","_cellbin.final.celltype.obs.csv at.ori/$a[0].all.at > at/$a[0].all.at\n"}' all.sn.list |sh &

for i in *at;do echo split -l 5000 -d -a 3 $i $i.;done|sh &
ls /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/03.distance/at/*all.at.* > all.list

perl -e 'while(<>){chomp;$a=$_;$b=$_;$a=~s/all.at.\S+/all.at/;@x=split(/\//);print "perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/03.distance/get.all.dis.pl $a $b \> /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/03.distance/out/$x[-1].out\n"}' all.list  > all.sh



