perl merge.pl merge.list list.filter pancancer.icar.all.cell.obs ori.filter.list 

ls /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/05.distance/at/*at > all.list

perl merge.pl merge.list list.filter pancancer.icar.all.cell.obs ori.filter.list

perl -e 'while(<>){chomp;$a=$_;@x=split(/\//);print "perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/05.distance/get.all.dis.pl $a $a \> /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/05.distance/out/$x[-1].out\n"}' all.list  > all.sh

cat out/*Lymphoid3*|perl -e 'while(<>){chomp;@a=split(/,/);print "$a[-3],$a[0],$a[-5],$a[1],$a[2]\n"}' > Lymphoid3.distance.stat.csv
cat out/*1_2*|perl -e 'while(<>){chomp;@a=split(/,/);print "$a[-3],$a[0],$a[-5],$a[1],$a[2]\n"}' > Lymphoid_1_2.distance.stat.csv

