cp /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/00.data/at.icm/* at &
for i in *at;do echo split -l 5000 -d -a 3 $i $i.;done|sh &
ls /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/03.dist/at/*all.at.* > all.list

perl -e 'while(<>){chomp;$a=$_;$b=$_;$a=~s/all.at.\S+/all.at/;@x=split(/\//);print "perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/03.dist/get.dist.pl $b $a \> /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/03.dist/out/$x[-1].out\n"}' all.list  > all.sh
perl -e 'while(<>){chomp;$a=$_;$b=$_;$a=~s/all.at.\S+/all.at/;@x=split(/\//);print "perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/03.dist/get.dist.pl $b $a \> /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/03.dist/out/$x[-1].out\n"}' all.list  > all.sh

perl filter.pl all.sn.list all.list

perl filter.all.pl all.sn.list all.list &

