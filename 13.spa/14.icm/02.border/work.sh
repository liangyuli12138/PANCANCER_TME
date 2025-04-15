cp ../01.plot/at/* at &
for i in *at;do echo split -l 5000 -d -a 3 $i $i.;done|sh &
ls /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/02.border/at/*all.at.* > all.list

perl -e 'while(<>){chomp;$a=$_;$b=$_;$a=~s/all.at.\S+/all.at/;@x=split(/\//);print "perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/02.border/get.dist.pl $b $a \> /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/02.border/out/$x[-1].out\n"}' all.list  > all.sh

perl -e 'for($i=0;$i<10000000;$i++){$c=`qstat`;@a=split(/\n/,$c);for($j=0;$j<@a;$j++){if($a[$j]=~/\s+t\s+/ || $a[$j]=~/\s+Eqw\s+/){@b=split(/\s+/,$a[$j]);$act="qdel $b[1]\nqsub -clear -cwd -l vf=0.3G,num_proc=1 -P P22Z10200N0431 -binding linear:1 -q st.q $b[3]"."h\n";`$act`}};sleep 600}' &

perl filter.pl all.sn.list all.list

perl filter.all.pl all.sn.list all.list &


