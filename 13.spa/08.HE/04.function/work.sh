perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$b=$a[0];$b=~s/\_TLS\S+//;$ha{$a[0]}=$a[1];$hb{$b}=1};open IN1,$ARGV[1];while(<IN1>){chomp;$id=$1 if(/at\/(\S+)\.at/);if(exists $hb{$id}){open IN,$_;open OUT0,">at/$id.at";open OUT1,">at/$id.input";print OUT1 "cell\n";print OUT0 "cell,cluster_TLS,group_TLS\n";while(<IN>){chomp;@c=split(/,/);if(exists $ha{$c[1]}){print OUT0 "$_,$ha{$c[1]}\n";print  OUT1 "$c[0]\n"}}}}' cluster.list at.list

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$b=$a[0];$b=~s/\_TLS\S+//;$hb{$b}=1};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;if(exists $hb{$a[0]}){print "$_\n"}}' cluster.list /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/08.marker/sn.list > sn.list

perl -e 'while(<>){chomp;@a=split;$o=`cat tmp.py`;$o=~s/aaaa/$a[0]/g; print "$o"}' sn.list >> diff.py 

echo /jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/diff.py TLS_C3 /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/out/TLS_C3.diff.csv > shell/s3.sh

