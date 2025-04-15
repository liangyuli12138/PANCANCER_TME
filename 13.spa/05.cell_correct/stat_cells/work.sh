perl -e '$o=`cat head.R`;print "$o";while(<>){chomp;@a=split;$s=$a[0];$o=`cat tmp.R`;$o=~s/outfix_sample/$s/g;print "$o"}' sn.list > stat.R
vi stat.R

perl -e 'while(<>){chomp;$id=$_;<>;$a=<>;$b=<>;$a=~s/\s+/\t/g;$n=$1 if($b=~/(\d+)/);print "$id\t$a\t$n\n"}' stats.filter.1030.xls > stats.filter.1030.txt

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;s/_nFeature_Spatial//;print "$ha{$_}\t$_\t";<IN1>;$t=<IN1>;@t=split(/\s+/,$t);print "$t[3]\t";<IN1>;<IN1>;$t=<IN1>;@t=split(/\s+/,$t);print "$t[3]\t";<IN1>;$x=<IN1>;$x=~s/^\s+//;chomp $x;print "$x\n"}' sn.list stats.filter.1217.xls > s1a.pancancer.ST.median.count.xls

ls /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/*/*final.obs|perl -e 'while(<>){chomp;open IN,$_;while(<IN>){chomp;@a=split(/,/);$n+=$a[-5];$m++}};$o=$n/$m;print "$o\n"'
850.254262210871

