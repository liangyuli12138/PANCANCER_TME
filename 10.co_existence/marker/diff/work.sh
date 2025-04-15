perl -e 'print "cell\n";open IN2,$ARGV[2];while(<IN2>){chomp;$hb{$_}=1};open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[1]}=$a[0]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[15]} && exists $hb{$a[8]}){print "$a[0],$ha{$a[15]}\n"}}' time.list pancancer.ref.0905.final.obs.csv sample.list > filter.at

perl -e 'for($i=1;$i<=5;$i++){print "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/10.co_existence/marker/diff/findmarker.py /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/10.co_existence/marker/diff/pancancer.ref.0905.final.h5ad TIME$i /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/10.co_existence/marker/diff/out/TIME$i.diff.csv\n"}' > all.diff.sh

