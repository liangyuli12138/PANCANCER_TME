perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);for($i=0;$i<@a;$i++){$hb{$a[$i]}=1}};while(<>){chomp;if(/\[(.+)\],$/){$a=$1;$a=~s/,$//;@b=split(/,/,$a);for($i=0;$i<@b;$i++){if(!exists $ha{$b[$i]} && !exists $hb{$b[$i]}){print "$b[$i],";$ha{$b[$i]}=1}};print "\n"}}' out.list plot.py |les 

