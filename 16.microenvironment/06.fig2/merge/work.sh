ls ../out/*UCell.score.csv|perl -e '$h=`cat head`;print "$h";while(<>){chomp;$id=$1 if(/out\/(\S+)\.UCell/);open IN,$_;<IN>;while(<IN>){chomp;s/\"/\"$id\_/;@a=split(/\,/);print "$a[0]";for($i=4;$i<@a;$i++){print ",$a[$i]"};print "\n"}}' > merge.ucell.score.ori.csv &

