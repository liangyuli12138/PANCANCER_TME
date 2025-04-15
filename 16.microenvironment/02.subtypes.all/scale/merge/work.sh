ls ../../out/*UCell.score.csv|perl -e '$h=`cat head`;print "$h";while(<>){chomp;$id=$1 if(/out\/(\S+)\.UCell/);open IN,$_;<IN>;while(<IN>){chomp;s/\"/\"$id\_/;@a=split(/\,/);print "$a[0],$a[4],$a[5],$a[6],$a[7]\n"}}' > hallmark.30.score.ori.csv &

sh scale.sh &

