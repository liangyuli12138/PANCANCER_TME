ls *Lymphoid3*|perl -e 'while(<>){chomp;$id=$1 if(/diff\.(\S+)\_Lym/);open IN,$_;<IN>;$n=0;while(<IN>){chomp;$n++;if($n<=3000){print "$_,$id\n"}}}' > lym3up.diff.gene.csv

ls *Lymphoid3*|perl -e 'while(<>){chomp;$id=$1 if(/diff\.(\S+)\_Lym/);open IN,$_;<IN>;$n=0;while(<IN>){chomp;if(/\,RPS/ || /,RPL/){next};$n++;if($n<=100){print "$_,$id\n"}}}'|cut -d "," -f 2  > lym3up100.diff.gene.csv.list

