perl -e 'for($i=1;$i<=3;$i++){$f="icm$i.diff.csv";open IN,$f;<IN>;$n=0;while(<IN>){chomp;@a=split(/,/);if($n<200 && $a[2]>1){print "ICM$i,$a[1]\n";$n++}}}' > diff.gene.list
