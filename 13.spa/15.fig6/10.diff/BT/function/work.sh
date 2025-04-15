perl -e 'print "gene\tcluster\n";open IN0,$ARGV[0];while(<IN0>){chomp;open IN,$_;$id=$1 if(/diff\.IGH\/diff\.(\S+)\.csv/);<IN>;$n=0;while(<IN>){chomp;if(/,RPS/ || /,RPL/){next};if($n<200){@a=split(/\,/);print "$a[1]\t","$id\n";$n++}}}' diff1.list > all.for.enrich.gene.list

