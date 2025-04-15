perl -e 'print ",sample,variable,score\n";$t=<>;chomp $t;@t=split(/\t/,$t);while(<>){chomp;@a=split(/\t/);for($i=1;$i<4;$i++){$n++;print "\"$n\",T$a[0],$t[$i],$a[$i]\n"}}' ori.xls > t.ibp.csv

