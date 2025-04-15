perl -e 'while(<>){chomp;@a=split(/\t/);for ($i=2;$i<@a;$i++){print "$a[0],$a[$i]\n"}}' go.list >> hallmarker.list 

