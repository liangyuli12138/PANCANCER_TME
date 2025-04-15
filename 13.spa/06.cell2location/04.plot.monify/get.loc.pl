$n=-1;while(<>){chomp;@a=split(/\s/);$n++;for($i=0;$i<@a;$i++){if($a[$i]==255){print "$a[$i]\t$i\t$n\n"}}}
