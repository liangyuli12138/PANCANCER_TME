perl -e 'while(<>){chomp;@a=split(/,/);if($a[0]=~/(\S+)\_(\d+\_\d+)$/){$x=$1;$y=$2;$ha{$x}.="$y,$a[1],$a[2]\n"}};for $i(keys %ha){$o="at/$i.tls.at";open OUT,">$o";print  OUT "cell,celltype,cluster\n$ha{$i}"}' immune.list.at 

