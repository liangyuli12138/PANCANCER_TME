perl -e 'while(<>){chomp;if(/^\#/ || /^geneID/){print "$_\n"}else{@a=split;$x=26459-$a[2];$x=$x<26457?$x:26457;$y=$a[1];$y=$y<26457?$y:26457;print "$a[0]\t$x\t$y\t$a[3]\t$a[4]\n"}}' B02324E3.gem.old > B02324E3.gem &

perl -e 'while(<>){@a=split;if($a[1]>$x){$x=$a[1]};if($a[2]>$y){$y=$a[2]}};print "$x\t$y\n"' B02324E3.gem.old
26459	26457


gzip B02324E3.gem &

perl -e 'print "cellbin\n";<>;while(<>){chomp;s/\"//g;@a=split(/,/);print "$a[0]\n"}' ../result/B02324E3/B02324E3_cellbin.filter.gene.list > ../result/B02324E3/B02324E3_cellbin.filter.gene.list.cellbin

