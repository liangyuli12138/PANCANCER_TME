perl -e 'while(<>){chomp;s/\t\t/\t/;s/\t+$//;@a=split(/\t/);undef %ha;for($i=1;$i<@a;$i++){if(!exists $ha{$a[$i]}){print "$a[0],$a[$i]\n"};$ha{$a[$i]}=1}}' gene.list > hallmarker.list

perl -e 'while(<>){chomp;@a=split;$t="$a[1].$a[0]";$o=`cat temp.R`;$o=~s/aaaa/$a[0]/g;$o=~s/bbbb/$t/g;open OUT,">shell/$t.R";print OUT "$o";open OUT,">shell/$t.sh";print OUT "/hwfssz4/BC_PUB/Software/07.User-defined/03.Animal_Plant/wubin/mamba/bin/Rscript $t.R\n"}' sn.list


