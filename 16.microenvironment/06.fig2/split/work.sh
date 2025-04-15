perl -e '$t=<>;while(<>){chomp;@a=split;$a[0]=~s/\"//g;if($a[0]=~/(\S+)\_(\d+\_\d+)/){$x=$1;$y=$2};s/$x\_//;$ha{$x}.="$_\n"};for $i(keys %ha){$o="$i.UCell.score.csv";open OUT,">split/$o";print OUT "$t","$ha{$i}"}'  merge.ucell.score.scale.csv &

perl -e 'while(<>){chomp;@a=split;$t="$a[1].$a[0]";$o=`cat temp.R`;$o=~s/aaaa/$a[0]/g;$o=~s/bbbb/$t/g;open OUT,">shell/$t.R";print OUT "$o";open OUT,">shell/$t.sh";print OUT "/hwfssz4/BC_PUB/Software/07.User-defined/03.Animal_Plant/wubin/mamba/bin/Rscript $t.R\n"}' sn.list

perl -e 'while(<>){chomp;print "mkdir fig2_ucell_REACTOME/$_\ncp out/*$_* fig2_ucell_REACTOME/$_\n"}' filter.list |sh &

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$ha{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;for $i(keys %ha){print "cp fig2_ucell_REACTOME/$_/*$i* fig2_ucell_REACTOME.filter/$_\n"}}' sample.list filter.list |sh

