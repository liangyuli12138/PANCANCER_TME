ls *list|perl -e 'while(<>){chomp;open IN,$_;$t="$_.pCD";open OUT,">$t";while(<IN>){if(/Migration/ || /max_score/ || /pDC/){print OUT "$_"}}}'

