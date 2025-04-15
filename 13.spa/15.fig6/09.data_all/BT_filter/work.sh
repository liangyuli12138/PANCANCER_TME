perl -e 'print "cell,difftype\n";while(<>){chomp;@a=split(/,/);if($a[-7] =~/Lymphoid_CD4/){$o="$a[-2].CD4_T";print "$a[0],$o\n"};if($a[-7] =~/Lymphoid_B/){$o="$a[-2].B";print "$a[0],$o\n"}}' pancancer.icar.background.cell.obs > bt.at


