for i in diff/*csv;do echo head -n 3001 $i \> contour.sub.diff/`basename $i`;done|sh

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];<IN1>;print "cell,merge_celltype\n";while(<IN1>){chomp;@a=split(/,/);if($ha{$a[-14]} eq "delete" || !exists $ha{$a[-14]}){next};if(exists $ha{$a[-14]}){print "$a[0],$ha{$a[-14]}","_","$a[-3]\n"}}' list pancancer.icar.contour.cell.obs > caf.at

