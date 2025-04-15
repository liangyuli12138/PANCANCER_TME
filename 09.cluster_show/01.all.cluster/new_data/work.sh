perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/\t/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];$t=<IN1>;chomp $t;print "$t,groups_sev\n";while(<IN1>){if(/Epithelium_Malig_Multi/){next};chomp;@a=split(/\,/);print "$_,$ha{$a[3]}\n"}' cell.groups_thi.list merge.all.all.level.filterB.six.at > merge.all.all.level.filterB.sev.at.old

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/\,/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];$t=<IN1>;chomp $t;print "$t\n";while(<IN1>){chomp;@a=split(/\,/);if(!/CD4/){print "$_\n"}elsif(exists $ha{$a[0]}){s/$a[1]/$ha{$a[0]}/;print "$_\n"}}' CD4.at merge.all.all.level.filterB.sev.at.old > merge.all.all.level.filterB.sev.at

