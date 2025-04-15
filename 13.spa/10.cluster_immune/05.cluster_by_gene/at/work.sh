perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$hc{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;$id=$1 if(/at\/(\S+)\.ori/);open IN,$_;while(<IN>){chomp;s/\"//g;@a=split(/,/);if(exists $hc{$a[2]}){$ip=$id."_".$a[1];$ha{$ip}="$a[2],$a[3]"}}};open IN2,$ARGV[2];<IN2>;print "cell\n";while(<IN2>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$a[0]\n"}}' immune.list ori.list merge_25chip_immune_area.obs > immune.list.input

#perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$hc{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;$id=$1 if(/at\/(\S+)\.ori/);open IN,$_;while(<IN>){chomp;s/\"//g;@a=split(/,/);if(exists $hc{$a[2]}){$ip=$id."_".$a[1];$ha{$ip}="$a[2],$a[3]"}}};open IN2,$ARGV[2];<IN2>;print "cell,celltype,cluster\n";while(<IN2>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$a[0],$ha{$a[0]}\n"}}' immune.list ori.list merge_25chip_immune_area.obs > immune.list.at
perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$hc{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;$id=$1 if(/at\/(\S+)\.ori/);open IN,$_;while(<IN>){chomp;s/\"//g;@a=split(/,/);if(exists $hc{$a[2]}){$ip=$id."_".$a[1];$ha{$ip}="$a[2],$a[3],$id"}}};open IN2,$ARGV[2];<IN2>;print "cell,celltype,cluster,batch\n";while(<IN2>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$a[0],$ha{$a[0]}\n"}}' immune.list ori.list merge_25chip_immune_area.obs > immune.list.at

perl -e '<>;while(<>){chomp;@a=split(/,/);if($a[0]=~/(\S+)\_\d+\_\d+/){$id=$1};$a[0]=~s/$id\_//;$ha{$id}.="$a[0]\n"};for $i(keys %ha){open OUT,">at/$i.input";print OUT "cell\n$ha{$i}"}' immune.list.at
perl -e '<>;while(<>){chomp;@a=split(/,/);if($a[0]=~/(\S+)\_\d+\_\d+/){$id=$1};$a[0]=~s/$id\_//;$ha{$id}.="$a[0],$a[1],$a[2]\n"};for $i(keys %ha){open OUT,">at/$i.at";print OUT "cell,celltype,tlsclu\n$ha{$i}"}' immune.list.at 

perl get.at.pl immune.list.at obs.list 

