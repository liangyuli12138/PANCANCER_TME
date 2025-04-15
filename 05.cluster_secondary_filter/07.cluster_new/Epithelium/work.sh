perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[2]}){$hb{$a[0]}=$ha{$a[2]}}};open IN2,$ARGV[2];while(<IN2>){chomp;@a=split(/,/);if($a[-1] eq "Epithelium" && exists $hb{$a[0]}){print "$a[0],$hb{$a[0]}\n"}}' mag.list all.monify.at pancancer.ref.0723.final.obs.csv > Epithelium.at


perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[2]}){$hb{$a[0]}=$ha{$a[2]}}else{if($a[-1]=~/Malig/){$hb{$a[0]}="Epithelium_Malig_Multi"}}};open IN2,$ARGV[2];while(<IN2>){chomp;@a=split(/,/);if($a[-1] eq "Epithelium" && exists $hb{$a[0]}){print "$a[0],$hb{$a[0]}\n"}}'  mag.list all.monify.at pancancer.ref.0723.final.obs.csv > Epithelium.all.at

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[2]}){$hb{$a[0]}=$ha{$a[2]}}else{if($a[-1]=~/Malig/){$hb{$a[0]}="Epithelium_Malig_Multi"}}};open IN2,$ARGV[2];while(<IN2>){chomp;@a=split(/,/);if($a[-1] eq "Epithelium" && exists $hb{$a[0]}){print "$a[0],$hb{$a[0]}\n"}elsif($a[-1] eq "Epithelium" && !exists $hb{$a[0]}){print "$a[0],Epithelium_TN\n"}}'  mag.list all.monify.at pancancer.ref.0723.final.obs.csv > Epithelium.all.all.at

