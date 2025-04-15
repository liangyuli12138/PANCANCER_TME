perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;if(/delete/){}else{@a=split;$ha{$a[0]}=$a[1]}};print "cell,groups_ref\n";open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[-1]}){print "$a[0],$ha{$a[-1]}\n"}}' annotation.list pancancer.ref.0728.final.B.umap.obs.txt > B.at

#perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=1};open IN1,$ARGV[1];while(<IN1>){chomp;if(exists $ha{$_}){print "$_\n"}}' /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/02.cluster/B/r.0.5/pancancer.ref.0728.final.B.umap.var.txt marker.list |perl -e 'while(<>){chomp;print "\"$_\","}'

perl -e 'open IN2,$ARGV[2];while(<IN2>){chomp;@a=split(/,/);if($a[-1] eq "4"){$ha{$a[0]}=1}};open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[-1] eq "3" || $a[-1] eq "20" || $a[-1] eq "14" || $a[-1] eq "9" || $a[-1] eq "19" || $a[-1] eq "24" || $a[-1] eq "30" || $a[-1] eq "32" || $a[-1] eq "21" || $a[-1] eq "2" || $a[-1] eq "7"){$ha{$a[0]}=1}};open IN1,$ARGV[1];<IN1>;print "cell,groups_refs\n";while(<IN1>){chomp;@a=split(/,/);if(/Plamsa/){if(exists $ha{$a[0]}){print "$a[0],Lymphoid_Plamsa_IGLC\n"}else{print "$a[0],Lymphoid_Plamsa_IGKC\n"}}else{print "$_\n"}}' pancancer.ref.final.final.pla.umap.obs.csv B.at.ori pancancer.ref.0807.final.B.umap.obs.txt > B.at


perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$b,"}' cell.list
perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$aa\n\n"}' cell.list

perl -e 'while(<>){chomp;@a=split(/\t/);for($i=0;$i<@a;$i++){print "\"$a[$i]\","}}' marker.list 

