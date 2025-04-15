perl -e 'print "cell,groups_ref\n";open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[-1] eq "0"){print "$a[0],Lymphoid_CD8_Tn\n"}elsif($a[-1] eq "1" || $a[-1] eq "2" || $a[-1] eq "4" || $a[-1] eq "6" || $a[-1] eq "11" || $a[-1] eq "15" || $a[-1] eq "10"){print "$a[0],Lymphoid_CD8_Teff\n"}elsif($a[-1] eq "3" || $a[-1] eq "8" || $a[-1] eq "9"){print "$a[0],Lymphoid_CD8_Tex\n"}elsif($a[-1] eq "7"){print "$a[0],Lymphoid_CD8_Tisg\n"}elsif($a[-1] eq "5" || $a[-1] eq "13" || $a[-1] eq "16"){print "$a[0],Lymphoid_CD8_Tm\n"}elsif($a[-1] eq "12"){print "$a[0],Lymphoid_CD8_Tstr\n"}}' pancancer.ref.0807.final.CD8_pca.umap.obs.txt > CD8.at

perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$b,"}' cell.list 

perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$aa\n\n"}' cell.list 

