perl -e 'print "cell,groups_ref\n";open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[-1] eq "0"){print "$a[0],Lymphoid_CD4_Treg\n"}elsif($a[-1] eq "1"){print "$a[0],Lymphoid_CD4_Tn\n"}elsif($a[-1] eq "2"){print "$a[0],Lymphoid_CD4_CTL\n"}elsif($a[-1] eq "3"){print "$a[0],Lymphoid_CD4_Tcm\n"}elsif($a[-1] eq "4"){print "$a[0],Lymphoid_CD4_Tstr\n"}elsif($a[-1] eq "5"){print "$a[0],Lymphoid_CD4_Tfh\n"}elsif($a[-1] eq "6"){print "$a[0],Lymphoid_CD4_Th17\n"}}' pancancer.ref.0807.final.CD4_pca.umap.obs.txt > CD4.at

perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$b,"}' cell.list 

