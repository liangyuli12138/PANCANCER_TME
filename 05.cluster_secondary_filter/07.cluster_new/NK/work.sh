perl -e 'print "cell,groups_ref\n";open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[-1] eq "10"){$ha{$a[0]}=1}};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);if($a[-1] eq "0" || $a[-1] eq "4" || $a[-1] eq "12"){print "$a[0],Lymphoid_MAIT\n"}elsif($a[-1] eq "6" || $a[-1] eq "1" || $a[-1] eq "9" || $a[-1] eq "13" || $a[-1] eq "5"){print "$a[0],Lymphoid_NK_CD56+\n"}elsif($a[-1] eq "7" || exists $ha{$a[0]} || $a[-1] eq "10" || $a[-1] eq "8"){print "$a[0],Lymphoid_NKT\n"}elsif(($a[-1] eq "3" || $a[-1] eq "2") && !exists $ha{$a[0]}){print "$a[0],Lymphoid_NK_CD16+\n"}elsif($a[-1] eq "11"){print "$a[0],Lymphoid_ILC\n"}}' pancancer.ref.0807.final.NK.umap.obs.txt.2.5 pancancer.ref.0807.final.NK.umap.obs.txt.1.5 > NK.at


perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$b,"}' cell.list |les

perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$aa\n\n"}' cell.list >> cluster.py 

