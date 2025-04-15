perl -e 'print "cell,groups_refs\n";<>;while(<>){chomp;s/TAM/Marco/g;s/RTM/Marco/g;@a=split(/,/);if($a[-2] eq "13"){print "$a[0],Myeloid_pDC\n"}elsif($a[-2] eq "8"){if($a[-1] eq "Myeloid_cDC1"){print "$a[0],Myeloid_cDC1\n"}elsif($a[-1] eq "Myeloid_cDC3"){print "$a[0],Myeloid_cDC3\n"}else{next}}elsif($a[-2] eq "11"){print "$a[0],Myeloid_Marco_SPP1\n"}elsif($a[-2] eq "9"){print "$a[0],Myeloid_Marco_C1QC\n"}elsif($a[-2] eq "6"){print "$a[0],Myeloid_Marco_LYVE1\n"}elsif($a[-2] eq "10" || $a[-2] eq "0"){print "$a[0],Myeloid_Mono\n"}else{print "$a[0],$a[-1]\n"}}' pancancer.ref.final.final.Myeloid.umap.old.obs.csv > Myeloid.at

perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$b,"}' cell.list

perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$aa\n\n"}' cell.list

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=1};open IN1,$ARGV[1];while(<IN1>){chomp;if(exists $ha{$_}){}else{print "$_\n"}}' pancancer.final.0314.var.csv.input.add marker.list 

