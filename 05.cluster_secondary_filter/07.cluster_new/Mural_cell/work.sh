perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;if(/delete/){}else{@a=split;$ha{$a[0]}=$a[1]}};print "cell,groups_ref\n";open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[-1]}){print "$a[0],$ha{$a[-1]}\n"}}' annotation.list pancancer.ref.0723.final.Mural_cell.umap.obs.txt > Mural_cell.at

perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$b,"}' cell.list 

perl -e '$tt=`cat tmp1.py`;while(<>){chomp;@a=split(/\t/);$a=$a[0];$b=$a;$aa=$tt;$b=~s/\//_/g;$b=~s/\s/_/g;$b=~s/\+/__/g;$b=~s/\-/_/g;;$aa=~s/aaaa/$a/g;$aa=~s/bbbb/$b/g;$n++;if($n%3==1){$b="[ax$b"}elsif($n%3==0){$b="ax$b]"}else{$b="ax$b"};print "$aa\n\n"}' cell.list 

perl -e 'while(<>){chomp;print "\"$_\","}' marker.list 


