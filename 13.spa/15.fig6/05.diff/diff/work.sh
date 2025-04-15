awk -F "," '$4>1 && $6<0.05 && $7>0.01' diff.Lymphoid_1_2.csv > diff.Lymphoid_1_2.csv.filter
awk -F "," '$4>1 && $6<0.05 && $7>0.01' diff.Lymphoid3.csv > diff.Lymphoid3.csv.filter

awk -F "," '$4>0.5 && $6<0.05 && $7>0.01' diff.Lymphoid_1_2.csv|head -n 200 > diff.Lymphoid_1_2.csv.filter
awk -F "," '$4>0.5 && $6<0.05 && $7>0.01' diff.Lymphoid3.csv |head -n 200 > diff.Lymphoid3.csv.filter

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/\,+$//;@a=split(/,/);for($i=1;$i<@a;$i++){$ha{$a[$i]}=$a[0]}};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);$ha{$a[1]}=$ha{$a[1]}?$ha{$a[1]}:"NA";print "$_,$ha{$a[1]}\n"}' gopb_gene.csv diff.Lymphoid_1_2.csv.filter |awk '!/NA/'|wc -l
perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/\,+$//;@a=split(/,/);for($i=1;$i<@a;$i++){$ha{$a[$i]}=$a[0]}};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);$ha{$a[1]}=$ha{$a[1]}?$ha{$a[1]}:"NA";print "$_,$ha{$a[1]}\n"}' gopb_gene.csv diff.Lymphoid3.csv.filter |awk '!/NA/'|wc -l

cat diff.Lymphoid_1_2.csv.filter diff.Lymphoid3.csv.filter |cut -d "," -f 2 > diff.Lymphoid.filter

