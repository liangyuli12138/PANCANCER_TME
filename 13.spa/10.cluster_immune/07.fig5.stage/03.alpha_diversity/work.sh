perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}="$a[-1]"};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/\,/);if(exists $ha{$a[0]}){print "$_,$ha{$a[0]}\n"}}' otu_table.csv l2m.stat.merge.csv > l2m.stat.merge.csv.sn

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];$t=<IN1>;print "$t";while(<IN1>){chomp;@a=split(/\,/);s/$a[0]/$ha{$a[0]}/;print "$_\n"}' l2m.stat.merge.csv.sn cell.merge.at > cell.merge.at.group

grep Lymphoid l2m.stat.merge.csv.sn.ori >  l2m.stat.merge.csv.sn 

awk -F "," '$2 !~/0/' l2m.stat.merge.csv.sn > l2m.stat.merge.csv.sn.filter

