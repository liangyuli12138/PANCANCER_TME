ls diff.*csv|perl -e 'while(<>){chomp;$t=$_;open IN,$_;<IN>;while(<IN>){chomp;@a=split(/,/);print "$t,$a[1],$a[-1]\n"}};' > all.cell.pct.xls

