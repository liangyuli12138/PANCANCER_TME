sed 's/Plamsa_cell/B_cell/' ../celltype_merge/all.cell.dist.csv.type > all.cell.dist.csv.type &

perl -e '$t=<>;print "$t";while(<>){chomp;@a=split(/,/);if($a[3]=~/ICM/){$o="ICM"}else{$o="other"};print "$a[0],$a[1],$a[2],$a[3],$a[4],$a[3],$o\n"}' all.cell.dist.csv.type > all.cell.dist.csv.type.ICM

