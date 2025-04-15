awk -F "," '$3>1 && $6<0.05' ../Lymphoid3up_vs_Lymphoid_12down_all_gene.0528.csv |perl -e 'while(<>){chomp;s/\"//g;@a=split(/\,/);print "$a[0]\tbulk_diff_lym3\n"}' > all.for.enrich.gene.list

