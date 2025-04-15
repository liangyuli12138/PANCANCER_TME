cat DESeq2_res1_down_FC2_Lym3_vs_12.csv DESeq2_res1_up_FC1_Lym3_vs_12.csv |perl -e 'while(<>){chomp;@a=split(/,/);if(/\"\"/){next};print "$a[0]\n"}' > merge.deseq2.list

