perl -e '$t=<>;print "$t";while(<>){chomp;if(/F_\d_N/){print "$_\n"}}' pancancer.ref.umap.0723.obs.txt > pancancer.ref.umap.0723.ka.obs.txt

