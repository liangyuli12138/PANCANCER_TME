awk -F "," '$4>0.5 && $7>0.2 {print $2}' diff.pancancer.icar.mig.cell.h5ad.Lymphoid3.csv > diff.pancancer.icar.mig.cell.h5ad.Lymphoid3.csv.gene.csv

