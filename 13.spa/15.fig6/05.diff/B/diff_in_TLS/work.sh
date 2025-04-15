head -n 101 diff.pancancer.icar.T.cell.h5ad.Lymphoid_1_2.csv |tail -n 100 |awk -F "," '{print "T.Lymphoid_1_2\t"$2}' > diff.T.csv
head -n 101 diff.pancancer.icar.T.cell.h5ad.Lymphoid3.csv |tail -n 100 |awk -F "," '{print "T.Lymphoid3\t"$2}' >> diff.T.csv

