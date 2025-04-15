awk -F "," '$4>1 && $6<0.05 && $7>0.01 && $2 !~/^RP/ && $2 !~/^MT-/' diff.Lymphoid_1_2.csv > diff.Lymphoid_1_2.csv.filter

cat diff.Lymphoid_1_2.csv.filter diff.Lymphoid3.csv.filter |cut -d "," -f 2 > diff.Lymphoid.csv

