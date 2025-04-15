perl -e 'print "\"\",\"id\",\"celltype\",\"max_score\",\"index\"\n";<>;while(<>){chomp;s/\"POLYGON.+\)\)\"\,//;@a=split(/,/);print "$a[0],\"$a[3]","_","$a[4]\",$a[-1],1,1\n"}' Merge_Live.B02324E3_meta_ori.csv > Merge_Live.B02324E3_meta_ori.monify.csv.list

