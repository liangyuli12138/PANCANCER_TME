perl -e 'while(<>){chomp;@a=split;print "perl get.new.at.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/04.plot/at/$a[0].all.at /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$a[0]","_region.csv > at/$a[0].all.at\n"}' all.sn.list |sh &

cut -f 1 all.sn.list |while read i;do echo perl get.color.at.pl colour.list at/$i.all.at \> at_change_colour/$i.all.at;done|sh &

