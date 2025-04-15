cut -f 1 sn.list |while read i;do echo perl get.cell.group.pl ../../../04.plot/at.group.zscore.r1/$i.group.at ../../../04.plot/at/$i.ori.at \> cell.group/$i.cell.group.csv;done|sh &


