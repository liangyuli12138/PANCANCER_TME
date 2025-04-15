cut -f 1 sn.list |while read i;do echo perl get.at.pl colour.list ../04.plot/at/$i.all.at \> at/$i.all.at;done|sh

