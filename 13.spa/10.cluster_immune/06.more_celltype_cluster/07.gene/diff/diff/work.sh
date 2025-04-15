for i in diff.*csv;do echo perl sort.pl $i \> filter/$i.filter;done|sh

