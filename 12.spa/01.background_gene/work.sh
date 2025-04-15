les gem.list |while read i;do echo perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/01.background_gene/stat.bg.pl $i \> /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/01.background_gene/filter.gem/`basename $i`.filter;done > all.sh

for i in filter.gem/*diff;do echo perl get.filter.gene.pl $i \> $i.rm;done|sh &

les sn.list |while read i;do echo perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/01.background_gene/filter.gem.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/01.background_gene/filter.gem/$i.cellbin.gef.gem.filter.diff.rm `ls /jdfssz1/ST_TSCBI/P22Z10200N0433/USER/wubin2/project/pancancer/09.spatial/00.data_0517/*/$i.cellbin.gef.gem` \> `ls /jdfssz1/ST_TSCBI/P22Z10200N0433/USER/wubin2/project/pancancer/09.spatial/00.data_0517/*/$i.cellbin.gef.gem`.gem;done > all.filter.sh

