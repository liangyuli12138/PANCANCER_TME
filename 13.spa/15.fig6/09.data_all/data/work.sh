perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/\_Lymphoid\S+//;$ha{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/\,/);if(/cell/){print "$_\n"}else{if(exists $ha{$a[0]}){print "$_\n"}}}' filter.cluster.list new.group.list.at > new.group.list.at.filter

cut -f 1 sn.list |while read i;do echo perl get.group.pl new.group.list.at.filter at100/$i.ex.at $i \> at100.group/$i.ex.at;done|sh &

