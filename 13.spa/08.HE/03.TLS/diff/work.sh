perl -e 'while(<>){chomp;@a=split;$x=$a[1].".".$a[0];print "mkdir -p out/$x\n"}' sn.list |sh

perl get.sh.pl sn.list 

