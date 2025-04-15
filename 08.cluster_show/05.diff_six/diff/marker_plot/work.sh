perl -e 'while(<>){chomp;s/\t+/\t/g;s/\t$//;if(/^$/){next};print "$_\n"}' marker.list > marker.list.sort

perl -e 'while(<>){chomp;@a=split;for($i=1;$i<@a;$i++){if(!exists $ha{$a[0]}{$a[$i]}){$ha{$a[0]}{$a[$i]}=1;$hb{$a[0]}.="$a[$i],"}}};for $i(keys %hb){print "$i\t$hb{$i}\n"}' marker.list.sort |sort -k 1,1 > marker.list.sort.uniq

cut -f 1 marker.list.sort.uniq |while read i;do echo mkdir -p find_plot/$i;done|sh

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){@a=split;print "$ha{$a[0]}\t$_"}' cell.sort.list marker.list.sort.uniq > marker.list.sort.uniq.group

perl get.py.pl marker.list.sort.uniq.group

cut -f 2 marker.list.sort.uniq.group |while read i;do echo -e "cd find_plot/$i\nqsub -cwd -l vf=50G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q find.$i.sh\ncd -";done|sh

