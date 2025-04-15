<>;
while(<>){
chomp;
@a=split(/,/);
undef %ha;
for($i=11;$i<@a-3;$i++)
{$ha{$a[$i]}=1};
$n=0;
print "$a[-2]";
for $j(sort {$b<=>$a} keys %ha)
{$n++;if($n>2){last};print "\t$j"}
print "\n"
}
