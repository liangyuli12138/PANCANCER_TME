open IN0,$ARGV[0];
while(<IN0>){
chomp;if(/data\/(\S+)\_cellbin/){$id=$1};
open IN,$_;
<IN>;
while(<IN>){
s/Lymphoid_NK_CD56\.\,nan/Lymphoid_NK_CD56\.\,NK_cell/g;s/Lymphoid_NK_CD16\.\,nan/Lymphoid_NK_CD16\.\,NK_cell/g;
chomp;@a=split(/,/);
$ha{$id}{$a[0]}="$a[-3],$a[-2]";
}}

open IN1,$ARGV[1];$t=<IN1>;chomp $t;print "$t,celltype,celltype_merge\n";
while(<IN1>){chomp;@a=split(/,/);print "$_,$ha{$a[1]}{$a[2]}\n"}
