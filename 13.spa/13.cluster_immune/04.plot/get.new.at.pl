open IN0,$ARGV[0];
while(<IN0>){
chomp;
s/\"//g;
@a=split(/,/);
if($a[2] ne "malignant" && $a[2] ne "other"){$ha{$a[1]}=$a[2]}
}

open IN1,$ARGV[1];
$t=<IN1>;chomp $t;
print "id","$t,","celltype\n";
while(<IN1>){
chomp;
@a=split(/,/);
if(exists $ha{$a[0]}){print "$_,$ha{$a[0]}\n"}
elsif($a[1] eq "malignant"){print "$_,malignant\n"}
else{print "$_,other\n"}
}

