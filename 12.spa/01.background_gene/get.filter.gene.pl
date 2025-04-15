while(<>){
chomp;
if(!/Lymphoid_Plamsa/ || !/^IG/){next};
@a=split;s/\{//g;s/\}//g;
if($a[4]<0.1){next};
@b=split(/;/,$a[-1]);
@c=split(/\|/,$b[0]);$m=int($a[5]*1.5);
if($a[4]>($c[2]*1.5)){
print "$_\t$m\n"}
}
