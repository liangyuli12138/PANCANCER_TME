$r=$ARGV[0];$r=~s/\.\.\/at\///;
open OUT,">at/$r";
$t=<>;print OUT "\"\",\"id\",\"celltype\"\n";
$n++;
while(<>){$n++;if(/other/){@a=split(/,/);print OUT "\"$n\",\"$a[0]\",\"others\"\n"}else{@a=split(/,/);print OUT "\"$n\",\"$a[0]\",\"TLS\"\n"}}
