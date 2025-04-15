$r=$ARGV[0];$r=~s/\.\.\/at\///;
open OUT,">at/$r";
$t=<>;print OUT "\"\",\"id\",\"celltype\"\n";
$n++;
while(<>){chomp;$n++;if(/other/){@a=split(/,/);print OUT "\"$n\",\"$a[0]\",\"others\"\n"}else{@a=split(/,/);$a[1]=~s/\S+\_TLS/TLS/;print OUT "\"$n\",\"$a[0]\",\"$a[1]\"\n"}}
