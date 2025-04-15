open IN0,$ARGV[0];
while(<IN0>){
chomp;if(/Boundary/){next};
$num=$1 if(/Size\:\s+(\d+)/);
$id=$1 if(/Cluster\s+(\d+)/);
$id="Cluster"."$id";
s/^.+\)\:\s+//;
s/\s+//g;
@a=split(/\,/);
for($i=0;$i<@a;$i++){
$ha{$a[$i]}=$id;
}
}

print "\"\",\"id\",\"celltype\"\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;$n++;
chomp;@b=split(/,/);
if(exists $ha{$b[0]}){print "\"$n\",\"$b[0]\",\"$ha{$b[0]}\"\n"}else{print "\"$n\",\"$b[0]\",\"other\"\n"};
}
