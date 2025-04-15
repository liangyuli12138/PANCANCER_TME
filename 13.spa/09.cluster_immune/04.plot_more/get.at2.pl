open IN0,$ARGV[0];
while(<IN0>){
if(/ori/){$n=1}else{$n=2};
chomp;if(/Boundary/){next};
$num=$1 if(/Size\:\s+(\d+)/);
$id=$1 if(/Cluster\S+\s+(\d+)/);
#$id="Cluster"."$id";
s/^.+\)\:\s+//;
s/\s+//g;
@a=split(/\,/);
for($i=0;$i<@a;$i++){
if($n==1){$ha{$a[$i]}="Cluster".$id};
if($n==2){$hb{$a[$i]}="Cluster".$id};
}
}

$id=$ARGV[2];

open OUT,">at100/$id.ori.at";
print OUT "\"\",\"id\",\"celltype\",cluster\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;$n++;
chomp;@b=split(/,/);
if(exists $ha{$b[0]}){print OUT "\"$n\",\"$b[0]\",\"$b[-3]\",$ha{$b[0]}\n"}else{print OUT "\"$n\",\"$b[0]\",\"other\",\"other\"\n"};
}

open OUT,">at100/$id.ex.at";
print OUT "\"\",\"id\",\"celltype\",cluster\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;$n++;
chomp;@b=split(/,/);
if(exists $hb{$b[0]}){print OUT "\"$n\",\"$b[0]\",\"$b[-3]\",$hb{$b[0]}\n"}else{print OUT "\"$n\",\"$b[0]\",\"other\",\"other\"\n"};
}
