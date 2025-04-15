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
if($n==1){$ha{$a[$i]}="Cluster_ori_".$id};
if($n==2){$hb{$a[$i]}="Cluster_ori_".$id};
}
}

open OUT,">D01972D1.ori.at";
print OUT "\"\",\"id\",\"celltype\"\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;$n++;
chomp;@b=split(/,/);
if(exists $ha{$b[0]}){print OUT "\"$n\",\"$b[0]\",\"$b[-3]\"\n"}else{print OUT "\"$n\",\"$b[0]\",\"other\"\n"};
}

open OUT,">D01972D1.ex.at";
print OUT "\"\",\"id\",\"celltype\"\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;$n++;
chomp;@b=split(/,/);
if(exists $hb{$b[0]}){print OUT "\"$n\",\"$b[0]\",\"$b[-3]\"\n"}else{print OUT "\"$n\",\"$b[0]\",\"other\"\n"};
}
