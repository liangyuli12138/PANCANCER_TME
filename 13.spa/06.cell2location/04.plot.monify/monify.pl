###输入需要修正的恶性上皮样品
open IN0,$ARGV[0];while(<IN0>){chomp;$hmg{$_}=1}
###输入需要修正的正常上皮样品
open IN1,$ARGV[1];while(<IN1>){chomp;$hnr{$_}=1}
###输入需要修正的浆细胞样品
open IN2,$ARGV[2];while(<IN2>){chomp;$hbk{$_}="\"Lymphoid_Plamsa_IGKC\"";$hbl{$_}="\"Lymphoid_Plamsa_IGLC\"";}
open IN7,$ARGV[7];while(<IN7>){chomp;$hvn{$_}="\"EC_Vein\""}
open IN8,$ARGV[8];while(<IN8>){chomp;$hmi{$_}="\"Epithelium_Malig_Migration\""}

###输入HE画出来的恶性上皮spot的XY坐标
open IN3,$ARGV[3];while(<IN3>){chomp;@a=split;$b=$a[1]."_".$a[2];$him{$b}=1}

###匹配输入文件样品id
$tu=$ARGV[4];$id=$1 if($tu=~/Merge\_\S+\.(\S+)\_meta/);

###输出全部信息和简化信息两个文件
open OUT0,">$ARGV[5]";
open OUT1,">$ARGV[6]";

###如果不属于上述三个样品，直接输出原细胞
if(!exists $hmg{$id} && !exists $hnr{$id} && !exists $hbl{$id} && !exists $hvn{$id} && !exists $hmi{$id}){
open IN4,$ARGV[4];
$t=<IN4>;chomp $t;
print OUT0 "$t,\"celltype.mo\",\"location\"\n";
print OUT1 "\"\",\"id\",\"celltype\",\"max_score\",\"index\"\n";
while(<IN4>){
chomp;
s/\"POLYGON.+\)\)\"\,//;
@a=split(/,/);
$ip=$a[3]."_".$a[4];
undef %hz;for($i=10;$i<@a-1;$i++){$hz{$i}=$a[$i]};for $j(sort {$hz{$b} <=> $hz{$a}} keys %hz){$max=$hz{$j};last}
if($max<0.05){$a[-1]="Unknown"}
print OUT0 "$_,$a[-1],$ip\n";print OUT1 "$a[0],\"$ip\",$a[-1],$max,0\n"}}
else{

###否则，读原细胞文件，对每一行，根据三个情况进行分析判断
open IN4,$ARGV[4];
$t=<IN4>;chomp $t;
print OUT0 "$t,\"celltype.mo\",\"location\"\n";
print OUT1 "\"\",\"id\",\"celltype\",\"max_score\",\"index\"\n";
$t=~s/\"contour\"\,//;@ti=split(/,/,$t);

while(<IN4>){
chomp;
s/\"POLYGON.+\)\)\"\,//;
@a=split(/,/);
$ip=$a[3]."_".$a[4];
$o=$_;
$z=0;
undef %hz;for($i=10;$i<@a-1;$i++){$hz{$i}=$a[$i]};for $j(sort {$hz{$b} <=> $hz{$a}} keys %hz){$max=$hz{$j};last}
if($max<0.05){$last="Unknown";
print OUT0 "$_,$last,$ip\n";print OUT1 "$a[0],\"$ip\",$last,$max,1\n";$z=1;next}

$xxx=0;
if(exists $hbl{$id} || exists $hvn{$id} || exists $hmi{$id}){
if($a[-1] eq $hbl{$id}){for($i=10;$i<@a-1;$i++){if($ti[$i] eq $hbl{$id}){$a[$i]=$a[$i]*0.3;$xxx=1}}}
if($a[-1] eq $hbk{$id}){for($i=10;$i<@a-1;$i++){if($ti[$i] eq $hbk{$id}){$a[$i]=$a[$i]*0.3;$xxx=1}}}
if($a[-1] eq $hvn{$id}){for($i=10;$i<@a-1;$i++){if($ti[$i] eq $hvn{$id}){$a[$i]=$a[$i]*0.5;$xxx=1}}}
if($a[-1] eq $hmi{$id}){for($i=10;$i<@a-1;$i++){if($ti[$i] eq $hmi{$id}){$a[$i]=$a[$i]*0.5;$xxx=1}}}

if($xxx==1){
undef %hz;
for($i=10;$i<@a-1;$i++){$hz{$i}=$a[$i]};
for $j(sort {$hz{$b} <=> $hz{$a}} keys %hz){
if($a[$j]<0.05){$type="Unknown"}else{$type=$ti[$j]};
print OUT0 "$_,$type,$ip\n";print OUT1 "$a[0],\"$ip\",$type,$a[$j],2\n";$z=1;last}
next}}

if(exists $hmg{$id}){
if(!exists $him{$ip} && $a[-1]=~/Malig/ && !/Migration/){print OUT0 "$_,\"Epithelium_Normal\",$ip\n";print OUT1 "$a[0],\"$ip\",\"Epithelium_Normal\",$max,3\n";$z=1;next}}

if(exists $hnr{$id}){
if(exists $him{$ip} && $a[-1]=~/Normal/){
undef %hz;
for($i=10;$i<@a-1;$i++){$hz{$i}=$a[$i]};
$y=0;
for $j(sort {$hz{$b} <=> $hz{$a}} keys %hz){$y++;if($y==2){print OUT0 "$_,$ti[$j],$ip\n";print OUT1 "$a[0],\"$ip\",$ti[$j],$max,4\n";$z=1;last}}
next;}}


if($z==0){print OUT0 "$_,$a[-1],$ip\n";print OUT1 "$a[0],\"$ip\",$a[-1],$max,5\n";$z=1}
}}

print "finish\n"
