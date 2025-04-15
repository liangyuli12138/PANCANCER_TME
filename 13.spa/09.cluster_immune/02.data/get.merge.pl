open  $listFile, '<', './sn.list' or die "无法打开list.file: $!";
while ( $line = <$listFile>) {
    chomp $line;
     ($id, $filePath) = split /\t/, $line;

$aa='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/'.$id.'_cellbin.final.celltype.obs.csv';
$bb='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/02.subtypes.all/out/'.$filePath.'.'.$id.'.UCell.score.csv';
$cc='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/04.hallmarker.all/out/'.$filePath.'.'.$id.'.UCell.score.csv';
$dd='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/02.data/at/'.$id.'.at';


open $file1, '<', $aa or die "无法打开文件1: $!";
open $file3, '<', $bb or die "无法打开文件3: $!";
open $file4, '<', $cc or die "无法打开文件4: $!";
open $output, '>', $dd or die "无法创建输出文件: $!";
open $file2, '<', '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/02.data/immune.list' or die "无>法打开文件2: $!";

undef %ha;undef %hb;undef %hc;

<$file3>;while(<$file3>){chomp;s/\"//g;@a=split(/,/);$ha{$a[0]}=$a[4]};
<$file4>;while(<$file4>){chomp;s/\"//g;@a=split(/,/);$hb{$a[0]}=$a[13]};
while(<$file2>){chomp;$hc{$_}="Immune"}


print $output "cell,region,antiTumor,hallmarkerEMT\n";
<$file1>;while(<$file1>){chomp;@a=split(/,/);$o=$hc{$a[-3]}?$hc{$a[-3]}:"other";if(exists $ha{$a[0]} && exists $hb{$a[0]}){print $output "$a[0],$o,$ha{$a[0]},$hb{$a[0]}\n"}}

}
