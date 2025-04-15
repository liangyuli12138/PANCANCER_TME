#!/usr/bin/perl

use strict;
use warnings;

# 读取list.file获取文件路径和关键id
open my $listFile, '<', './sn.list' or die "无法打开list.file: $!";
while (my $line = <$listFile>) {
    chomp $line;
    my ($id, $filePath) = split /\t/, $line;
#}
#close $listFile;

my $aa='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/'.$id.'_cellbin.final.celltype.obs.csv';
my $bb='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/02.subtypes.all/out/'.$filePath.'.'.$id.'.UCell.score.csv';
my $cc='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/04.hallmarker.all/out/'.$filePath.'.'.$id.'.UCell.score.csv';
my $dd='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/02.data/at/'.$id.'.at';

# 打开文件1
open my $file1, '<', $aa or die "无法打开文件1: $!";
# 打开文件3
open my $file3, '<', $bb or die "无法打开文件3: $!";
# 打开文件4
open my $file4, '<', $cc or die "无法打开文件4: $!";
# 创建输出文件
open my $output, '>', $dd or die "无法创建输出文件: $!";

open my $file2, '<', '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/02.data/sn.list' or die "无>法打开文件2: $!";

# 读取文件2的每一行的celltype并存储到数组中
my @celltypes;
while (my $line = <$file2>) {
    chomp $line;
    push @celltypes, $line;
}
close $file2;

# 写入输出文件的表头
print $output "index,region,Anti.tumor,Hallmarker_EMT\n";

# 遍历文件1的每一行
while (my $line = <$file1>) {
    chomp $line;
    my @fields = split /,/, $line;
    my $celltype = $fields[7];  # 文件1的celltype在第8列

    # 初始化输出行的值
    my $index = $fields[0];
    my $region = "other";
    my $antiTumor = "";
    my $hallmarkerEMT = "";

    # 判断文件1的celltype是否在文件2中存在
    if (grep { $_ eq $celltype } @celltypes) {
        $region = "Immune";
    }

    # 读取文件3，匹配Anti.tumor.microenvironment_UCell列的数据
    while (my $line3 = <$file3>) {
        chomp $line3;
        my @fields3 = split /,/, $line3;
        if ($fields3[0] eq $index) {
            $antiTumor = $fields3[4];
            last;
        }
    }
    # 将文件指针重置到文件开头，以便下次读取
    seek $file3, 0, 0;

    # 读取文件4，匹配HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_UCell列的数据
    while (my $line4 = <$file4>) {
        chomp $line4;
        my @fields4 = split /,/, $line4;
        if ($fields4[0] eq $index) {
            $hallmarkerEMT = $fields4[13];
            last;
        }
    }
    # 将文件指针重置到文件开头，以便下次读取
    seek $file4, 0, 0;

    # 输出结果到输出文件
    print $output "$index,$region,$antiTumor,$hallmarkerEMT\n";
}

# 关闭文件句柄
close $file1;
close $file3;
close $file4;
close $output;
}
close $listFile;

