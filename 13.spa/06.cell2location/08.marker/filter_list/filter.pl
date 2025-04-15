open OUT0,">$ARGV[1]";
open OUT1,">$ARGV[2]";
open IN,$ARGV[0];
$t=<IN>;
print OUT0 "cell\n";
print OUT1 "cell,celltype_sc\n";
while(<IN>){
chomp;
@a=split(/,/);
#if($a[-2] eq "B_cell" || $a[-2] eq "DC" || $a[-2] eq "EC" || $a[-2] eq "Epithelium" || $a[-2] eq "Fibroblast" || $a[-2] eq "Marcophage" || $a[-2] eq "Mural_cell" || $a[-2] eq "NK_cell")
$x=0;
if($a[-2] eq "B_cell"){$o="Lymphoid_B";$x=1};
if($a[-2] eq "DC" || $a[-2] eq "Marcophage" || $a[-2] eq "Myeloid_other"){$o="Myeloid";$x=1};
if($a[-2] eq "EC"){$o="EC";$x=1};
if($a[-2] eq "Epithelium"){$o="Epithelium";$x=1};
if($a[-2] eq "Fibroblast"){$o="Fibroblast";$x=1};
if($a[-2] eq "Mural_cell"){$o="Mural_cell";$x=1};
if($a[-2] eq "NK_cell" || $a[-2] eq "T_cell"){$o="T/NK";$x=1};
if($a[-2] eq "Plamsa_cell"){$o="Plamsa_cell";$x=1};
if($a[-3] eq "Myeloid_Mast"){$o="Mast";$x=1};

if($x==1){print OUT0 "$a[0]\n";print OUT1 "$a[0],$o\n";}

}
