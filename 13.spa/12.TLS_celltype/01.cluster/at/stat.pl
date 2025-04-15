while(<>){chomp;@a=split(/,/);if(@a<6){if(!exists $ha{$a[-1]}){$n++;$ha{$a[-1]}=1}}};print "$n\n"
