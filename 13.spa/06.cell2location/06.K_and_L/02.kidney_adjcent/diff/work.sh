perl -e 'while(<>){chomp;@a=split(/\//);@b=split(/\_/,$a[-1]);print "mkdir -p diff/$b[0]\n"}'  all.obs.list |sh

perl -e 'while(<>){chomp;open IN,$_;$n=0;$t=$_;@z=split(/\//,$t);@b=split(/\_/,$z[-1]);<IN>;undef %ha;while(<IN>){chomp;@a=split(/\,/);if(!exists $ha{$a[-1]}){$n++;$ha{$a[-1]}=1}};$x=$t;$t=~s/obs\.txt/h5ad/;for($i=0;$i<$n;$i++){print "/hwfssz4/BC_PUB/Software/07.User-defined/03.Animal_Plant/wubin/mamba/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/06.K_and_L/02.kidney_adjcent/diff/findmarker.py $t $i /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/06.K_and_L/02.kidney_adjcent/diff/diff/$b[0]/diff.$i.csv\n"}}'  all.obs.list > all.diff.sh


