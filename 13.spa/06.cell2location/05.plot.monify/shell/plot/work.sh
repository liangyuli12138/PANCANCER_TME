cut -f 1 sn.list |while read i;do echo qsub -cwd -l vf=30G,num_proc=1 -P P22Z10200N0433 -binding linear:1 plot.$i.sh;done|sh

