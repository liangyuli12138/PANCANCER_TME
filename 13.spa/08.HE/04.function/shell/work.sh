for i in {Bladder.B01613B1.sh,Bladder.B01615B2.sh,Breast.D01872C5.sh,Breast.D01972D1.sh,Breast.D01972D6.sh,Cervical.D01872C6.sh};do echo qsub -cwd -l vf=100G,num_proc=1 -P P22Z10200N0431_super -binding linear:1 -q st_supermem.q $i;done|sh

