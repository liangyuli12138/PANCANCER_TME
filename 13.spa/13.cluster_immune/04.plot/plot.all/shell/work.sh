for i in {plot.B01613A5.sh,plot.B02324A1.sh,plot.B02324E5.sh,plot.D01872D1.sh,plot.D01872D2.sh};do echo qsub -clear -cwd -l vf=10G,num_proc=1 -P P22Z10200N0431 -binding linear:1 -q st.q $i;done|sh

