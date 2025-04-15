import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt


atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Colorectal.D01872D2_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Colorectal.D01872D2')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Colorectal.D01872D2_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Colorectal.SS200000929BL_D2_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Colorectal.SS200000929BL_D2')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Colorectal.SS200000929BL_D2_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Colorectal.B01615B5_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Colorectal.B01615B5')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Colorectal.B01615B5_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Colorectal.D01872D1_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Colorectal.D01872D1')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Colorectal.D01872D1_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Breast.D01972D6_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Breast.D01972D6')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Breast.D01972D6_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Breast.D01972D1_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Breast.D01972D1')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Breast.D01972D1_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Breast.D01872C5_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Breast.D01872C5')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Breast.D01872C5_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Thyroid.D01872C4_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Thyroid.D01872C4')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Thyroid.D01872C4_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Thyroid.D01972D2_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Thyroid.D01972D2')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Thyroid.D01972D2_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Thyroid.D01972C1_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Thyroid.D01972C1')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Thyroid.D01972C1_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Gastric.D01972B6_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Gastric.D01972B6')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Gastric.D01972B6_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Gastric.B01615B1_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Gastric.B01615B1')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Gastric.B01615B1_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Gastric.D01872D3_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Gastric.D01872D3')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Gastric.D01872D3_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Live.B02324E3_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Live.B02324E3')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Live.B02324E3_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Live.B02324B5_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Live.B02324B5')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Live.B02324B5_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Live.B01613C6_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Live.B01613C6')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Live.B01613C6_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Kidney.B02324E4_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Kidney.B02324E4')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Kidney.B02324E4_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Kidney.B02324A5_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Kidney.B02324A5')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Kidney.B02324A5_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Kidney.B02324A1_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Kidney.B02324A1')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Kidney.B02324A1_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Bladder.B01613A5_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Bladder.B01613A5')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Bladder.B01613A5_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Bladder.B01615B2_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Bladder.B01615B2')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Bladder.B01615B2_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Bladder.B01613B1_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Bladder.B01613B1')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Bladder.B01613B1_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Ovarian.B02324E5_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Ovarian.B02324E5')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Ovarian.B02324E5_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Ovarian.B01613D1_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Ovarian.B01613D1')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Ovarian.B01613D1_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Cervical.D01872D4_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Cervical.D01872D4')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Cervical.D01872D4_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Cervical.D01972C4_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Cervical.D01972C4')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Cervical.D01972C4_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Cervical.D01872C6_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Cervical.D01872C6')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Cervical.D01872C6_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Endometrial.B01613B2_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Endometrial.B01613B2')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Endometrial.B01613B2_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Endometrial.B01317E6_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Endometrial.B01317E6')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Endometrial.B01317E6_meta_ori.monify.max.png")
plt.close()

atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Endometrial.B01615B3_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('Endometrial.B01615B3')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_Endometrial.B01615B3_meta_ori.monify.max.png")
plt.close()

