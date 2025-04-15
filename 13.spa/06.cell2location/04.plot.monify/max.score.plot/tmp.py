atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_aaaa_meta_ori.monify.csv.list")
blist =  atlist["max_score"]
plt.hist(blist, bins=100, alpha=0.7, range=(0,2))
plt.title('aaaa')
plt.xlabel('Cell2location.max.score')
plt.ylabel('Cell_number')
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/ori_csv/Merge_aaaa_meta_ori.monify.max.png")
plt.close()

