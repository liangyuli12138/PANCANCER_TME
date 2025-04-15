import tifffile
import cv2
import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt

mask_path = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/04.mask.tif.from.qingdao/mask.tif/B01615B3_mask.tif" # cellseg_mask.tif
ssdna_path = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/05.cell_correct/regist_tif/B01615B3_regist.tif" # regist.tif
#out_path = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/05.cell_correct/regist_tif/test/test.out.tif" # mask_filter.tif
bright = 0 # avg_bright threshold of cells
bright = float(bright)

# load mask and ssDNA image
mask = tifffile.imread(mask_path)
mask = mask.astype(np.uint8)
ssdna_img = tifffile.imread(ssdna_path)
out_img = np.zeros(ssdna_img.shape,dtype = np.uint16)

num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(mask, connectivity=4)

list = []
for i in range(1,num_labels):
    # calculate cell area
    area = stats[i][4]
    if area < 80 or area > 3000:
        continue
    else:
        x1 = stats[i][0]
        x2 = stats[i][0] + stats[i][2]
        y1 = stats[i][1]
        y2 = stats[i][1] + stats[i][3]
        curr_mask=labels[y1:y2,x1:x2].copy()
        curr_mask[curr_mask!=i] = 0
        curr_mask[curr_mask==i] = 65535
        # calculate cell bright
        and_img = cv2.bitwise_and(ssdna_img[y1:y2,x1:x2].astype(np.uint16), curr_mask.astype(np.uint16))
        num_bright = float(np.sum(and_img))
        avg_bright = float(num_bright)/float(area)
        list.append(avg_bright)
        if avg_bright < bright:
            continue
        else:
            cut_img=labels[y1:y2,x1:x2].copy()
            id = np.where(cut_img==i)
            y= id[0] + y1
            x = id[1] + x1
            out_img[y,x] = 65535
plt.hist(list, bins=1000, alpha=0.7)
plt.ylim(0, 2500)
plt.xlim(0, 30000)
plt.savefig("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/05.cell_correct/bright_filter_test/result/B01615B3.bright.stat.png")
#cv2.imwrite("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/05.cell_correct/regist_tif/test/test.uint16-3.tif", out_img)
