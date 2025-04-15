import tifffile
import cv2
import numpy as np
import pandas as pd
import sys
import os

mask_path = sys.argv[1] # cellseg_mask.tif
ssdna_path = sys.argv[2] # regist.tif
out_path = sys.argv[3] # mask_filter.tif
bright = sys.argv[4] # avg_bright threshold of cells
bright = float(bright)

# load mask and ssDNA image
mask = tifffile.imread(mask_path)
mask = mask.astype(np.uint8)
ssdna_img = tifffile.imread(ssdna_path)
out_img = np.zeros(ssdna_img.shape,dtype = np.uint8)

num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(mask, connectivity=4)

list = []
for i in range(1,num_labels):
    # calculate cell area
    area = stats[i][4]
    if area < 0 or area > 2500:
        continue
    else:
        x1 = stats[i][0]
        x2 = stats[i][0] + stats[i][2]
        y1 = stats[i][1]
        y2 = stats[i][1] + stats[i][3]
        curr_mask=labels[y1:y2,x1:x2].copy()
        curr_mask[curr_mask!=i] = 0
        curr_mask[curr_mask==i] = 255
        # calculate cell bright
        and_img = cv2.bitwise_and(ssdna_img[y1:y2,x1:x2].astype(np.uint8), curr_mask.astype(np.uint8))
        num_bright = float(np.sum(and_img))
        avg_bright = float(num_bright)/float(area)
        if avg_bright < bright:
            continue
        else:
            list.append(avg_bright)
            cut_img=labels[y1:y2,x1:x2].copy()
            id = np.where(cut_img==i)
            y= id[0] + y1
            x = id[1] + x1
            out_img[y,x] = 255

cv2.imwrite(out_path, out_img)
