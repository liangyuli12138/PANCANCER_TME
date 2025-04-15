import os
import datetime
import argparse
import glob
import subprocess
import pandas as pd

# cell seg correct demo
from cellbin.modules.cell_labelling import CellLabelling
from stio.matrix_loader import MatrixLoader
from cellbin.image import Image

###### Version and Date
PROG_VERSION = 'StereoCell Python Wrapper: v1.1.0'
PROG_DATE = '2023-08-25'

###### Usage
usage="""StereoCell Wrapper v1.1.0"""

parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("--version", action="version", version=PROG_VERSION)
parser.add_argument("-m",dest="cell_mask_path")
parser.add_argument("-g", "--gem", dest="gem_path")
parser.add_argument("-o", "--output",  dest="output_dir")
parser.add_argument("-d", "--distance",  dest="distance")
args = parser.parse_args()

print(PROG_VERSION)
print(PROG_DATE)
print(args)
cell_mask_path=args.cell_mask_path
gem_path=args.gem_path
output_path=args.output_dir
distance=args.distance

# read cell mask
i = Image()
i.read(image=cell_mask_path)
mask = i.image
mask.shape

# read gene matrix
new_file = gem_path

# 修正
cl = CellLabelling(
    mask,
    new_file
)
cl.set_process(10)
correct_mask, exp_matrix = cl.run_fast(distance=20)

# mask写出
correct_mask_path = os.path.join(output_path, "cell_mask_corect.compression.tif")
Image.write_s(correct_mask, correct_mask_path, compression=True)
correct_mask_path = os.path.join(output_path, "cell_mask_corect.tif")
Image.write_s(correct_mask, correct_mask_path, compression=False)

# gem写出
correct_out_gem = os.path.join(output_path, "cell_mask_corect.gem")
exp_matrix.to_csv(correct_out_gem,sep='\t',index_col=0)

