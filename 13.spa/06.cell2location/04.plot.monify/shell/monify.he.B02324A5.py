import numpy as np
import matplotlib.pyplot as plt
import cv2
import geopandas as gpd
from shapely.geometry import Polygon
from PIL import Image

he_path = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/00.he/B02324A5_malignant.jpg"

Image.MAX_IMAGE_PIXELS = None
mat = Image.open(he_path)
img = np.array(mat)

img = img[:,:,0]-img[:,:,1]

img[img<250] = 0
img[img>250] = 255

np.savetxt("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/04.plot.monify/cell.monify/B02324A5.img",img,fmt='%d')

