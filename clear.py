import glob
import os

path = '*.png'
files = glob.glob(path)
for src in files:
    os.remove(src)

path = '*_bb'
files = glob.glob(path)
for src in files:
    os.remove(src)

    
