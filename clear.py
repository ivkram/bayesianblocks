import glob
import os

path = 'r vs g light curves (bayesian)/'+ '*.png'
files = glob.glob(path)
for src in files:
    os.remove(src)
    
