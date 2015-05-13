import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def reading_data(filename):
    data = np.loadtxt("../data/CLARA/"+filename)
    initpos=data[:,0]
    index_clean = np.where(~np.isnan(initpos))
    data = data[index_clean[0],:]
    kz = data[:,5]
    x = data[:,6]

    return kz, x

def viewing_angle(filename, angle):
    kz, x = reading_data(filename)
    k = angle/90.0
    index = np.where((abs(kz)<(k+0.1)) & (abs(kz)>(k-0.1)))
    nbins =40
    x = x[index]
    hist, bins = np.histogram(x, bins=nbins)

    return hist, bins
 

filename = sys.argv[1]
angle = float(sys.argv[2])

hist, bins = viewing_angle(filename, angle)

for i in range(len(hist)):
	print hist[i], bins[i]
