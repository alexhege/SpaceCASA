#update array file
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors

from pylab import *

sc8rmse = np.loadtxt("rmses8.txt")
sc16rmse = np.loadtxt("rmses16.txt")
sc24rmse = np.loadtxt("rmses24.txt")
sc32rmse = np.loadtxt("rmses32.txt")

length = len(sc16rmse)

#plot(array(range(length))+1, sc8rmse, 'purple', label='8 SC')
plot(array(range(length))+1, sc16rmse, label='16 SC')
plot(array(range(length))+1, sc24rmse, label='24 SC')
plot(array(range(length))+1, sc32rmse, label='32 SC')


legend()

title("Comparing RMSE across Number of Spacecraft")
ylabel("RMSE")
xlabel("Minutes of Integration")

savefig('rmseCompare2.png')
