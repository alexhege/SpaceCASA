#!/usr/bin/env python
################################################################################
# PROGRAM: CME.py
################################################################################

'''
Alex Hegedus 8/10/17
alexhege@umich.edu
This code is used to generate many CASA truth images for RELIC to test
currently it makes a disk of emission, and a gaussian feature somewhere random
in the disk

run with casa, no arguments taken.  Edit key params in top of code

'''
#import the required python modules
import time
import os,sys
import shutil
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random
import math

### import casa modules
from rmtables import rmtables
from importvla import importvla
from importuvfits import importuvfits
from listobs import listobs
from applycal import applycal
from setjy import setjy
from gaincal import gaincal
from fluxscale import fluxscale
from accum import accum
#from plotxy import plotxy
from plotcal import plotcal
from clean import clean
from split import split
from clearcal import clearcal
from imstat import imstat
from imregrid import imregrid
from exportfits import exportfits
import casac
from taskinit import *
import numpy as np
import random
import math
from casa import imhead
from casa import *
from casac import *
from concat import concat
import signal
import time

import matplotlib.image as plimg
import scipy.ndimage.interpolation as spndint

from pylab import *
import random

# from joblib import Parallel, delayed
# import multiprocessing

# num_cores = multiprocessing.cpu_count()

## import image
#resolution
imSize = 512

#brightness of background gaussian, on top of which a smaller point src gaussian feature will be with (relStr-1)*baseStr brightness
baseStr = 10000. #jy


#different strength of gaussian compared to background disk
brightnesses = [1.1, 1.3, 1.5, 2.0, 4.0]

#how many trials to make per relative brightness
numTrials = 100

#brightestPix = 500000. #brightest pixel in jansky

ra = '03h00m00.0s'
dec = '-05d00m00.0s'

freq = 10e6 #observing frequency in Hz
imWidth = 400.*10 #arcseconds width of image
imageName = 'DRAGN-%3.3fMHz.truth'%(freq/1e6)

#width of image in radians
width = imWidth/206265. # arcsec to radians #res/60. * np.pi/180.
dres = width/imSize

phaseCenter = me.direction('J2000', ra, dec)



#Save all in folder
dir1=os.path.expandvars('$PWD')
# print dir1

datafolder = os.path.join(dir1, 'MCImages')
# print datafolder

if not os.path.exists(datafolder):
    os.makedirs(datafolder)

os.chdir(datafolder)



#Make gaussian features

#input just for
def createTruth(i):

    currfolder = os.path.join(relStrFolder, 'run_'+str(i))
    # print currfolder

    if not os.path.exists(currfolder):
        os.makedirs(currfolder)

    os.chdir(currfolder)

    ia.fromshape(imageName,shape=[imSize,imSize,1,1],overwrite=True)
    #adding components to  the empty image file
    gauss1 = me.shift(phaseCenter, offset=qa.toangle(str(width/4.)+'rad'), pa=qa.toangle(str(90)+'deg'))
    #gauss2 = me.shift(phaseCenter, offset=qa.toangle(str(width/4.)+'rad'), pa=qa.toangle(str(270)+'deg'))
    pt1 = me.shift(phaseCenter, offset=qa.toangle(str(random.random()*width/3.)+'rad'), pa=qa.toangle(str(random.random()*360)+'deg'))
    #pt2 = me.shift(gauss2, offset=qa.toangle(str(random.random()*width/8.)+'rad'), pa=qa.toangle(str(random.random()*360)+'deg'))

    SizeGauss = qa.toangle(str(2*width/3.)+'rad')
    SizePoint = qa.toangle(str(width/12.)+'rad')
    PA = qa.toangle(str(0)+'deg')


    # #factors to make flux of the max pixel baseStr, instead of the integrated gaussian flux being basestr
    # cl.addcomponent(dir=gauss1,
    #                 flux=baseStr*pi*(imWidth/8)**2/41, fluxunit='Jy', freq=freq/1e6,
    #                 shape="Gaussian", majoraxis=SizeGauss, minoraxis=SizeGauss, positionangle=PA)
    #
    # #factors to make flux of the entire disk baseStr, instead of the integrated flux being basestr
    # cl.addcomponent(dir=gauss2,
    #                 flux=baseStr*pi*(imWidth/8)**2/15, fluxunit='Jy', freq=freq/1e6,
    #                 shape="disk", majoraxis=SizeGauss, minoraxis=SizeGauss, positionangle=PA)

    cl.addcomponent(dir=phaseCenter,
                    flux=baseStr*pi*(imWidth/3)**2/15, fluxunit='Jy', freq=freq/1e6,
                    shape="disk", majoraxis=SizeGauss, minoraxis=SizeGauss, positionangle=PA)

    #factors to make flux of the entire disk baseStr, instead of the integrated flux being basestr
    # cl.addcomponent(dir=gauss2,
    #                 flux=baseStr*pi*(imWidth/8)**2/15, fluxunit='Jy', freq=freq/1e6,
    #                 shape="disk", majoraxis=SizeGauss, minoraxis=SizeGauss, positionangle=PA)

    cl.addcomponent(dir=pt1,
                    flux=baseStr*(relStr-1)*pi*(imWidth/24)**2/41, fluxunit='Jy', freq=freq/1e6,
                    shape="Gaussian", majoraxis=SizePoint, minoraxis=SizePoint, positionangle=PA)

    # cl.addcomponent(dir=pt2,
    #                 flux=baseStr*(relStr-1)*pi*(imWidth/48)**2/41, fluxunit='Jy', freq=freq/1e6,
    #                 shape="Gaussian", majoraxis=SizePoint, minoraxis=SizePoint, positionangle=PA)


    cs=ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])
    cell_rad=dres #qa.convert(qa.quantity("1arcmin"),"rad")['value']
    cs.setincrement([-cell_rad,cell_rad],'direction')
    cs.setreferencevalue([phaseCenter['m0']['value'], phaseCenter['m1']['value']],
                          type="direction")
    cs.setrestfrequency(freq)

    # add important header keywords
    imhead(imagename=imageName,mode="put",hdkey="object",hdvalue="DRAGN")
    imhead(imagename=imageName,mode="put",hdkey="imtype",hdvalue='Intensity')
    imhead(imagename=imageName,
          mode="put",hdkey="observer",hdvalue="simulation")
    imhead(imagename=imageName,
          mode="put",hdkey="date-obs",hdvalue="2023/03/15/00:00:00")
    imhead(imagename=imageName,mode="put",hdkey="reffreqtype",hdvalue='TOPO')
    imhead(imagename=imageName,
           mode="put",hdkey="restfreq",hdvalue=str(freq))
    imhead(imagename=imageName,mode='list')
    cs.setreferencevalue(str(freq)+'Hz', 'spectral')
    Telescope='VLA' #or else it breaks and whines
    cs.settelescope(Telescope)
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")

    ia.modify(cl.torecord(),subtract=False)


    ia.close()
    cl.done()
    qa.done()
    me.done()
    cs.done()

    os.chdir(relStrFolder)

    return True



for relStr in brightnesses:

    relStrFolder = os.path.join(datafolder, 'relStr_'+str(relStr))
    print relStrFolder

    if not os.path.exists(relStrFolder):
        os.makedirs(relStrFolder)

    os.chdir(relStrFolder)

    inputs = range(numTrials)

    # results = Parallel(n_jobs=num_cores)(delayed(createTruth)(i) for i in inputs)

    for i in inputs:
        temp = createTruth(i)


    os.chdir(datafolder)

os.chdir(dir1)
