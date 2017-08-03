#!/usr/bin/env python
################################################################################
# PROGRAM: CME.py
################################################################################

'''
This simulation code has been put together to read orbit files and simulate space based radio arrays

instrution on how to run the program
within casa run this command:

%run corrRelicGif.py


'''

#Parameters

#True: simulate manual correlation to form visibilities (slow)
# False: use CASA sm.predict to form visibilities (much faster)
corrMan = False

freq = 10e6 #observing frequency in Hz

antennaDiameter = 6.0 #meters

#iterations of cleaning to do while imaging
cleanIters = 50

#create and move to output folder, absolute or relative path (to where you called casa if in casa) detected
outDir = 'RelicOut'

#will be interpolated to imsize x imsize, absolute or relative path detected
imName = 'cyga_21cm.png'

#resolution to extrpolate input image to, what CASA will use.  power of 2 better
imSize = 1024

#what casa truth image copy of input truth png will be named
imageName = 'DRAGN-%3.3fMHz.truth'%(freq/1e6)

#True: form image & compute RMSE at each time step
#False: only form image & compute RMSE at end
gifMode = True


imWidth = 800.  #arcseconds width of image
brightestPix = 500000. #brightest pixel in jansky

#location of image brightness in the sky, format correctly for CASA!
ra = '03h00m00.0s'
dec = '-05d00m00.0s'


#formatted orbit file, absolute or relative path detected
#'/Users/hegedus/Downloads/Relic_Pipeline/inert_traj_061216_reconfig.txt'
orbitFile = 'inert_traj_061216_reconfig.txt'

#indexes into time column of spacecraft positions in the orbitFile
ranges =  ((0, 10), (20, 30)) # or ((0, 2), (23, 28))

timeStep = 60. #timestep in seconds between each sample in orbit file

#choose number of spacecraft to test, takes off random antennae if less than number in orbit file
wantedSC = 32

#positional error with sigma dTau nanoseconds of light travel time
posErr = True
dTau = 5e-9

thermalNoise = True

#make number in Jy or 'galactic' for simple power law approx of galactic noise
tNoise = 'galactic'

##################

#import the required python modules
import time
import os,sys
import shutil
from collections import defaultdict
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



print 'done with imports'

############


#end Params, do preliminary calcs
casaFactor = 10
imWidth = imWidth*casaFactor #factor to allow casa to work, ignore
#makes image ten times bigger with 10 times shorter baselines as a workaround to CASA
#not liking 600 km baselines

img = plimg.imread(imName).astype(np.float32)
pfile = open(orbitFile)


if tNoise == 'galactic':
    Tgal = 9.e7*(freq/10e6)**(-2.477)
    tNoise = Tgal


dir1=os.path.expandvars('$PWD')
#create and move to output folder
if outDir[0] == '/':
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    os.chdir(outDir)

#assume relative path
else:
    #Save all in folder
    datafolder = os.path.join(dir1, outDir)
    # print datafolder
    if not os.path.exists(datafolder):
        os.makedirs(datafolder)

    os.chdir(datafolder)


#width of image in radians
width = imWidth/206265. # arcsec to radians #res/60. * np.pi/180.
c = 3e8
kwavnum = 2*np.pi*freq
wavelen = c/freq
kb = 1.380648e-23
#calculate width of each pixel
dres = width/imSize

#position error in seconds
poserr = 3e8*dTau #3. #3000. meters
dv = poserr/wavelen
dphase = 2*np.pi*dv





##############read in image


print 'now creating CASA truth image'


phaseCenter = me.direction('J2000', ra, dec)



dims = np.shape(img)
print 'dims of png is ' + str(dims)
avimg = img
if len(dims) == 2:
    avimg = img
    d1 = float(max(dims))
else:
    d3 = min(3,dims[2])
    d1 = float(max(dims))
    avimg = np.average(img[:,:,:d3],axis=2)


avimg -= np.min(avimg)
avimg *= brightestPix/np.max(avimg)

zoomimg = spndint.zoom(avimg,float(imSize)/d1)
zdims = np.shape(zoomimg)

for i in range(zdims[0]):
    for j in range(zdims[1]):
        if zoomimg[i][j] < 0.0:
            zoomimg[i][j] = 0.0

zoomimg -= np.min(zoomimg)
zoomimg *= brightestPix/np.max(zoomimg)


if zdims[0] != zdims[1]:
    newarr = np.zeros((imSize, imSize))
    d2 = min(zdims)
    if zdims[0] == d2:
        for i in range(d2):
            newarr[imSize/2 - d2/2 + i][:] = zoomimg[i][:]
    elif zdims[1] == d2:
        for i in range(d2):
            newarr[:][imSize/2 - d2/2 + i] = zoomimg[:][i]
    zoomimg = newarr

z = zoomimg.copy().T
z= np.fliplr(z)  #these operations flip to CASA style of storing data
#which starts at lower left corner, and goes up by columns left to right


casaArr = z.reshape((imSize, imSize, 1, 1))

ia.fromarray(imageName, pixels=casaArr, overwrite=True)

#simple point src test
# ia.fromshape(imageName,shape=[imSize,imSize,1,1],overwrite=True)
# delta2 = qa.toangle(str(1)+'deg')
# offangle2 = qa.toangle(str(45)+'deg')
# NewDir2 = me.shift(phaseCenter, offset=delta2, pa=offangle2)
# cl.addcomponent(dir=NewDir2, flux=7e9, fluxunit='Jy', freq=freq, shape='point')
# cl.addcomponent(dir=NewDir2,
#                 flux=1e10, fluxunit='Jy', freq=freq,
#                 shape="Gaussian", majoraxis=qa.toangle(str(2)+'deg'), minoraxis=qa.toangle(str(1)+'deg'), positionangle=qa.toangle(str(45)+'deg'))

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

# ia.modify(cl.torecord(),subtract=False)


pix = ia.getchunk()
pix = pix.reshape((imSize, imSize))

zoomimg = pix


ia.close()
cl.done()
qa.done()
me.done()

print 'Now importing antenna positions'
#####import antenna positions
#77200 times, 1 min each


times = []
scpos = []

#get into arrays, preprocess strings a bit
for line in pfile.readlines():
    line = line.split(',')
    times.append(line[0])
    line[-1] = line[-1][:-2]
    scpos.append(line[1:])


numsc = np.shape(scpos)[1]/3


snips = []
timesnips = []

if len(np.shape(ranges)) == 1:
  snips += scpos[ranges[0]:ranges[1]]
  timesnips += times[ranges[0]:ranges[1]]
else:
  for r in ranges:
    snips += scpos[r[0]:r[1]]
    timesnips += times[r[0]:r[1]]

scpos = snips
times = timesnips

print 'length of scpos '+str(len(scpos))

#move into np array, convert from strings to float with map
numbl = numsc*(numsc - 1)/2
nscpos = np.zeros((len(scpos), numsc*3))
for i in range(len(scpos)):
  nscpos[i] = np.array(map(float, scpos[i]))


#reshape array so first indexed by spacecraft

length = len(scpos) #number of positions for each s/c
allscpos = np.zeros((numsc, length, 3))
for i in range(numsc):
  allscpos[i] = nscpos[:, i*3:(i+1)*3]

##take away random antenna
takeOff=0
if wantedSC > 0 and wantedSC < numsc:
    takeOff = int(numsc - wantedSC)

while takeOff != 0:
    ind = random.randint(0, numsc-1)
    allscpos = np.delete(allscpos, ind, axis=0)
    takeOff -= 1
    numsc -= 1

numbl = numsc*(numsc - 1)/2

allscpos *= 1e3 / casaFactor #3 #convert to meters, but /10 from imSize change

#one time things to compute
dishDiameter = np.ones(numsc)*antennaDiameter

mounts = []
antnames = []

for i in range(numsc):
    mounts.append('X-Y')
    antnames.append('S'+str(i))

###### MS stuff

RelicMS = 'RELIC-%3.3fMHz_time_'%(freq/1e6)+times[0]+'.ms'

RelicMergeMS = 'RELICMERGE-%3.3fMHz_time_'%(freq/1e6)+times[0]+'.ms'

phasescpos = np.zeros((numsc,length, 3))

#here loop over length for each time sample

visibilities = np.zeros((numbl, length), dtype=np.complex128)
casavisibilities = np.zeros((numbl, length), dtype=np.complex128)
casaabsq = np.zeros((numbl, length))
baselines = np.zeros((numbl, length, 3), dtype=np.float64)




print 'reprojecting baselines'
srcdir = phaseCenter
ra = srcdir['m0']['value'] #returns radians
dec = srcdir['m1']['value'] #returns radians
res = imSize

x = np.cos(dec) * np.cos(ra)
y = np.cos(dec) * np.sin(ra)
z = np.sin(dec)
s = np.array([x, y, z])


#get x and y axes perp to s
x = np.cross(s, [0, 0, 1])
y = np.cross(s, x)

#normalize since cross product isn't already 1 for some reason
# this is
px = x/np.linalg.norm(x)
py = y/np.linalg.norm(y)
pz = s/np.linalg.norm(s)

for l in range(length):

    for i in range(numsc):
        phasescpos[i][l][0] = np.dot(allscpos[i][l], px)
        phasescpos[i][l][1] = np.dot(allscpos[i][l], py)
        phasescpos[i][l][2] = np.dot(allscpos[i][l], pz)

    k=0
    for i in range(numsc):
      for j in range(i+1, numsc):
        baselines[k][l] =  phasescpos[j][l] - phasescpos[i][l]
        k+=1


#find longest projected baseline for cleaned image resolution later

norms = np.zeros((numbl, length))

for i in range(numbl):
    for j in range(length):
        norms[i][j] = np.linalg.norm(baselines[i][j][:2])

largestBL = np.amax(norms)


print 'starting loop of visibility calculations & MS creation'
#$

print corrMan

RelicCurrMSs = []
jpgNames = []
for l in range(length):

    RelicCurrMS = 'RELIC-%3.3fMHz_time_'%(freq/1e6)+times[l]+'.ms'
    RelicCurrMSs.append(RelicCurrMS)

    sm.open(RelicCurrMS)
    print 'opened'
    refloc=me.observatory('VLA')

    Ant_x = allscpos[:, 0, 0]
    Ant_y = allscpos[:, 0, 1]
    Ant_z = allscpos[:, 0, 2]

    #set the antenna configuration
    sm.setconfig(telescopename=Telescope,
                 x=Ant_x,y=Ant_y,z=Ant_z,
                 dishdiameter=dishDiameter,
                 mount=mounts,
                 antname=antnames,
                 coordsystem='local',referencelocation=refloc)



    sm.setspwindow(spwname='HF',
                   freq=str(freq/1e6)+'MHz',
                   deltafreq='12200Hz', freqresolution='12200Hz', nchannels=1)
    # Set the feeds to be X Y.  Haven't figured it out.
    print 'Setting the feed polarization ...'
    sm.setfeed('perfect X Y')

    #
    print 'Setting the source ...'

    srcdir = phaseCenter
    sm.setfield(sourcename='CME1', sourcedirection=srcdir); #CMEDir); ##Change

    #
    print 'Setting integration times ...'
    rday1 = 60017.83
    rday1 += l/60./24 *(timeStep/60.)#add 1 min per length
    refday=me.epoch(rf='UTC',v0=str(rday1)+'d')
    sm.settimes(integrationtime=str(timeStep)+'s',usehourangle=False,referencetime=refday)

    #
    # Don't bother with autocorrelations
    sm.setauto(0.0);
    print 'Observing ...'
    sm.observe(sourcename='CME1',
               spwname='HF',
               observer='SunRISE',
               starttime='0s', stoptime=str(1*timeStep)+'s')
    # Fourier transform model image into data set


    print ' now on time ' + str(l)




    if corrMan:
        print 'doing manual correlation on step '+str(l)
        print corrMan

        #which s/c is farthest down along pz?
        minZ = np.where(phasescpos[:,l,2] == min(phasescpos[:,l,2]))[0][0]

        #get delays in seconds, from diff in z from minZ
        delays = np.zeros(numsc)
        for i in range(numsc):
            delays[i] = phasescpos[i][l][2] - phasescpos[minZ][l][2]
            delays[i] = delays[i]/c

        #change later for more frequencies
        fourierCoeff = np.zeros(numsc, dtype=np.complex128)

        #print np.shape(zoomimg)
        d1, d2 = np.shape(zoomimg)

        #loop over other pixels, adding in apporpriate delays for each sc fourierCoeff
        #like calculation of phase center, but t in front of variables for temporary, don't keep
        for i in range(d1):
            for j in range(d2):
                if zoomimg[i][j] == 0.0:
                    continue

                s = cs.toworld([i,j,0,0],'s')['string']

                tra = qa.convert(qa.toangle(s[0]), 'rad')['value']
                tdec = qa.convert(qa.toangle(s[1]), 'rad')['value']

                # print tra, tdec, i, j

                #unit direction pointing in direction of source/phase center
                tx = np.cos(tdec) * np.cos(tra)
                ty = np.cos(tdec) * np.sin(tra)
                tz = np.sin(tdec)


                ts = np.array([tx, ty, tz])


                # #get x and y axes perp to s

                tx = np.cross(ts, [0,0,1])
                ty = np.cross(ts, tx)

                #normalize since cross product isn't already 1 for some reason
                # this is
                tx = tx/np.linalg.norm(tx)
                ty = ty/np.linalg.norm(ty)
                tz = ts/np.linalg.norm(ts)


                tphasescpos = np.zeros((numsc,3))

                for k in range(numsc):
                    tphasescpos[k][0] = np.dot(allscpos[k][l], tx)
                    tphasescpos[k][1] = np.dot(allscpos[k][l], ty)
                    tphasescpos[k][2] = np.dot(allscpos[k][l], tz)

                #which s/c is farthest down along pz?
                tminZ = np.where(tphasescpos[:,2] == min(tphasescpos[:,2]))[0][0]

                tdelays = np.zeros(numsc)
                for k in range(numsc):
                    tdelays[k] = tphasescpos[k][2] - tphasescpos[tminZ][2]
                    tdelays[k] = tdelays[k]/c


                #phase of furthest down z in this projection, delay from this one
                phase = np.exp(1j*2*np.pi*random.random())
                C = zoomimg[i][j]*phase
                for k in range(numsc):
                    fourierCoeff[k] = C*np.exp(-1j*kwavnum*(tdelays[k] - delays[k]))
                    # print 'i j k four', i, j, k, C*np.exp(-1j*kwavnum*(tdelays[k] - delays[k]))
                    # print 'delays t and p phase shift', tdelays[k], delays[k], np.exp(-1j*kwavnum*(tdelays[k] - delays[k]))

                k=0
                for m in range(numsc):
                  for n in range(m+1, numsc):
                    toAdd = fourierCoeff[n]*fourierCoeff[m].conj()
                    toAdd = toAdd/np.sqrt(abs(toAdd))
                    visibilities[k][l] += toAdd
                    #print 'k, toAdd, vis ', k, toAdd, visibilities[k]
                    k+=1



    tb.open(RelicCurrMS, nomodify=False)



    actual_uvw = tb.getcol("UVW")

    baselines2 = np.reshape(baselines[:,l,:], (numbl,3), order='F').T


    tb.putcol("UVW", baselines2)

    print 'Predicting the visibilities ...'

    sm.predict(imagename=imageName)




    if corrMan:

        data=tb.getcol("DATA")
        cdata=tb.getcol("CORRECTED_DATA")

        #"Now replacing CASA data with self made correlated data, comparing it"
        for i in range(numbl):
            # for j in range(l, l+1):
            #"check casa uvw and calculated"
            #print actual_uvw[:,i], baselines[:,i]

            q = data[0][0][i]/visibilities[i][l]

            # print data[0][0][i], visibilities[i][l], q, abs(q), q/abs(q)


            casaabsq[i][l] = abs(q)
            casavisibilities[i][l] = q/abs(q) #data[0][0][i]


            data[0][0][i] = visibilities[i][l]
            data[1][0][i] = visibilities[i][l]

            cdata[0][0][i] = visibilities[i][l]
            cdata[1][0][i] = visibilities[i][l]


        tb.putcol("DATA", data)
        tb.putcol("CORRECTED_DATA", cdata)


    if posErr:
        data=tb.getcol("DATA")
        cdata=tb.getcol("CORRECTED_DATA")
        #"Now replacing CASA data with self made correlated data, comparing it"
        for i in range(numbl):
            for j in range(0, 1):
                phaseErr = np.random.normal(0, dphase*np.sqrt(2./3)) #sqrt 2/3 : 2 from 1 for each sc pos, div by 3 to go from total pos error to only z direction
                data[0][0][j*numbl+i] *= exp(1.j*phaseErr)
                data[1][0][j*numbl+i] *= exp(1.j*phaseErr)

        tb.putcol("DATA", data)
        tb.putcol("CORRECTED_DATA", cdata)

        print("corrupting with phase Noise")


    #change flags, being stupid
    # f = tb.getcol("FLAG")
    # for i in range(numbl):
    #     for j in range(length):
    #         f[0][0][j*numbl+i] = False
    #         f[1][0][j*numbl+i] = False
    #
    #
    # tb.putcol("FLAG", f)



    tb.close()
            #sm.open(CMEms)

        ###### end newwww code

    ## add errors (noise, gains, polarization leakage, bandpass) to the vsiibility data
    if thermalNoise:
        sm.setnoise(mode='simplenoise',simplenoise=str(tNoise)+'Jy')
        print("corrupting with thermal Noise")
        sm.corrupt()

    #Clean here?

    sm.close()

    concat(vis=RelicCurrMS,
           concatvis=RelicMergeMS,
           freqtol='',
           dirtol='',
           respectname=False,
           timesort=False,
           copypointing=True,
           visweightscale=[])


    #now image

    #if not taking images every time step and not on last step, don't image
    if gifMode == False and l != (length - 1):
        continue

    print 'cleaning dirty image'

    #assume 600 km baseline is largest, lambda/D/4

    cellsizerad = wavelen/largestBL/4.
    csas = cellsizerad*180./np.pi*3600

    #csas = cellsizerad

    imsize = int(imWidth/csas) #int(1024/60.*3600./csas)
    print 'pre imsize is ' + str(imsize)
    if imsize < 1:
        imsize = 1024

    #use next power of 2
    imsize = int(pow(2, math.ceil(math.log(imsize, 2))))
    print 'imsize used is ' + str(imsize)

    print 'cell size for ' + RelicMergeMS + ' is ' + str(csas) +' arcsec'

    try:
        clean(vis=RelicMergeMS,imagename=RelicMergeMS+'.dirty'+str(l),
            outlierfile='',
            field='',spw='',
            selectdata=False,
            mode='mfs',nterms=1,
            gridmode='widefield', wprojplanes = -1, facets = 1,
            niter=cleanIters,gain=0.1,threshold='0.0mJy',
            interactive=False,
            mask=[],
            imsize=[imsize, imsize], #[1024, 1024],
            cell=[str(csas)+'arcsec', str(csas)+'arcsec'] ,#cell=['1arcsec', '1arcsec'],
            phasecenter='',
            stokes='I',
            weighting='briggs',robust=-0.5,
            uvtaper=False,
            usescratch=False,
            allowchunk=False
            )

        viewer(infile=RelicMergeMS+'.dirty'+str(l)+'.image',displaytype='raster',
              outfile='RELICMERGE-%3.3fMHz_time_'%(freq/1e6)+str(l)+'.jpg',outformat='jpg',
              gui=False)


        jpgNames.append('RELICMERGE-%3.3fMHz_time_'%(freq/1e6)+str(l)+'.jpg')

    except:
        print 'Clean wanted to crash on ' + RelicMergeMS
        pass


#plot vis diff
if corrMan:
    print 'making comparison plots of correlated vs casa simulated visibilities'

    relErrors = abs(casavisibilities - visibilities)/casavisibilities

    polarfig = figure()

    ang = arcsin(casavisibilities.flatten().imag)
    absq = casaabsq.flatten()
    todel = []
    for i in range(len(ang)):
        if math.isnan(abs(ang[i])):
            todel.append(i)
            continue
        polar([0, ang[i]],[0, absq[i]],marker='o')

    for i in range(len(todel)-1, -1, -1):
        ang = np.delete(ang, todel[i])
        absq = np.delete(absq, todel[i])

    title("Ratios of CASA Computed and Correlated Visibilities, 1.0 at 0 deg is perfect")
    # xlabel('Relative Error')
    # ylabel('Bin Count')

    savefig('VisRelDiff%i.png'%numsc)

    hist(ang, bins=20)

    title("Histogram of Visibility Ratios, 0 deg is perfect")
    #

    savefig('VisRelDiffHist%i.png'%numsc)

    plt.close(polarfig)

    #plot hist of absq
    histplot = figure()
    hist(absq, bins=20)
    title("Histogram of abs(Visibility) Ratios, 1.0 is perfect")
    savefig('VisAbsRatioHist%i.png'%numsc)
    plt.close(histplot)


#get RMSE

print 'calculating RMSE'
rmses = []

imregrid(imagename=RelicMergeMS+'.dirty'+str(l)+'.image', output="dirtyReshaped.image", template="DRAGN-10.000MHz.truth")

ia.close()
ia.open(imageName)
pix = ia.getchunk()
truth = pix.copy() #zpix
truth = truth.reshape((imSize, imSize))
# zoomimg = spndint.zoom(truth,float(imsize)/imSize)

for l in range(length):
    if gifMode == False and l != (length - 1):
        continue

    ia.close()

    imregrid(imagename=RelicMergeMS+'.dirty'+str(l)+'.image', output="dirtyReshaped.image", template=imageName, overwrite = True)

    ia.open('dirtyReshaped.image')
    pix = ia.getchunk()
    pix = pix.reshape((imSize, imSize))
    # p = np.flipud(pix.T)
    rmse = np.sqrt(((pix - truth)**2).mean())
    rmses.append(rmse)
    ia.close()

if gifMode == True:
    f = figure()
    plot(array(range(length))+1, rmses)
    title('RMSE as a function of Time with %i Spacecraft'%numsc)
    ylabel("RMSE")
    xlabel("Minutes used")
    savefig('rmses%i.png'%numsc)
    plt.close(f)


np.savetxt('rmses'+str(numsc)+'.txt', np.array(rmses).flatten(), fmt='%.7f')



if gifMode:

    GifCommand = "convert -delay 15 -loop 0 "+ ' '.join(map(str, jpgNames))+ " RELICMERGE-10.000MHz_time.gif"
    print 'to create gif, copy paste this command into terminal in output folder (requires imagemagick plugin for gif conversion)'
    print GifCommand


print 'done, now going back to ' +dir1
os.chdir(dir1)
