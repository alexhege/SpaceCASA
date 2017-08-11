#!/usr/bin/env python
################################################################################
# PROGRAM: CME.py
################################################################################

'''
Alex Hegedus 7/10/17
alexhege@umich.edu
This code runs a single instance of the monte carlo analysis for RELIC,
attempting to localize a feature on a disk

This code to be used with parallelRelic.sh

casa --nologger --nologfile -c corrRelicMC.py ${relStr} ${runNum} ${wantedSC}

'''
#import the required python modules
import matplotlib
import matplotlib.pyplot as plt
import time
import os,sys
import shutil
from collections import defaultdict

#import matplotlib.pyplot as plt
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
matplotlib.use('Agg')
import matplotlib.pyplot as plt





#call as casa corrRelicMS.py relStr runNum

if __name__ == "__main__":

    relStr = sys.argv[7]

    runNum = sys.argv[8]

    wantedSC = int(sys.argv[9])

    Telescope = 'VLA'
    scriptmode = True

    # for relStr in [1.1, 1.3, 1.5, 2.0, 4.0]:
    #     for runNum in range(5):
    #         for wantedSC in [16, 24, 32]:


    #Save all in folder
    dir1=os.path.expandvars('$PWD')
    print dir1

    orbitpath = os.path.join(dir1, 'inert_traj_061216_reconfig.txt')


    imageFolder = os.path.join(dir1, 'MCImages')
    print imageFolder
    os.chdir(imageFolder)

    strfolder = os.path.join(imageFolder, 'relStr_'+str(relStr))
    print strfolder
    os.chdir(strfolder)

    runfolder = os.path.join(strfolder, 'run_'+str(runNum))
    print runfolder
    os.chdir(runfolder)

    SCfolder = os.path.join(runfolder, 'numSC_'+str(wantedSC))
    if not os.path.exists(SCfolder):
        os.makedirs(SCfolder)
    print SCfolder
    os.chdir(SCfolder)



    ## import image
    #resolution
    freq = 10e6 #observing frequency in Hz
    imageName = '../DRAGN-%3.3fMHz.truth'%(freq/1e6)
    ia.open(imageName)
    pix = ia.getchunk()
    truth = pix.copy() #zpix
    imSize = shape(pix)[0]
    truth = pix.reshape((imSize, imSize))

    ia.close()


    ra = '03h00m00.0s'
    dec = '-05d00m00.0s'

    # ra = '19h59m28.3566s'
    # dec = '40d44m02.096s'

    antennaDiameter = 6.0 #meters

    imWidth = 400.*10 #arcseconds width of image

    orbitFile = orbitpath #'/Users/hegedus/Downloads/Relic_Pipeline/inert_traj_061216_reconfig.txt'
    #what the sed script looks to change to modify integration time, and time intervals used
    ranges =  ((0, 5)) #((2457690.5, 2457690.51))



    #width of image in radians
    width = imWidth/206265. # arcsec to radians #res/60. * np.pi/180.
    c = 3e8
    kwavnum = 2*np.pi*freq
    wavelen = c/freq
    kb = 1.380648e-23
    #calculate width of each pixel
    dres = width/imSize

    #position error in seconds
    posErr = True
    dTau = 5e-9
    poserr = 3e8*dTau #3. #3000. meters
    dv = poserr/wavelen
    dphase = 2*np.pi*dv


    phaseCenter = me.direction('J2000', ra, dec)

    zoomimg = truth



    pfile = open(orbitFile)
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

    ##take away random antenna arg 3 is numsc, so 32-32 is 0 full
    takeOff = numsc - wantedSC
    while takeOff != 0:
        ind = random.randint(0, numsc-1)
        allscpos = np.delete(allscpos, ind, axis=0)
        takeOff -= 1
        numsc -=1

    numbl = numsc*(numsc - 1)/2


    allscpos *= 1e2 #3 #convert to meters

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

    Ant_x = allscpos[:, 0, 0]
    Ant_y = allscpos[:, 0, 1]
    Ant_z = allscpos[:, 0, 2]

    print 'starting loop ' + str(length)
    #$

    RelicCurrMSs = []

    ia.open(imageName)
    csys = ia.coordsys()
    ia.close()
    for l in range(length):

        RelicCurrMS = 'RELIC-%3.3fMHz_time_'%(freq/1e6)+times[l]+'.ms'
        RelicCurrMSs.append(RelicCurrMS)

        sm.open(RelicCurrMS)
        refloc=me.observatory('VLA')

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
        rday1 += l/60./24 #add 1 min per length
        refday=me.epoch(rf='UTC',v0=str(rday1)+'d')
        sm.settimes(integrationtime='60s',usehourangle=False,referencetime=refday)

        #
        # Don't bother with autocorrelations
        sm.setauto(0.0);
        print 'Observing ...'
        sm.observe(sourcename='CME1',
                   spwname='HF',
                   observer='SunRISE',
                   starttime='0s', stoptime=str(1*60)+'s')
        # Fourier transform model image into data set



        #insert man made correlation visibilities here, from time series code NEWWWW


        ra = srcdir['m0']['value'] #returns radians
        dec = srcdir['m1']['value'] #returns radians

        res = imSize
        print 'res is ' + str(res)
    #real world sky width
    #assuming eash point is armin squared, as in CME code

    #https://casaguides.nrao.edu/index.php/N891_simdata_(CASA_3.4)


        x = np.cos(dec) * np.cos(ra)
        y = np.cos(dec) * np.sin(ra)
        z = np.sin(dec)
        s = np.array([x, y, z])


        #get x and y axes perp to s
        x = np.cross([0, 0, 1], s)
        y = np.cross(s, x)

        #normalize since cross product isn't already 1 for some reason
        # this is
        px = x/np.linalg.norm(x)
        py = y/np.linalg.norm(y)
        pz = s/np.linalg.norm(s)

        #numsc = len(Ant_x)

        #allscpos = np.array([Ant_x, Ant_y, Ant_z]).T #now 6x3



        #calculate new uvw either way
        # for l in range(length):

        for i in range(numsc):
            phasescpos[i][l][0] = np.dot(allscpos[i][l], px)
            phasescpos[i][l][1] = np.dot(allscpos[i][l], py)
            phasescpos[i][l][2] = np.dot(allscpos[i][l], pz)

        k=0
        for i in range(numsc):
          for j in range(i+1, numsc):
            baselines[k][l] =  phasescpos[j][l] - phasescpos[i][l]
            k+=1

        # correlate = True
        #
        # if correlate:


        print ' now on time ' + str(l)

        #which s/c is farthest down along pz?
        minZ = np.where(phasescpos[:,l,2] == min(phasescpos[:,l,2]))[0][0]

        #get delays in seconds, from diff in z from minZ
        delays = np.zeros(numsc)
        for i in range(numsc):
            delays[i] = phasescpos[i][l][2] - phasescpos[minZ][l][2]
            delays[i] = delays[i]/c

        #change later for more frequencies
        fourierCoeff = np.zeros(numsc, dtype=np.complex128)

        print np.shape(zoomimg)
        d1, d2 = np.shape(zoomimg)



        #loop over other pixels, adding in apporpriate delays for each sc fourierCoeff
        #like calculation of phase center, but t in front of variables for temporary, don't keep
        for i in range(d1):
            for j in range(d2):
                if zoomimg[i][j] == 0.0:
                    continue

                s = csys.toworld([i,j,0,0],'s')['string']

                tra = qa.convert(qa.toangle(s[0]), 'rad')['value']
                tdec = qa.convert(qa.toangle(s[1]), 'rad')['value']

                # print tra, tdec, i, j

                #unit direction pointing in direction of source/phase center
                tx = np.cos(tdec) * np.cos(tra)
                ty = np.cos(tdec) * np.sin(tra)
                tz = np.sin(tdec)


                ts = np.array([tx, ty, tz])


                # #get x and y axes perp to s

                tx = np.cross(ts, [0, 0, 1])
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




            ############## do baselines & cross correlation, multiply fourierCoeff

                    # k=0
                    # for i in range(numsc):
                    #   for j in range(i+1, numsc):
                    #     visibilities[k][l] =  fourierCoeff[j]*fourierCoeff[i].conj()
                    #     k+=1

                #load in visibilities into MS Data columns


            #sm.close()
    #$ loop it in!
        tb.open(RelicCurrMS, nomodify=False)



        actual_uvw = tb.getcol("UVW")

        baselines2 = np.reshape(baselines[:,l,:], (numbl,3), order='F').T


        tb.putcol("UVW", baselines2)

        print 'Predicting the visibilities ...'

        sm.predict(imagename=imageName)

            #exit()



            # if correlate:
            #sm.predict(imagename=imageName)
        data=tb.getcol("DATA")
        cdata=tb.getcol("CORRECTED_DATA")

        #"Now replacing CASA data with self made correlated data, comparing it"
        for i in range(numbl):
            # for j in range(l, l+1):
            #"check casa uvw and calculated"
            #print actual_uvw[:,i], baselines[:,i]

            q = data[0][0][i]/visibilities[i][l]

            #print data[0][0][i], visibilities[i][l], q, abs(q), q/abs(q)


            casaabsq[i][l] = abs(q)
            casavisibilities[i][l] = q/abs(q) #data[0][0][i]
            # if abs(q) == 0.:
            #     casavisibilities[i][l] = 0.


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

        #approx val of galactic noise
        Tgal = 9.e7*10.**(-2.477)
        ## add errors (noise, gains, polarization leakage, bandpass) to teh vsiibility data
        sm.setnoise(mode='simplenoise',simplenoise=str(Tgal)+'Jy')
        #sm.setnoise(mode='simplenoise',simplenoise='100Jy')


        print("corrupting the visibility data")
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

        print 'cleaning dirty image'

        #assume 600 km baseline is largest, lambda/D/4
        cellsizerad = wavelen/(600e3)/4.*10 #10 from spatial change from large baselines
        csas = cellsizerad*180./np.pi*3600

        beamwidth = 1.22*wavelen/(600e3)*10.

        #csas = cellsizerad

        imsize = int(imWidth/csas) #int(1024/60.*3600./csas)
        print 'pre imsize is ' + str(imsize)
        if imsize < 1:
            imsize = 1024
            print 'imsize used is ' + str(imsize)
        #csas = min(csas, 1.)
        imsize=256

        print 'cell size for ' + RelicMergeMS + ' is ' + str(csas) +' arcsec'

        for iters in [0, 100, 500, 1000]:
            try:
                clean(vis=RelicMergeMS,imagename=RelicMergeMS+'.'+str(iters)+'dirty'+str(l),
                    outlierfile='',
                    field='',spw='',
                    selectdata=False,
                    mode='mfs',nterms=1,
                    gridmode='widefield', wprojplanes = -1, facets = 1,
                    niter=iters,gain=0.1,threshold='0.0mJy',
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

                viewer(infile=RelicMergeMS+'.'+str(iters)+'dirty'+str(l)+'.image',displaytype='raster',
                      outfile='RELICMERGE-%3.3f_%ditersMHz_time_'%(freq/1e6, iters)+str(l)+'.jpg',outformat='jpg',
                      gui=False)
            except:
                print 'Clean wanted to crash on ' + RelicMergeMS
                pass


    #plot vis diff

    print 'now working on final figures\n'

    relErrors = abs(casavisibilities - visibilities)/casavisibilities
    # y = range(numbl)
    # x = range(length)
    # [X, Y] = np.meshgrid(x,y)
    #
    # cs=contourf(X, Y, relErrors)
    # plot(casavisibilities.flatten())

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
    xlabel("abs(casa vis/correlated vis)")
    savefig('VisAbsRatioHist%i.png'%numsc)
    plt.close(histplot)

    #get RMSE


    zoomimg = spndint.zoom(truth,float(imsize)/imSize)

    for iters in [0, 100, 500, 1000]:
        rmses = []
        for l in range(length):
            ia.close()
            ia.open(RelicMergeMS+'.'+str(iters)+'dirty'+str(l)+'.image')
            pix = ia.getchunk()
            pix = pix.reshape((imsize, imsize))
            p = np.flipud(pix.T)
            rmse = np.sqrt(((p - zoomimg)**2).mean())
            rmses.append(rmse)
            ia.close()


        fig1 = figure()
        plot(array(range(length))+1, rmses)
        title('RMSE as a function of Time with %i Spacecraft and cleaned for %i iterations'%(numsc, iters))
        ylabel("RMSE")
        xlabel("Minutes used")
        savefig('rmses%iiters%i.png'%(numsc, iters))

        plt.close(fig1)


        np.savetxt('rmses'+str(numsc)+'iters%d.txt'%iters, np.array(rmses).flatten(), fmt='%.7f')


        #now analyze pointing errors

        truthMaxPos = imstat(imageName)['maxposf'].split()
        maxRA = truthMaxPos[0][:-1]
        maxDec = truthMaxPos[1][:-1]

        tFeatDir = me.direction('J2000', maxRA, maxDec)

        maxRArad = tFeatDir['m0']['value']
        maxDecrad = tFeatDir['m1']['value']

        dirtyMaxPos = imstat(RelicMergeMS+'.'+str(iters)+'dirty'+str(length-1)+'.image')['maxposf'].split()
        dmaxRA = dirtyMaxPos[0][:-1]
        dmaxDec = dirtyMaxPos[1][:-1]

        dFeatDir = me.direction('J2000', dmaxRA, dmaxDec)

        dmaxRArad = dFeatDir['m0']['value']
        dmaxDecrad = dFeatDir['m1']['value']

        diff = me.separation(tFeatDir, dFeatDir)['value']*pi/180. #degrees to rad


        #do imfit
        xmaxp, ymaxp = imstat(RelicMergeMS+'.'+str(iters)+'dirty'+str(length-1)+'.image')['maxpos'][:2]
        boxhalfwidth=10
        bx = str(max(xmaxp - boxhalfwidth, 0)) + ', ' + str(max(ymaxp - boxhalfwidth, 0)) + ', ' + str(min(xmaxp + boxhalfwidth, imsize-1)) + ', ' + str(min(ymaxp + boxhalfwidth, imsize-1))
        ##addcomponent
        try:
            imfitout = imfit(imagename=RelicMergeMS+'.'+str(iters)+'dirty'+str(length-1)+'.image', box = bx)
            raimout = imfitout['results']['component0']['shape']['direction']['m0']['value']
            decimout = imfitout['results']['component0']['shape']['direction']['m1']['value']

        except:
            print 'failed imfit on truth image ' + RelicMergeMS+'.'+str(iters)+'dirty'+str(length-1)+'.image'
            raimout = dmaxRArad
            decimout = dmaxDecrad

        # raimout = imfitout['results']['component0']['shape']['direction']['m0']['value']
        # decimout = imfitout['results']['component0']['shape']['direction']['m1']['value']

        fitDir = me.direction('J2000', str(raimout)+'rad', str(decimout)+'rad')

        diffFit = me.separation(tFeatDir, fitDir)['value']*pi/180. #degrees to rad


        toWrite = defaultdict(list)
        toWrite[RelicMergeMS].append('true rad ra dec, dirty rad ra dec, diff rad, fitDiff rad, beamwidth/pixsize rad')
        toWrite[RelicMergeMS].append([maxRArad, maxDecrad])
        toWrite[RelicMergeMS].append([dmaxRArad, dmaxDecrad])
        toWrite[RelicMergeMS].append([diff])
        toWrite[RelicMergeMS].append([diffFit])
        toWrite[RelicMergeMS].append([beamwidth])

        f = open('posDiffIters%d.txt'%iters, 'w')
        for i in range(len(toWrite[RelicMergeMS])):
            f.write(str(toWrite[RelicMergeMS][i])+'\n')

        f.close()

        # np.savetxt('posDiff.txt', toWrite[RelicMergeMS])





    os.chdir(dir1)

                # convert -delay 15 -loop 0 $(for ((a=0; a<= 1; a+= 1)); do printf -- "%s_%s.jpg " RELICMERGE-10.000MHz_time $a; done;) RELICMERGE-10.000MHz_time.gif


                # import subprocess
                # subprocess.Popen("convert -delay 15 -loop 0 $(for ((a=0; a<= ; a+= 1)); do printf -- \" \%s_\%s.png \" $rootName $a; done;) ${rootName}.gif")
