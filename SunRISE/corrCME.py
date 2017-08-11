#!/usr/bin/env python
################################################################################
# PROGRAM: CME.py
################################################################################

'''
Alex Hegedus 8/10/17
alexhege@umich.edu
This simulation code has been put together to simulate CME emission for SunRISE mission

Read about parameters in argparse in code

instrution on how to run the program
within casa run this command:

%run corrCME.py -posFile SunRISE_relTrajs_sixSC.csv -componentFile components_new.dat -startPos ${startPos} -correlate True -numSC ${numsc}

or in terminal
casa --nologger --nologfile --nogui -c corrCME.py -posFile SunRISE_relTrajs_sixSC.csv -componentFile components_new.dat -startPos ${startPos} -correlate True -numSC ${numsc}
'''
#frame of reference is J2000/eme2000
#Unless you are doing milliarcsecond astronomy, you can ignore that bias and rotation.
#For most applications, J2000/FK5=EME2000=ICRF=ICRF2.

#import the required python modules
import time
import os,sys
import shutil
from collections import defaultdict
import matplotlib
#matplotlib.use('Agg')
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

from pylab import *
################################################################################
#defien input parameters
################################################################################
##define input options for running the code

def _getParser():
    argParser = argparse.ArgumentParser(description= __doc__,
                                        formatter_class=argparse. \
                                        RawDescriptionHelpFormatter )

    argParser.add_argument( '-posFile', required = True,
                            default = None, type = str, \
                            help = """ input posfile in ECI format """ )


    argParser.add_argument( '-componentFile', required = True, type=str,\
                              help="text file which contains the info to build the truth image.")

    argParser.add_argument( '-startPos', required = False,
                            default = 0, type = int,\
                            help = """ index into time of posFile to test """ )

    argParser.add_argument( '-correlate', required = False,
                            default = False, type = bool,\
                            help = """ #True: simulate manual correlation to form visibilities (slow)
                            # False: use CASA sm.predict to form visibilities (much faster)""" )

    argParser.add_argument( '-thermalNoise', required = False,
                            default = True, type = bool,\
                            help = """ turns on galactic thermal noise """ )

    argParser.add_argument( '-phaseNoise', required = False,
                            default = True, type = bool,\
                            help = """ turns on positional and clock bias noise """ )

    argParser.add_argument( '-numSC', required = False,
                            default = 6, type = int,\
                            help = """ randomly removes sc if less than 6 """ )

    argParser.add_argument( '--nologger', required = False,
                            action='store_true', \
                            help = """ For casa no log window""" )

    argParser.add_argument( '--nologfile', required = False,
                            action='store_true', \
                            help = """ For casa no log window""" )

    argParser.add_argument( '--nogui', required = False,
                            action='store_true', \
                            help = """ For casa no log window""" )

    argParser.add_argument( '-c', required = False,
                            action='store_true', \
                            help = """ For casa run as script""" )


    args = argParser.parse_known_args()

    return args[0]
################################################################################

def CMEinitial(datFile):
    '''
    This funciton retunrs a numpy array which initialized teh CME
    for each frequency
    '''

    c = []
    f = open(datFile, 'r')
    for line in f.readlines():
        line = line.split()
        c.append(map(float, line))

    c = np.array(c)

    #c[3:6, :] = components[15:18, :]



    return c#[:4, :]

    #Return 3 freqs around 8 MHz
    return components[15:18, :]

    #return every eighth to cut down on data
    return components[::8]



###############################################################################
################################################################################

def CMEform(comp):
    '''
    Form the truth image, We assume that the CME moves outward in frequency.

    '''
    ##
    delta = comp[1]
    delta1 = qa.toangle(str(delta)+'deg')
    offangle = comp[2]
    offangle1 = initoffangle #qa.toangle(str(offangle)+'deg')
    NewDir=me.shift(CMEDir, offset=delta1, pa=offangle1)
    major = comp[3]
    SizeMajor=qa.toangle(str(major)+'deg')
    minor = comp[4]
    SizeMinor=qa.toangle(str(minor)+ 'deg')
    PA=(90 + initoffangle['value'])%360. #comp[5]
    PA=qa.toangle(str(PA)+'deg')
    Flux=comp[6]
    resolution = comp[7]
    imSize = int(comp[8])

    delta2 = qa.toangle(str(delta)+'deg')
    offangle2 = qa.toangle(str(45)+'deg')
    NewDir2 = me.shift(CMEDir, offset=qa.toangle(str(2)+'deg'), pa=qa.toangle(str(280)+'deg'))
    #Construct an empty casa image from a shape
    ia.fromshape(imageName,shape=[imSize,imSize,1,1],overwrite=True)
    #adding components to  the empty image file
    cl.addcomponent(dir=NewDir,
                    flux=Flux/10, fluxunit='Jy', freq=comp[0],
                    shape="Gaussian", majoraxis=SizeMajor, minoraxis=SizeMinor, positionangle=PA)

    # cl.addcomponent(dir=CMEDir, flux=1e10, fluxunit='Jy', freq=comp[0],
    #                 shape='point')
    #
    # cl.addcomponent(dir=NewDir, flux=5e9, fluxunit='Jy', freq=comp[0],
    #                 shape='point')
    #
    # cl.addcomponent(dir=NewDir2, flux=7e9, fluxunit='Jy', freq=comp[0],
    #                 shape='point')

    cs=ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])
    cell_rad=qa.convert(qa.quantity("1.0arcmin"),"rad")['value']
    cs.setincrement([-cell_rad,cell_rad],'direction')
    cs.setreferencevalue([CMEDir['m0']['value'], CMEDir['m1']['value']],
                          type="direction")
    cs.setrestfrequency(comp[0])
    # add important header keywords
    imhead(imagename=imageName,mode="put",hdkey="object",hdvalue="Model CME")
    imhead(imagename=imageName,mode="put",hdkey="imtype",hdvalue='Intensity')
    imhead(imagename=imageName,
          mode="put",hdkey="observer",hdvalue="simulation")
    imhead(imagename=imageName,
          mode="put",hdkey="date-obs",hdvalue="2023/03/15/00:00:00")
    imhead(imagename=imageName,mode="put",hdkey="reffreqtype",hdvalue='TOPO')
    imhead(imagename=imageName,
           mode="put",hdkey="restfreq",hdvalue=str(comp[0]))
    imhead(imagename=imageName,mode='list')
    cs.setreferencevalue(str(comp[0])+'MHz', 'spectral')
    cs.settelescope('VLA')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(),subtract=False)

    pix = ia.getchunk()
    pix = pix.reshape((imSize, imSize))

    ia.close()
    return pix


################################################################################

def antPos(posfile):
    '''
    This function loads the snapshots of teh spacecrafts
    '''
    #antPos = np.loadtxt(posfile, skiprows=1, delimiter=',')

    antPos = np.loadtxt(posfile, skiprows=2, delimiter=',', usecols = (1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25))


    return antPos
################################################################################

def CMEsim(CMEms,comp,truthImage, img, Ant_x, Ant_y, Ant_z, args, dishSize=6, zBLErrs = None, clockBLErrs = None):
    '''
    This function will simulate the observations and creates an MS file
    '''

    correlate = args.correlate
    phaseNoise = args.phaseNoise
    thermalNoise = args.thermalNoise

    numsc = len(Ant_x)

    #one time things to compute
    dishDiameter = np.ones(numsc)*dishSize

    mounts = []
    antnames = []

    for i in range(numsc):
        mounts.append('X-Y')
        antnames.append('S'+str(i))

    sm.open(CMEms)
    refloc=me.observatory('VLA')
    #set the antenna configuration
    sm.setconfig(telescopename='VLA',
                 x=Ant_x,y=Ant_y,z=Ant_z,
                 dishdiameter=dishDiameter,
                 mount=mounts,
                 antname=antnames,
                 coordsystem='local',referencelocation=refloc)



    sm.setspwindow(spwname='HF',
                   freq=str(comp[0])+'MHz',
                   deltafreq='12200Hz', freqresolution='12200Hz', nchannels=1)
    # Set the feeds to be X Y.  Haven't figured it out.
    print 'Setting the feed polarization ...'
    sm.setfeed('perfect X Y')

    #
    print 'Setting the source ...'

    srcdir = CMEDir #me.direction('J2000', '00h00m00.0s', '00d00m00.0s')
    sm.setfield(sourcename='CME1', sourcedirection=srcdir); #CMEDir); ##Change

    #
    print 'Setting integration times ...'
    sm.settimes(integrationtime='.027s',usehourangle=False,referencetime=refday)

    #
    # Don't bother with autocorrelations
    sm.setauto(0.0);
    print 'Observing ...'
    sm.observe(sourcename='CME1',
               spwname='HF',
               observer='SunRISE',
               starttime='0s', stoptime='.027s')
    # Fourier transform model image into data set

    freq = comp[0]*1e6

    delta = comp[1]
    delta1 = qa.toangle(str(delta)+'deg')
    offangle = comp[2]
    offangle1 = initoffangle #qa.toangle(str(offangle)+'deg')
    NewDir=me.shift(CMEDir, offset=delta1, pa=offangle1)

    au = 1.496e11 #AU in m

    #compute baselines to put into ms

    ra = srcdir['m0']['value'] #returns radians
    dec = srcdir['m1']['value'] #returns radians

    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    s = np.array([x, y, z])

    #get x and y axes perp to s
    #get correct ones though!!!!

    x = np.cross(s, [0, 0, 1])
    y = np.cross(s, x)

    #normalize since cross product isn't already 1 for some reason
    # this is
    px = x/np.linalg.norm(x)
    py = y/np.linalg.norm(y)
    pz = s/np.linalg.norm(s)

    numsc = len(Ant_x)

    numbl = numsc*(numsc - 1)/2

    allscpos = np.array([Ant_x, Ant_y, Ant_z]).T #now 6x3

    phasescpos = np.zeros((numsc,3))

    for i in range(numsc):
        phasescpos[i][0] = np.dot(allscpos[i], px)
        phasescpos[i][1] = np.dot(allscpos[i], py)
        phasescpos[i][2] = np.dot(allscpos[i], pz)


    baselines = np.zeros((numbl, 3))
    k=0
    for i in range(numsc):
      for j in range(i+1, numsc):
        baselines[k] =  phasescpos[j] - phasescpos[i]
        k+=1

    baselines = baselines.T

    tb.open(CMEms, nomodify=False)
    actual_uvw = tb.getcol("UVW")
    tb.putcol("UVW", baselines)


    #insert man made correlation visibilities here, from time series code NEWWW

    if correlate:
        res = int(comp[8])
        print 'res ra dec ', res, ra, dec
        #real world sky width
        #assuming eash point is armin squared, as in CME code
        width = res/60. * np.pi/180.
        c = 3e8
        kwavnum = 2*np.pi*freq
        wavelen = c/freq
        kb = 1.380648e-23
        #calculate width of each pixel
        dres = width/res

        #make cs here to get ra decscs=ia.coordsys()
        ia.fromshape('temp.truth',shape=[res,res,1,1],overwrite=True)
        cs=ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])
        cell_rad=qa.convert(qa.quantity("1.0arcmin"),"rad")['value']
        cs.setincrement([-cell_rad,cell_rad],'direction')
        cs.setreferencevalue([CMEDir['m0']['value'], CMEDir['m1']['value']],
                              type="direction")
        cs.setrestfrequency(comp[0])
        cs.setreferencevalue(str(comp[0])+'MHz', 'spectral')

        ####

        #https://casaguides.nrao.edu/index.php/N891_simdata_(CASA_3.4)


        #which s/c is farthest down along pz?
        minZ = np.where(phasescpos[:,2] == min(phasescpos[:,2]))[0][0]

        #get delays in seconds, from diff in z from minZ
        delays = np.zeros(numsc)
        for i in range(numsc):
            delays[i] = phasescpos[i][2] - phasescpos[minZ][2]
            delays[i] = delays[i]/c

        #change later for more frequencies
        fourierCoeff = np.zeros(numsc, dtype=np.complex128)

        visibilities = np.zeros(numbl, dtype=np.complex128)

        print np.shape(img)
        d1, d2 = np.shape(img)

        #loop over other pixels, adding in apporpriate delays for each sc fourierCoeff
        #like calculation of phase center, but t in front of variables for temporary, don't keep
        for i in range(d1):
            for j in range(d2):
                if img[i][j] == 0.0:
                    continue

                #try casa way

                s = cs.toworld([i,j,0,0],'s')['string']

                tra = qa.convert(qa.toangle(s[0]), 'rad')['value']
                tdec = qa.convert(qa.toangle(s[1]), 'rad')['value']

                #print tra, tdec, i, j

                #unit direction pointing in direction of source/phase center
                tx = np.cos(tdec) * np.cos(tra)
                ty = np.cos(tdec) * np.sin(tra)
                tz = np.sin(tdec)


                ts = np.array([tx, ty, tz])


                # #get x and y axes perp to s

                tx = np.cross(ts, [0, 0, 1]) #, ts)
                ty = np.cross(ts, tx)

                #normalize since cross product isn't already 1 for some reason
                # this is
                tx = tx/np.linalg.norm(tx)
                ty = ty/np.linalg.norm(ty)
                tz = ts/np.linalg.norm(ts)

                tphasescpos = np.zeros((numsc,3))

                for k in range(numsc):
                    tphasescpos[k][0] = np.dot(allscpos[k], tx)
                    tphasescpos[k][1] = np.dot(allscpos[k], ty)
                    tphasescpos[k][2] = np.dot(allscpos[k], tz)

                #which s/c is farthest down along pz?
                tminZ = np.where(tphasescpos[:,2] == min(tphasescpos[:,2]))[0][0]

                tdelays = np.zeros(numsc)
                for k in range(numsc):
                    tdelays[k] = tphasescpos[k][2] - tphasescpos[tminZ][2]
                    tdelays[k] = tdelays[k]/c


                #phase of furthest down z in this projection, delay from this one
                phase = np.exp(1j*2*np.pi*np.random.random())
                C = img[i][j]*phase#*prefact
                for k in range(numsc):
                    fourierCoeff[k] = C*np.exp(-1j*kwavnum*(tdelays[k] - delays[k]))#*(1+(delta/360.*2*np.pi)))

                k=0
                for m in range(numsc):
                  for n in range(m+1, numsc):
                    toAdd = fourierCoeff[n]*fourierCoeff[m].conj()
                    toAdd = toAdd/np.sqrt(abs(toAdd))
                    visibilities[k] += toAdd
                    k+=1



        ############## do baselines & cross correlation, multiply fourierCoeff


        ia.close()


    print 'Predicting the visibilities ...'
    print("iFence['model'",comp[0])
    sm.predict(imagename=truthImage)
    #sm.predict(truthImage)

    qs = zeros(numbl, dtype=np.complex128)
    dqs = zeros(numbl, dtype=np.complex128)

    pure = zeros(numbl, dtype=np.complex128)

    if correlate:

        #prepare plots

        data=tb.getcol("DATA")
        cdata=tb.getcol("CORRECTED_DATA")

        absq = zeros(numbl)
        angq = zeros(numbl, dtype=np.complex128)



        ##

        print "Now replacing CASA data with self made correlated data, comparing it"
        for i in range(numbl):
            #"check casa uvw and calculated"
            #print actual_uvw[:,i], baselines[:,i]

            q = data[0][0][i]/visibilities[i]

            pure[i] = data[0][0][i]

            qs[i]=  q
            absq[i] = abs(q)
            angq[i] = arcsin(q.imag)

            #print data[0][0][i], visibilities[i], q, abs(q), q/abs(q)
            #print data[1][0][i], visibilities[i]

            data[0][0][i] = visibilities[i]
            data[1][0][i] = visibilities[i]

            cdata[0][0][i] = visibilities[i]
            cdata[1][0][i] = visibilities[i]


        tb.putcol("DATA", data)
        tb.putcol("CORRECTED_DATA", cdata)


    if phaseNoise == True:
        lambda1 = 3e8/freq

        data=tb.getcol("DATA")
        cdata=tb.getcol("CORRECTED_DATA")
        #"Now replacing CASA data with self made correlated data, comparing it"
        for i in range(numbl):
            dTau = zBLErrs[i] + clockBLErrs[i]
            poserr = 3e8*dTau #3. #3000. meters
            dv = poserr/lambda1
            dphase = 2*np.pi*dv
            #print dphase*180./np.pi
            #phaseErr = np.random.normal(0, dphase*np.sqrt(2./3)) #sqrt 2/3 : 2 from 1 for each sc pos, div by 3 to go from total pos error to only z direction
            data[0][0][i] *= np.exp(1.j*dphase)
            data[1][0][i] *= np.exp(1.j*dphase)
            cdata[0][0][i] *= np.exp(1.j*dphase)
            cdata[1][0][i] *= np.exp(1.j*dphase)

        tb.putcol("DATA", data)
        tb.putcol("CORRECTED_DATA", cdata)

        print("corrupting with phase Noise")

    if thermalNoise:

        Tgal = 9.e7*(freq/10e6)**(-2.477)
        tNoise = Tgal
        sm.setnoise(mode='simplenoise',simplenoise='1000000Jy')
        print("corrupting the visibility data with thermal noise")
        sm.corrupt()

    sm.close()


    if correlate:
        dvis=tb.getcol("DATA")[0][0][:]
        for i in range(numbl):
            dqs[i] = dvis[i]/pure[i]


    return (qs, dqs)


################################################################################

def CMEConcatModel(inputCME, outCMEms):
    '''
    this function concats nearby frequency channels
    inputCME: a list of CME files which will be concatenated
    outCME: the output file
    '''
    concat(vis=inputCME,
           concatvis=outCMEms,
           freqtol='',
           dirtol='',
           respectname=False,
           timesort=False,
           copypointing=True,
           visweightscale=[])



    return

################################################################################
def CMEuvmodelfit(MS, args, truth = True, qs = None, dqs = None, largestBL = 1.5e4):

    '''
    This function fits a uvmodel fit to a MS and compares the fit with actual
    values from the truth image.
    MS: the name of teh MS file
    sourceInitial: a text file which contains initial estimates of the component to be fitted
    truthPos: actual position of the component in the truth image

    '''
    print 'start of uvmodelfit'

    res = 1024

    #No breaky if not truth
    raimout = -1.
    decimout = -1.
    majimout = -1.
    minimout = -1.
    paimout = -1.
    diff = defaultdict(list)

    ##get the info for truth image

    if truth:
###
        xmax, ymax = imstat('CME-'+MS[9:18]+'.truth')['maxposf'].split(',')[:2]

        raTruth = qa.toangle(xmax)['value']*np.pi/180.
        decTruth = qa.toangle(ymax)['value']*np.pi/180.

        print imstat('CME-'+MS[9:18]+'.truth')['maxposf']

        ##Imstat truth stuff into log??
        CMEimstatst = imstat('CME-'+MS[9:18]+'.truth')
        xmaxp, ymaxp = CMEimstatst['maxpos'][:2]
        boxhalfwidth=100
        bx = str(max(xmaxp - boxhalfwidth, 0)) + ', ' + str(max(ymaxp - boxhalfwidth, 0)) + ', ' + str(min(xmaxp + boxhalfwidth, res-1)) + ', ' + str(min(ymaxp + boxhalfwidth, res-1))
        ##addcomponent
        try:
            imfitout = imfit(imagename='CME-'+MS[9:18]+'.truth', box = bx)
        except:
            print 'failed imfit on truth image ' + MS
            return

        raimout = imfitout['results']['component0']['shape']['direction']['m0']['value']
        decimout = imfitout['results']['component0']['shape']['direction']['m1']['value']
        majimout = imfitout['results']['component0']['shape']['majoraxis']['value']
        minimout = imfitout['results']['component0']['shape']['minoraxis']['value']
        paimout = imfitout['results']['component0']['shape']['positionangle']['value']

        ##ADDing imfitout
        diff[MS].append([raimout, decimout])
        diff[MS].append([majimout, minimout])
        diff[MS].append([paimout, 0.0])
        ###


    if not truth:
        #Hegedus add in initial pos prediction from brightest point in dirty image
        if args.correlate:
            #create baseline plots
            #plots
            print 'blah1'
            polarfig = figure()
            print 'figure created'

            absq = abs(qs.flatten())
            angq = arcsin(qs.flatten().imag)

            todel = []
            for i in range(len(angq)):
                if math.isnan(abs(angq[i])):
                    todel.append(i)
                    print i, ' wah bad angle'
                    continue
                polar([0, angq[i]],[0, absq[i]],marker='o')

            for i in range(len(todel)-1, -1, -1):
                angq = np.delete(angq, todel[i])
                absq = np.delete(absq, todel[i])

            title("Ratios of CASA Computed and Correlated Visibilities\n 1.0 at 0 deg is perfect")
            # xlabel('Relative Error')
            # ylabel('Bin Count')

            savefig('VisRelDiff%s.png'%MS, bbox_inches='tight', pad_inches=1)

            hist(angq, bins=4)

            title("Histogram of Visibility Ratios, 0 deg is perfect")
            #

            savefig('VisRelDiffHist%s.png'%MS)

            plt.close(polarfig)

            #plot hist of absq
            histplot = figure()
            hist(absq, bins=20)
            title("Histogram of abs(Visibility) Ratios, 1.0 is perfect")
            xlabel("abs(casa vis/correlated vis)")
            savefig('VisAbsRatioHist%s.png'%MS)
            plt.close(histplot)

            #
            #create dirty baseline plots
            #plots
            polarfig = figure()

            absq = abs(dqs.flatten())
            angq = arcsin(dqs.flatten().imag)

            todel = []
            for i in range(len(angq)):
                if math.isnan(abs(angq[i])):
                    todel.append(i)
                    print i, ' wah bad angle'
                    continue
                polar([0, angq[i]],[0, absq[i]],marker='o')

            for i in range(len(todel)-1, -1, -1):
                angq = np.delete(angq, todel[i])
                absq = np.delete(absq, todel[i])

            title("Ratios of CASA Computed and Noisy Correlated Visibilities \n 1.0 at 0 deg is perfect")
            # xlabel('Relative Error')
            # ylabel('Bin Count')

            savefig('NoisyVisRelDiff%s.png'%MS, bbox_inches='tight', pad_inches=1)

            hist(angq)

            title("Histogram of Noisy Visibility Ratios, 0 deg is perfect")
            #

            savefig('NoisyVisRelDiffHist%s.png'%MS)

            plt.close(polarfig)

            #plot hist of absq
            histplot = figure()
            hist(absq, bins=20)
            title("Histogram of Noisy abs(Visibility) Ratios, 1.0 is perfect")
            savefig('NoisyVisAbsRatioHist%s.png'%MS)
            plt.close(histplot)


        print 'cleaning dirty image'

        hz=float(MS[11:16])*10**6
        wavelen = 3e8/hz
        cellsizerad = wavelen/largestBL/16.
        csas = cellsizerad*180./np.pi*3600

        imWidth = 1024.*60 #arcsec, 1024 res built in

        #csas = cellsizerad

        imsize = int(imWidth/csas) #int(1024/60.*3600./csas)
        print 'pre imsize is ' + str(imsize)
        if imsize < 1:
            imsize = 1024

        #use next power of 2
        imsize = int(pow(2, math.ceil(math.log(imsize, 2))))/2
        print 'imsize used is ' + str(imsize)

        print 'cell size for ' + MS + ' is ' + str(csas) +' arcsec'

        try:
            clean(vis=MS,imagename=MS+'.dirty',
                outlierfile='',
                field='',spw='',
                selectdata=False,
                mode='mfs',nterms=1,
                gridmode='', #widefield', wprojplanes = -1, facets = 1,
                niter=0,gain=0.1,threshold='0.0mJy',
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
        except:
            print 'Clean wanted to crash on ' + MS
            pass

        print 'done cleaning, now getting brightest dirty point'




        try:
            CMEimstats=imstat(MS+'.dirty.image')
        except:
            print "Error, Cannot read image, exiting on " + str(MS)
            return

        if type(CMEimstats) == type(None) or type(CMEimstats) == type(True):
            print "Error Couldn't clean to image, exiting on " + MS
            return


        print CMEimstats['maxposf']

        #Get importuvfits
        xmaxp, ymaxp = CMEimstats['maxpos'][:2]
        print 'max x y of dirty is ' + str(xmaxp) +', '+ str(ymaxp)
        boxhalfwidth = int(imsize/40) #[25, 25, 10]
        # boxhalfwidth = boxhalfwidths[ind]
        print 'boxhalfwidth ', boxhalfwidth

        bx = str(max(xmaxp - boxhalfwidth, 0)) + ', ' + str(max(ymaxp - boxhalfwidth, 0)) + ', ' + str(min(xmaxp + boxhalfwidth, imsize-1)) + ', ' + str(min(ymaxp + boxhalfwidth, imsize-1))
        print 'box is ' + bx

        try:
            imfitout = imfit(imagename=MS+'.dirty.image', box = bx)
        except:
            print "Error, Cannot fit dirty image, exiting on " + str(MS)
            xmax, ymax = CMEimstats['maxposf'].split(',')[:2]

            raTruth = qa.toangle(xmax)['value']*np.pi/180.
            decTruth = qa.toangle(ymax)['value']*np.pi/180.
            diff[MS].append([raTruth, decTruth])
            diff[MS].append([-1., -1.])
            diff[MS].append([0.0, 0.0])
            np.savetxt('uvfits_'+MS+'.txt', diff[MS])

            return

        try:
            raimout = imfitout['results']['component0']['shape']['direction']['m0']['value']
        except:
            print "Error, Cannot fit dirty image get component 0 or output, exiting on " + str(MS)
            xmax, ymax = CMEimstats['maxposf'].split(',')[:2]

            raTruth = qa.toangle(xmax)['value']*np.pi/180.
            decTruth = qa.toangle(ymax)['value']*np.pi/180.
            diff[MS].append([raTruth, decTruth])
            diff[MS].append([-1., -1.])
            diff[MS].append([0.0, 0.0])
            np.savetxt('uvfits_'+MS+'.txt', diff[MS])
            return

        raimout = imfitout['results']['component0']['shape']['direction']['m0']['value']
        decimout = imfitout['results']['component0']['shape']['direction']['m1']['value']
        majimout = imfitout['results']['component0']['shape']['majoraxis']['value']
        minimout = imfitout['results']['component0']['shape']['minoraxis']['value']
        paimout = imfitout['results']['component0']['shape']['positionangle']['value']

        cl.close()
        diff[MS].append([raimout, decimout])
        diff[MS].append([majimout, minimout])
        diff[MS].append([paimout, 0.0])


    np.savetxt('uvfits_'+MS+'.txt', diff[MS])

    return

###########################################################################

if __name__ == "__main__":

    """
    how to run the code
    """

    refday=me.epoch(rf='UTC',v0='60017.83d')
    import argparse
    ##define the input argument
    args = _getParser()

    antPos = antPos(args.posFile)

    components = CMEinitial(args.componentFile)
    CMEDir = me.direction('J2000', '21h00m00.0s', '-20d00m00.0s')

    initoffangle = qa.toangle(str(360*np.random.random())+'deg')

    print 'initoffangle is ' + str(initoffangle)

    #Calculate phase noises
    #Calculate possible phase errors for this possible
    # baseline 1-2 (or 0-1) is mod 0 in the 10 long cycle of baselines.  Changing u of this by 3 meters
    #4.7 ns per sc pos error, clock error is 5.4 ns
    clockRMS = [4.8e-9, 2.6e-9, 4.5e-9, 5.4e-9, 4.4e-9, 3.7e-9]

    scPosErr = 4.7e-9
    scClockErr = 5.4e-9
    numSC = 6
    numbl = numSC*(numSC-1)/2
    zErrs = np.random.normal(0, scPosErr*np.sqrt(1./3), numSC)
    #clockErrs = np.random.normal(0, scClockErr, numSC)
    clockErrs = np.zeros(numSC)
    for i in range(numSC):
        clockErrs[i] = np.random.normal(0, clockRMS[i])

    zBLErrs = np.zeros(numSC*.5*(numSC-1))
    clockBLErrs = np.zeros(numSC*.5*(numSC-1))
    q=0
    for i in range(numSC):
        for k in range(i+1, numSC):
            zBLErrs[q] = zErrs[i] - zErrs[k]
            clockBLErrs[q] = clockErrs[i] - clockErrs[k]
            q += 1

    #close previosuly exising components
    cl.done()
    qa.done()
    me.done()
    countForm = 0

    pix_comp = []
    ## creating a truth image
    for i in range(len(components)):
        imageName = 'CME-%06.3fMHz.truth'%(components[i,0])
        print imageName
        cl.done()
        qa.done()
        me.done()
        pix = CMEform(components[i,:])
        pix_comp.append(pix)
        countForm = countForm +1

        outfile=os.path.join('.', imageName.strip('truth')+'jpg')
	#if the truth image already exists remove it

       # if os.path.exists(outfile):
     	#	os.remove(outfile)
        viewer(infile=imageName,displaytype='raster',
              outfile=outfile,outformat='jpg',
              gui=False)

    countSim = 0


    ## simulating the image
    ## first loop through the antPos file, to get the relative posiiton of teh antennas
    #for i in range(len(antPos)):
    startpos = args.startPos ###

    for j in [startpos]: ###range(1):
        ## antenna positions in m

        Ant_x = [antPos[j,4]*1E3 ,antPos[j,7]*1E3 ,antPos[j,10]*1E3 ,antPos[j,13]*1E3 , antPos[j,16]*1E3, antPos[j,19]*1E3]
        Ant_y = [antPos[j,5]*1E3 ,antPos[j,8]*1E3 ,antPos[j,11]*1E3 ,antPos[j,14]*1E3 , antPos[j,17]*1E3, antPos[j,20]*1E3]
        Ant_z = [antPos[j,6]*1E3 ,antPos[j,9]*1E3 ,antPos[j,12]*1E3 ,antPos[j,15]*1E3 , antPos[j,18]*1E3, antPos[j,21]*1E3]

        ##take away random antenna
        takeOff=0
        wantedSC = args.numSC
        if wantedSC > 0 and wantedSC < numSC:
            takeOff = int(numSC - wantedSC)

        while takeOff != 0:
            ind = np.random.randint(0, numSC-1)
            Ant_x = np.delete(Ant_x, ind)
            Ant_y = np.delete(Ant_y, ind)
            Ant_z = np.delete(Ant_z, ind)
            takeOff -= 1
            numSC -= 1

        numbl = numSC*(numSC-1)/2
        ## define the list of nearby frequency channels to image
        inputCME = []
        freq=[]
        quotients = zeros((4, numbl), dtype=np.complex128)
        dquotients = zeros((4, numbl), dtype=np.complex128)
        counter = 0
        for i in range(len(components)):

            truthImage = 'CME-%06.3fMHz.truth'%(components[i,0])

            sm.close()
            me.done()
            vp.reset()
            #start up the empty measurement set
            CMEms='CME-sim1-%06.3fMHz_time_%dsec.ms'%(components[i,0],antPos[j,0])
            print("components[i,0],antPos[j,0]", components[i,0],antPos[j,0])
            #call the CMEsim function to simulate the observations

            q, dq = CMEsim(CMEms,components[i,:],truthImage, pix_comp[i], Ant_x, Ant_y, Ant_z, args, dishSize=6, zBLErrs = zBLErrs, clockBLErrs = clockBLErrs)
            quotients[counter, :] = q
            dquotients[counter, :] = dq
	    #perform uvmodel fit on the MS file
            # try_one(CMEuvmodelfit(CMEms,args.CMEInitialPrefix,niter=10, compType='G'),20)
            CMEuvmodelfit(CMEms, args, truth = True)
            inputCME.append(CMEms)
            #concat three nearby frequency channels
            freq.append(components[i,0])
            counter+=1
            if len(inputCME) == 4 :
                freqVal = (freq[0] + freq[1] + freq[2] + freq[3])/4.
                print(freqVal)
                #if the concatenated file exists remove it
                outCMEms = 'CME_CONCAT_%sMHz_%dsec.ms'%(str(freqVal), antPos[j,0])
                #outCMEmscopy = 'CME_CONCAT_%sMHz_%dsec_control.ms'%(str(freqVal), antPos[i,0])
                if os.path.exists(outCMEms):
                    shutil.rmtree(outCMEms)

                countSim = countSim+1

                CMEConcatModel(inputCME, outCMEms)
                #CMEConcatModel(inputCME, outCMEmscopy)
                print 'entering uvmodelfit'
                CMEuvmodelfit(outCMEms, args, truth=False, qs = quotients, dqs = dquotients)
                print 'exiting uvm'
                counter = 0
                inputCME=[]
                freq = []


            sm.close()
