'''
Alex Hegedus 8/10/17
alexhege@umich.edu

This code is used to analyze the data left behind from SunRISE trials

run with casa, no arguments taken.

Structure of codes assumes you used multcasa.sh and corrCME.py to run analysis

'''

#update array file
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import casac
from taskinit import *
from casa import *
from casac import *




numTrials = 104

d1=os.path.expandvars('$PWD')

vals = {}

secs = ['0', '1920', '3840']
poss = ['startPos_0', 'startPos_32', 'startPos_64']
freqs = ['00.300', '00.311', '00.322', '00.345', 'uvfits_CME_CONCAT_0.3195MHz_',  \
        '08.000', '08.330', '08.660', '08.990', 'uvfits_CME_CONCAT_8.495MHz_', \
        '14.257', '14.780', '15.294', '16.406', 'uvfits_CME_CONCAT_15.18425MHz_']
numSC = [4,5,6]

for sec,pos in zip(secs, poss):

    posdir = os.path.join(d1,pos)
    vals[pos] = {}

    for num in numSC:

        numdir = os.path.join(posdir, 'numSC_'+str(num))
        vals[pos][num] = {}

        #skip where I dont have data, numsc 4 pos 64
        # if pos == poss[-1] and num == numSC[0]:
        #     continue

        for i in range(len(freqs)): # in freqs:

            freq = freqs[i]

            vals[pos][num][freq] = {}
            vals[pos][num][freq]['brightra'] = []
            vals[pos][num][freq]['brightdec'] = []
            vals[pos][num][freq]['truera'] = []
            vals[pos][num][freq]['truedec'] = []
            vals[pos][num][freq]['measra'] = []
            vals[pos][num][freq]['measdec'] = []
            vals[pos][num][freq]['diffra'] = []
            vals[pos][num][freq]['diffdec'] = []
            vals[pos][num][freq]['overalldiff'] = []
            #Adding new stuff
            vals[pos][num][freq]['overalldiffb'] = []
            vals[pos][num][freq]['majdiff'] = []
            vals[pos][num][freq]['mindiff'] = []
            vals[pos][num][freq]['padiff'] = []
            vals[pos][num][freq]['bdiffra'] = []
            vals[pos][num][freq]['bdiffdec'] = []
            vals[pos][num][freq]['truemaj'] = []
            vals[pos][num][freq]['truemin'] = []
            vals[pos][num][freq]['truepa'] = []

            for trial in range(1, numTrials+1): #13):
                trdir = os.path.join(numdir, 'trial_'+str(trial))

                if freq[:2] != 'uv': #fits_CME_CONCAT_8.36466666667MHz_60sec.ms.txt':
                    fdir = os.path.join(trdir, 'uvfits_CME-sim1-' + freq + 'MHz_time_' + sec + 'sec.ms.txt')
                else:
                    fdir = os.path.join(trdir, freq + sec + 'sec.ms.txt') #'uvfits_CME_CONCAT_1.01866666667MHz_60sec.ms.txt') #'uvfits_CME_CONCAT_8.36466666667MHz_60sec.ms.txt')



                try:
                    f = open(fdir, 'r')
                except IOError:
                    print 'missing ' + fdir
                    vals[pos][num][freq]['truera'].append(-1)
                    vals[pos][num][freq]['truedec'].append(-1)
                    vals[pos][num][freq]['truemaj'].append(-1)
                    vals[pos][num][freq]['truemin'].append(-1)
                    vals[pos][num][freq]['truepa'].append(-1)
                    continue


                #truthra, truthdec/-1,-1 ; truth maj min pa/-1, -1, -1; imfit dirty ra dec ; imfit dirty maj min pa ; bright ra dec ; true/bright ra dec ;  true diff/ 0 ; bright radec

                # truera,truedec = map(float, f.readline().split())
                # truemaj, truemin = map(float, f.readline().split())
                # truepa = map(float, f.readline().split())[0]
                measra, measdec = map(float, f.readline().split())
                measmaj, measmin = map(float, f.readline().split())
                measpa = map(float, f.readline().split())[0]
                # brightra,brightdec = map(float, f.readline().split())

                if freq[:2] =='uv': #fits_CME_CONCAT_8.36466666667MHz_60sec.ms.txt':
                    #print 'trial '+str(trial) + ' shape truera = ' + str(vals[pos][num][freqs[1]]['truera'])
                    truera = vals[pos][num][freqs[i-2]]['truera'][trial-1]
                    truedec = vals[pos][num][freqs[i-2]]['truedec'][trial-1]
                    truemaj = vals[pos][num][freqs[i-2]]['truemaj'][trial-1]
                    truemin = vals[pos][num][freqs[i-2]]['truemin'][trial-1]
                    truepa = vals[pos][num][freqs[i-2]]['truepa'][trial-1]

                    if truera == -1:
                        vals[pos][num][freq]['overalldiff'].append(-1)
                        continue


                    vals[pos][num][freq]['truera'].append(truera)
                    vals[pos][num][freq]['truedec'].append(truedec)
                    vals[pos][num][freq]['truemaj'].append(truemaj)
                    vals[pos][num][freq]['truemin'].append(truemin)
                    vals[pos][num][freq]['truepa'].append(truepa)

                    vals[pos][num][freq]['measra'].append(measra)
                    vals[pos][num][freq]['measdec'].append(measdec)

                    # bdiffra = abs(truera - brightra)
                    # bdiffdec = abs(truedec - brightdec)
                    diffra = abs(truera - measra)
                    diffdec = abs(truedec - measdec)

                    if truera > np.pi:
                        truera -= 2*np.pi
                        diffra = abs(truera - measra)

                    diffmaj = abs(truemaj - measmaj)
                    diffmin = abs(truemin - measmin)
                    diffpa = abs(truepa - measpa)

                    diff = np.sqrt((diffra*np.cos(truedec))**2 + diffdec**2)

                    #casadiff calc
                    #ang dist calc
                    tFeatDir = me.direction('J2000', str(truera)+'rad', str(truedec)+'rad')
                    mFeatDir = me.direction('J2000', str(measra)+'rad', str(measdec)+'rad')
                    diff = me.separation(tFeatDir, mFeatDir)['value']#*np.pi/180.


                    vals[pos][num][freq]['diffra'].append(diffra)
                    vals[pos][num][freq]['diffdec'].append(diffdec)
                    vals[pos][num][freq]['majdiff'].append(diffmaj)
                    vals[pos][num][freq]['mindiff'].append(diffmin)
                    vals[pos][num][freq]['padiff'].append(diffpa)

                    if diff < 1.:
                        vals[pos][num][freq]['overalldiff'].append(diff)

                    else:
                        vals[pos][num][freq]['overalldiff'].append(1.0)
                else:
                    vals[pos][num][freq]['truera'].append(measra)
                    vals[pos][num][freq]['truedec'].append(measdec)
                    vals[pos][num][freq]['truemaj'].append(measmaj)
                    vals[pos][num][freq]['truemin'].append(measmin)
                    vals[pos][num][freq]['truepa'].append(measpa)



allra = []
alldec = []

print 'now making figs'

#localization goal for .3 8 and 15 MHz in degrees

#assume cme is 2/3 wide as it is far from the sun
#want to localize to 1/4 that width

angDist = np.array([2.5, 0.75, 0.5])


successCutoff = angDist*2./3*1/4.

for sec,pos in zip(secs, poss):


    for k in range(len(freqs[4::5])):

        freq = freqs[4::5][k]
        cutOff = successCutoff[k]


        #plot all sc hist
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #plt.margins(0.1, 0.1)
        ax.tick_params(axis='x', pad=10)

        for j in range(len(numSC)):

            num = numSC[j]

            # #skip where I dont have data, numsc 4 pos 64
            # if pos == poss[-1] and num == numSC[0]:
            #     continue

            overalldiff = np.array(vals[pos][num][freq]['overalldiff'])
            #Histogram of positional differences

            bins = np.linspace(0., np.amax(overalldiff), max(int(100*np.amax(overalldiff)), 10))
            try:
                plt.hist(overalldiff, bins, alpha=0.5, label=str(num)+' Spacecraft') #, bins, range = (0, .05))
            except:
                print 'couldnt do hist on ' + str(pos) + ' ' + str(freq)


        plt.xlabel('absolute difference (degrees)')
        plt.ylabel('number of occurrences')

        if freq[:2] != 'uv': #fits_CME_CONCAT_8.36466666667MHz_60sec.ms.txt':
            plt.title('Positional difference histogram at time '+ sec+ ', \nfrequency ' + freq + 'MHz', y=1.08)
        else:
            #plt.title('Positional difference histogram at time '+ sec+ ', \nconcatenated frequency ' + freq[18:23] + 'MHz', y=1.08)
            plt.title('Positional difference histogram at pos '+ pos+ '\n Concatenated Frequency ' + freq[18:23] + 'MHz', y=1.05)

        ax.tick_params(axis='x', pad=10)

        plt.legend(loc='upper right')

        fig.savefig('scdiff_'+pos+'_'+freq+'MHz_hist.png', bbox_inches='tight', pad_inches=1) #, bbox_inches='tight', pad_inches=0.0)
        plt.clf()

        plt.close('all')

        #now do success plots

        succRatios = []
        for j in range(len(numSC)):

            num = numSC[j]

            #skip where I dont have data, numsc 4 pos 64
            # if pos == poss[-1] and num == numSC[0]:
            #     succRatios.append(0.)
            #     continue

            fig = plt.figure()

            objects = ('Success', 'Failure')
            y_pos = np.arange(len(objects))



            overalldiff = np.array(vals[pos][num][freq]['overalldiff'])

            #do spacecraft success plots
            succ=0
            fail=0
            for i in range(len(overalldiff)):
                if overalldiff[i] <= cutOff:
                    succ+=1
                else:
                    fail+=1

            succRatios.append(float(succ)/numTrials)



            performance = [succ, fail]

            plt.bar(y_pos[0], performance[0], align='center', alpha=0.5, color='blue')
            plt.bar(y_pos[1], performance[1], align='center', alpha=0.5, color='red')

            plt.xticks(y_pos, objects)

            plt.title("Number of trials where error < cutoff %2.4f deg \n #S/C = %d, Freq = %s MHz, Position = %s"%(cutOff, num, freq[18:23], pos))


            plt.savefig('successBars_%s_numsc%d_freq_%s.png'%(pos, num, freq[18:23]))

            plt.close(fig)


        #plot ratio over different numsc
        fig = plt.figure()

        plt.bar(numSC, succRatios, align='center')

        plt.title("Success Ratio over # Spacecraft \n Freq = %s MHz, Cutoff Error = %2.4f deg, Position = %s"%(freq[18:23], cutOff, pos))
        plt.xlabel("Number of Spacecraft Used")
        plt.ylabel("Success Ratio")
        plt.savefig('successRatio_%s_freq_%s.png'%(pos, freq[18:23]))
        plt.close(fig)








import json

json.dump(vals, open("data1.txt", 'w'))
