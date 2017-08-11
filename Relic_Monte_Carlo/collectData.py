'''
Alex Hegedus 8/10/17
alexhege@umich.edu

This code is used to analyze the data left behind from RELIC trials

run with regular python, no arguments taken.

Structure of codes assumes you used genMC.py to generate truth images
and corrRelicMC.py to run the analysis

'''


import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from matplotlib.mlab import griddata




dir1=os.path.expandvars('$PWD')
print dir1

outFolder = os.path.join(dir1, 'outIm')
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

vals = {}

imageFolder = os.path.join(dir1, 'MCImages')
print imageFolder
os.chdir(imageFolder)

sc32succ = []
sc16succ = []
strs = [1.1, 1.3, 1.5, 2.0, 4.0]

for relStr in [1.1, 1.3, 1.5, 2.0, 4.0]:
    strfolder = os.path.join(imageFolder, 'relStr_'+str(relStr))
    print strfolder
    os.chdir(strfolder)

    vals[relStr] = {}

    for iters in [0, 100, 500, 1000]:
        vals[relStr][iters] = {}
        for numsc in [8, 16, 24, 32]:
            vals[relStr][iters][numsc] = {}
            vals[relStr][iters][numsc]["trueRA"] = []
            vals[relStr][iters][numsc]["trueDec"] = []
            vals[relStr][iters][numsc]["diff"] = []


            for runNum in range(100):

                runfolder = os.path.join(strfolder, 'run_'+str(runNum))
                print runfolder
                os.chdir(runfolder)

                SCfolder = os.path.join(runfolder, 'numSC_'+str(numsc))
                if not os.path.exists(SCfolder):
                    os.makedirs(SCfolder)
                print SCfolder
                os.chdir(SCfolder)




                posFile = os.path.join(SCfolder, 'posDiffIters%d.txt'%iters)

                try:
                    f = open(posFile, 'r')
                except IOError:
                    print 'missing ' + posFile
                    continue

                #true rad ra dec, dirty rad ra dec, diff rad, fitDiff rad, beamwidth/pixsize rad
                # [0.7862346124414664, -0.08673616854718207]
                # [0.7866535641839975, -0.09276642170878999]
                # [0.006044672293648232]
                # [0.005680357107761912]
                # [0.0006100000000000001]
                f.readline()
                trueCoord = f.readline().split()
                trueRA = float(trueCoord[0][1:-1])
                trueDec = float(trueCoord[1][:-1])

                f.readline()
                f.readline()
                fitDiff = float(f.readline()[1:-2])
                beamWidth = float(f.readline()[1:-2])

                vals[relStr][iters][numsc]["trueRA"].append(trueRA)
                vals[relStr][iters][numsc]["trueDec"].append(trueDec)
                vals[relStr][iters][numsc]["diff"].append(fitDiff)

            #now make plot out of that
            fig = plt.figure()
            #ax = fig.gca(projection='3d')
            # ax = Axes3D(fig)

            minRA = min(vals[relStr][iters][numsc]["trueRA"])
            minRA -= abs(.01*minRA)
            maxRA = max(vals[relStr][iters][numsc]["trueRA"])
            maxRA += abs(.01*maxRA)
            stepRA = abs(maxRA - minRA)/256.
            minDec = min(vals[relStr][iters][numsc]["trueDec"])
            minDec -= abs(.01*minDec)
            maxDec = max(vals[relStr][iters][numsc]["trueDec"])
            maxDec += abs(.01*maxDec)
            stepDec = abs(maxDec - minDec)/256.

            # Make data.
            X = np.arange(minRA, maxRA, stepRA)
            Y = np.arange(minDec, maxDec, stepDec)

            Z = griddata(vals[relStr][iters][numsc]["trueRA"], vals[relStr][iters][numsc]["trueDec"], vals[relStr][iters][numsc]["diff"], X, Y)

            X, Y = np.meshgrid(X, Y)

            plt.contourf(X, Y, Z)
            plt.colorbar()

            #make x and y axes same scale
            plt.axis('equal')

            os.chdir(dir1)

            plt.xlabel('RA location of truth hotspot (radians)')
            plt.ylabel('Dec location of truth hotspot (radians)')
            plt.title('Positional Difference in recovered feature position in radians\n RelStr = %2.2f, #S/C = %d, clean iters = %d'%(relStr, numsc, iters))

            plt.savefig('outIm/diffDist_relStr%2.2f_numsc%d_iters%d.png'%(relStr, numsc, iters))

            plt.close(fig)
            #Z = np.zeros(np.shape(X))


            fig = plt.figure()
            ax = Axes3D(fig)

            #scatter plot
            for i in range(len(vals[relStr][iters][numsc]["diff"])):
                currRA = vals[relStr][iters][numsc]["trueRA"][i]
                currDec = vals[relStr][iters][numsc]["trueDec"][i]

                ax.scatter(currRA, currDec, vals[relStr][iters][numsc]["diff"][i])

            plt.xlabel('RA location of truth hotspot (radians)')
            plt.ylabel('Dec location of truth hotspot (radians)')
            plt.title('Positional Difference in recovered feature position in radians\n RelStr = %2.2f, #S/C = %d, clean iters = %d'%(relStr, numsc, iters))

            #plt.tight_layout(pad=1.2, h_pad=1.2, w_pad=1.2)

            plt.savefig('outIm/ScatterdiffDist_relStr%2.2f_numsc%d_iters%d.png'%(relStr, numsc, iters), bbox_inches='tight', pad_inches=1)

            plt.close(fig)


            #success Ratios
            succ=0
            fail=0
            for i in range(len(vals[relStr][iters][numsc]["diff"])):
                if vals[relStr][iters][numsc]["diff"][i] <= beamWidth:
                    succ+=1
                else:
                    fail+=1

            if numsc == 16 and iters == 0:
                sc16succ.append(succ/100.)

            if numsc == 32 and iters == 0:
                sc32succ.append(succ/100.)

            fig = plt.figure()

            objects = ('Success', 'Failure')
            y_pos = np.arange(len(objects))
            performance = [succ, fail]

            plt.bar(y_pos[0], performance[0], align='center', alpha=0.5, color='blue')
            plt.bar(y_pos[1], performance[1], align='center', alpha=0.5, color='red')

            plt.xticks(y_pos, objects)

            plt.title("Number of trials where error < beamwidth %f rad \n RelStr = %2.2f, #S/C = %d, clean iters = %d"%(beamWidth, relStr, numsc, iters))


            plt.savefig('outIm/successBars_relStr%2.2f_numsc%d_iters%d.png'%(relStr, numsc, iters))

            plt.close(fig)


#plot scpacecraft over different relStr plots, success Rotate
fig = plt.figure()

xs = np.arange(len(sc32succ))

plt.bar(xs, sc16succ, align='center')

plt.xticks(xs, strs)

plt.title("Success Ratio for 16 Spacecraft")
plt.xlabel("Relative Strength of Feature (SNR)")
plt.ylabel("Success Ratio")
plt.savefig('outIm/successRatio_numsc%d.png'%(16))
plt.close(fig)

fig = plt.figure()

plt.bar(xs, sc32succ, align='center')

plt.xticks(xs, strs)

plt.title("Success Ratio for 32 Spacecraft")
plt.xlabel("Relative Strength of Feature (SNR)")
plt.ylabel("Success Ratio")
plt.savefig('outIm/successRatio_numsc%d.png'%(32))
plt.close(fig)
