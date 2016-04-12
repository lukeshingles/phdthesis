#!/usr/bin/env python3
import struct
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from plotcommon import *

everynmodels = 1

timetp1 = []
mcoretp1 = []

evrootfolder = '/Users/lukes/Dropbox/Papers (first author)/2015 He-enhanced IMAGB Stars/generateplots/evolution/'

modelname = 'm3z0006y24'

# get time at TP1
with open(evrootfolder + modelname + "/iscvn.dat", mode='r') as file:
    for i in range(6):
        row = file.readline().split()
    row = file.readline().split()
    timetp1 = float(row[1])+34
    mintime = timetp1 - 50
    timeselectedtp = float(file.readline().split()[1])
    maxtime = timeselectedtp + 150

#listtimetp1[1] -= 2400

# read ev file
listTime = []
listceinner = []
listmhcore = []
listmhecore = []
listconvtime = [[]]
listconvinner = [[]]
listconvouter = [[]]
listL = []
listLH = []
listLHe = []
listLC = []
with open(evrootfolder + modelname + "/evall.dat", mode='rb') as file:
    print('Processing ' + modelname + "...")
    fileContent = file.read()
    print('Starting mass: %.3f' % struct.unpack("d", fileContent[4:12])[0])

    #b is the byte number for the current model
    for b in range(20,len(fileContent)-16,268*everynmodels): #268 means every model
#            nmod = struct.unpack("i", fileContent[b:b+4])[0]
        T = struct.unpack("d", fileContent[b+4:b+12])[0]
        Toffset = (T-timetp1) / 1.0e0 # scaled relative time

        if T <= maxtime:
            if T >= mintime:
                convinner = getVariable(fileContent,b,9)
                convouter = getVariable(fileContent,b,8)
                if (convouter - convinner) > 1e-5: # convection zone exists below the envelope
                    # end of one convection zone and start of another
                    if len(listconvinner[-1]) > 0 and ((convouter < listconvinner[-1][-1] or convinner > listconvouter[-1][-1])):
                        # when line plotting, close up top and bottom of convective zones at start and end
#                        listconvtime[-1].append(listconvtime[-1][-1])
#                        listconvinner[-1].append(listconvinner[-1][-1])
#                        listconvouter[-1].append(listconvinner[-1][-1])
#                        listconvtime.append([listTime[-1]])
#                        listconvinner.append([convinner])
#                        listconvouter.append([convinner])
                        listconvtime.append([])
                        listconvinner.append([])
                        listconvouter.append([])
                    listconvtime[-1].append(Toffset)
                    listconvinner[-1].append(convinner)
                    listconvouter[-1].append(convouter)
                elif len(listconvinner[-1]) > 0:
                    #signal end the current conv. zone by making a new empty one
                    listconvtime.append([])
                    listconvinner.append([])
                    listconvouter.append([])

                mtot = getVariable(fileContent,b,60)
#                listTime.append(T)
                listTime.append(Toffset)
                listceinner.append(getVariable(fileContent,b,7))
                listmhcore.append(getVariable(fileContent,b,5))
                listmhecore.append(getVariable(fileContent,b,10))
                listL.append(getVariable(fileContent,b,15))
                listLH.append(getVariable(fileContent,b,16))
                listLHe.append(getVariable(fileContent,b,17))
                listLC.append(getVariable(fileContent,b,18))
                
                # keep the base of the convective envelope fixed if it disappears
                if len(listceinner) > 1 and (mtot - listceinner[-1]) < 1e-2:
                    listceinner[-1] = listceinner[-2]
        else:
            break

assert len(listconvtime) == len(listconvinner) == len(listconvouter)
print("%d time samples" % len(listTime))

listtimeranges = [[listTime[0],listTime[980]],[listTime[980],listTime[-990]],[listTime[-990],listTime[-1]]]

fig, axes = plt.subplots(2, len(listtimeranges), figsize=(5,4), sharex=False, sharey=False, tight_layout={"pad":0.6, "w_pad":0.0, "h_pad":0.0})

modellabel = modelname[1] + ' M$_\odot$ $Y=0.' + modelname[-2:] + '$'

for nx in range(len(axes[0])): #column within the row
    for ny in range(len(axes)): #row
        plt.setp(axes[ny][nx].get_xticklabels(), fontsize=fsticklabel)
        plt.setp(axes[ny][nx].get_yticklabels(), fontsize=fsticklabel)
    #    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=2.0))
        axes[ny][nx].xaxis.set_major_locator(ticker.MultipleLocator(base=[25,10000,50][nx]))
        axes[ny][nx].set_xlim(xmin=listtimeranges[nx][0],xmax=listtimeranges[nx][1])
        axes[ny][nx].get_yaxis().set_label_coords(-0.12,0.5)
        axes[ny][nx].get_xaxis().get_major_formatter().set_useOffset(False)
        if nx != 0:
            plt.setp(axes[ny][nx].get_yticklabels(), visible=False)
        if ny != len(axes)-1:
            plt.setp(axes[ny][nx].get_xticklabels(), visible=False)
        for axis in ['top','bottom','left','right']:
            axes[ny][nx].spines[axis].set_linewidth(framewidth)

    axes[0][nx].set_yscale('linear')
    axes[0][nx].yaxis.set_minor_locator(ticker.MultipleLocator(base=0.1))
    axes[0][nx].set_ylim(ymin=min(listmhecore)*0.999,ymax=max(listceinner)*1.002)

    axes[0][nx].plot(listTime, listmhcore,  color='red', lw=1.0, label='H-exhausted core')
    axes[0][nx].plot(listTime, listmhecore, color='blue', lw=1.0, label='He-exhausted core')
    axes[0][nx].fill_between(listTime, listceinner, [max(listceinner)*2 for x in listceinner], facecolor=(0.55,0.80,0.64), interpolate=True, edgecolor='None')
    axes[0][nx].plot(listTime, listceinner, color=(0.09,0.5,0.16), markersize=markersize*0.5, marker='o', markeredgewidth=0, linestyle='None')

    for n in range(len(listconvinner)-1):
        axes[0][nx].fill_between(listconvtime[n], listconvinner[n], listconvouter[n], facecolor=(0.48,0.73,0.54), interpolate=True, edgecolor='None')
        axes[0][nx].plot(listconvtime[n], listconvinner[n], color=(0.09,0.5,0.16), markersize=markersize*0.5, marker='o', markeredgewidth=0, linestyle='None')
        axes[0][nx].plot(listconvtime[n], listconvouter[n], color=(0.09,0.5,0.16), markersize=markersize*0.5, marker='o', markeredgewidth=0, linestyle='None')

    axes[1][nx].set_yscale('log')
    axes[1][nx].set_ylim(ymin=1e0,ymax=max(listLHe)*1.1)
    labellist = ['$\mathrm{L}_\mathrm{H}$','$\mathrm{L}_\mathrm{He}$']
    for (yList,color,dashes,label) in zip([listLH,listLHe],['red','blue'],[(),(4,1)],labellist):
        axes[1][nx].plot(listTime, yList, color=color, marker='None', lw=linewidth, dashes=dashes, label=label)

axes[-1][1].set_xlabel(r'Time [yr]', fontsize=fs)
axes[-1][1].set_xticks(axes[-1][1].get_xticks()[1:-1])

axes[0][1].legend(loc='upper left',handlelength=1,frameon=False,numpoints=1,prop={'size':fs-2}, markerscale=3)
axes[0][0].set_ylabel(r'$\mathrm{M}_\mathrm{r} \,[\mathrm{M}_\odot]$', fontsize=fs)
axes[1][0].set_ylabel(r'$\mathrm{L} \,[\mathrm{L}_\odot]$', fontsize=fs)
axes[0][0].get_yaxis().set_label_coords(-0.26,0.5)
axes[1][0].get_yaxis().set_label_coords(-0.26,0.5)
axes[1][1].legend(loc='upper left',handlelength=1.5,frameon=False,numpoints=1,prop={'size':fs-1})

    #plt.locator_params(axis = 'y', nbins = 30)

fig.savefig('../chapter1/fig-tpcycle.pdf',format='pdf')
plt.close()
