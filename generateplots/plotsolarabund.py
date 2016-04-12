#!/usr/bin/env python3
import math
import struct
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from plotcommon import *

elements = ('','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P',
        'S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu',
        'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru',
        'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La',
        'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
        'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At',
        'Rn','Fr','Ra','Ac','Th','Pa','U')

listatomicnumber = []
listabundance = []

highlightedelements = ['H','C','O','Fe','Sr','Ba','Pb','Th']
highlightedatomicnumbers = [elements.index(x) for x in highlightedelements]
highlightedelementyposition = [0.0 for x in highlightedelements]

fig, ax = plt.subplots(1, 1, sharey=True, figsize=(5,3), tight_layout={"pad":0.2,"w_pad":0.0,"h_pad":0.0})

with open('solar_abund_ARAAedit.dat', mode='r') as file:
    for line in file.readlines():
        row = line.split()
        if len(row) > 1:
            abundance = float(row[3])
            atomicnumber = int(row[0])
#            if abundance > -5:
            listatomicnumber.append(atomicnumber)
            listabundance.append(abundance)

            for highlightedatomicnumber in highlightedatomicnumbers:
                if abs(atomicnumber-highlightedatomicnumber) <= 1: #is the element or a close neighbour
                    index = highlightedatomicnumbers.index(highlightedatomicnumber)
                    if abundance > highlightedelementyposition[index]:
                        highlightedelementyposition[index] = abundance

tcindex = listatomicnumber.index(43) #Tc
pmindex = listatomicnumber.index(61) #Pm
# skip Tc and Pm
for ((startindex,stopindex)) in [[0,tcindex],[tcindex+1,pmindex],[pmindex+1,len(listatomicnumber)]]:
    ax.plot(listatomicnumber[startindex:stopindex], listabundance[startindex:stopindex], marker='o', markersize=markersize*1.5, markeredgewidth=0, color='blue', lw=linewidth)

plt.setp(plt.getp(ax, 'xticklabels'), fontsize=fsticklabel)
plt.setp(plt.getp(ax, 'yticklabels'), fontsize=fsticklabel)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(framewidth)
#    ax.annotate(modellabel, xy=(0.05, 0.93), xycoords='axes fraction', horizontalalignment='left', verticalalignment='top', fontsize=fs)
    #ax.set_yscale('log')

for (x,y,symbol) in zip(highlightedatomicnumbers,highlightedelementyposition,highlightedelements):
    ax.annotate(symbol, xy=(x, y - 0.0 * (x % 2)), xycoords='data', textcoords='offset points', xytext=(0,10), horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=fs-1.5)

ax.set_xlabel("Atomic number", fontsize=fs)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=5))
ax.set_xlim(xmin=0,xmax=95)

ax.set_ylabel(r'log(X/H) + 12', fontsize=fs)
ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.5))
ax.set_ylim(ymin=-0.5,ymax=13.5)

fig.savefig('../chapter1/fig-solarabundances.pdf',format='pdf')
plt.close()