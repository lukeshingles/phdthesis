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
        'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac')

highlightedelements = ['Fe','Rb','Y','Ba','Ce','Pb'] #'He','N','F',
highlightedatomicnumbers = [elements.index(x) for x in highlightedelements]
highlightedelementyposition = [0.0 for x in highlightedelements]

fig, ax = plt.subplots(1, 1, sharey=True, figsize=(5,3), tight_layout={"pad":0.2,"w_pad":0.0,"h_pad":0.0})

dppnspath = '/Users/lukes/Dropbox/Papers (first author)/2015 He-enhanced IMAGB Stars/generateplots/ppns320species/'

modelnames = ['m1.7z0006y24a4pmz001s320','m6z0006y24a0s320']
modellabels = ['1.7 M$_\odot$, Z=0.0006','6.0 M$_\odot$, Z=0.0006']
for (modelname,modellabel,color) in zip(modelnames,modellabels,['red','blue']):
    listatomicnumber = []
    listabundance = []
    with open(dppnspath + modelname + '/final.dat', mode='r') as file:
        for line in file.readlines():
            row = line.split()
            if len(row) > 1 and not line.startswith('#'):
                abundance = float(row[4])
                atomicnumber = int(row[1])
                if abundance > -5 and atomicnumber != 84:
                    listatomicnumber.append(atomicnumber)
                    listabundance.append(abundance)
                    for highlightedatomicnumber in highlightedatomicnumbers:
                        if abs(atomicnumber-highlightedatomicnumber) <= 1: #is the element or a close neighbour
                            index = highlightedatomicnumbers.index(highlightedatomicnumber)
                            if abundance > highlightedelementyposition[index]:
                                highlightedelementyposition[index] = abundance
    
    tcindex = listatomicnumber.index(elements.index('Tc'))
    pmindex = listatomicnumber.index(elements.index('Pm'))
    # skip Tc and Pm
    for ((startindex,stopindex)) in [[0,tcindex],[tcindex+1,pmindex],[pmindex+1,len(listatomicnumber)]]:
        ax.plot(listatomicnumber[startindex:stopindex], listabundance[startindex:stopindex], marker='o', markersize=markersize*1.5, markeredgewidth=0, color=color, lw=linewidth, label=(modellabel if startindex==0 else ''))

plt.setp(plt.getp(ax, 'xticklabels'), fontsize=fsticklabel)
plt.setp(plt.getp(ax, 'yticklabels'), fontsize=fsticklabel)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(framewidth)

for (x,y,symbol) in zip(highlightedatomicnumbers,highlightedelementyposition,highlightedelements):
    ax.annotate(symbol, xy=(x, y - 0.0 * (x % 2)), xycoords='data', textcoords='offset points', xytext=(0,10), horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=fs-1.5)

#ax.set_yscale('log')
ax.set_xlabel('Atomic number', fontsize=fs)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=5))
ax.set_xlim(xmin=25,xmax=85)

ax.set_ylabel(r'[X/Fe]', fontsize=fs)
ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.5))
ax.set_ylim(ymin=-0.1,ymax=3.4)
ax.legend(loc='upper left',handlelength=2,frameon=False,numpoints=1,prop={'size':fs-1})

fig.savefig('../chapter1/fig-agbfinalabund.pdf',format='pdf')
plt.close()