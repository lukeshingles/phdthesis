#!/usr/bin/env python3
from collections import OrderedDict
import os
import struct

yieldfilelist = []
for (dirpath, dirnames, filenames) in os.walk("shinglesetal2015onlinedata/"):
    yieldfilelist.extend(filter(lambda fn: fn[-4:]=='.txt', filenames))
    break

with open('yieldtables.tex', 'w') as fileout:
    fileout.write('\n')

    prevmass = '-1'
    for yieldfilename in yieldfilelist:
        mass = yieldfilename[10:].split('z')[0]
#        if mass != prevmass:
#            fileout.write(r'\section{' + mass + ' \Msun models}' + '\n')
        prevmass = mass
    
        modelproperties = []
        with open(dirpath + yieldfilename,'r') as fyields:
            print("reading " + yieldfilename)
            for i in range(3):
                modelproperties.append(fyields.readline().split())
            fyields.readline()
            
            #tablecaption = yieldfilename
            tablecaption = '{Yields for ' + mass + r' \Msun, $Z=0.0006$, $Y=' + modelproperties[0][10].rstrip(',') + '$'
            if mass == '3':
                mpmz = float(modelproperties[0][13])
                if mpmz > 1e-6:
                    tablecaption += ', $\M{pmz}=' + '{0:.3f}'.format(mpmz) + '$'
                else:
                    tablecaption += ', no PMZ.'
            tablecaption += '}'
#            fileout.write(r'\begin{longtable}{r r r r r r r}' + '\n')
#            fileout.write(r'\caption{' + tablecaption + '.}' + '\n')
            fileout.write(r'\clearpage\ctable[cap = ' + tablecaption + ', caption = ' + tablecaption + ', label = {tab:' + yieldfilename[:-4] + '}]{r r r r r r r}{}{' + '\n')
            fileout.write(r'\hline' + '\n')
            fileout.write(r'El&   $Z$&    log $\epsilon$(X)&[X/H]&   [X/Fe]&M$_\mathrm{yield}$&X$_\mathrm{yield}$\\' + '\n')
            fileout.write(r'\hline' + '\n')
            
            for line in fyields:
                if not line.startswith("#") and not line.startswith('   '): #and not line.startswith(' bi  83'):
                    row = line.split()
                    row[0] = row[0].capitalize()
                    if row[1] == '1':
                        row[0] = 'H'
                    fileout.write("&".join(row) + r'\\' + '\n')

                    if row[1] == '42':
                        fileout.write(r'\hline' + '\n')
                        fileout.write(r'}' + '\n')
                        fileout.write(r'\ctable[cap = , caption = ' + tablecaption + ', label = {tab:' + yieldfilename[:-4] + '-cont}, continued]{r r r r r r r}{}{' + '\n')
                        fileout.write(r'\hline' + '\n')
                        fileout.write(r'El&   $Z$&    log $\epsilon$(X)&[X/H]&   [X/Fe]&M$_\mathrm{yield}$&X$_\mathrm{yield}$\\' + '\n')
                        fileout.write(r'\hline' + '\n')

            fileout.write(r'\hline' + '\n')
            fileout.write(r'}' + '\n\n')
#            fileout.write(r'\end{longtable}' + '\n')
