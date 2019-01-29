#!/usr/bin/env python
#The script to run over all model variations
#First make a soft link to the datafiles directory in the current one, like `ln -s  ../datafiles`

from math import sqrt

jobName='testNewNLOvfns'

def valErr(cnt, err1, err2 = -999):
    if err2 == -999:
        return [cnt, cnt+err1, cnt-err1]
    else:
        return [cnt, cnt+err1, cnt+err2]
#def valErr(cnt, err):

#For the template file fName replace the parameters with errors
def replace(fName, sh, sign): #0=cnt, 1=up, 2=dn
    with open(fName, "rt") as fin:
        tag = ''
        if sign == 1:
            tag = 'u'
        elif sign == 2:
            tag = 'd'
        fOutName = jobName +'/'+ fName+str(sh)+tag
        with open(fOutName, "wt") as fout:
            #print( sh - 1)
            for line in fin:
                for i,p in enumerate(pars):
                    if i == sh-1:
                        #print( 'radek', '@'+p[0], str(p[1][sign]))
                        line = line.replace('@'+p[0], str(p[1][sign]))
                    else:
                        line = line.replace('@'+p[0], str(p[1][0]))

                line = line.replace('@outFile', fOutName+'_dir'+'/'+'out'+'.root')

                fout.write(line)
        return fOutName


#For the template file fName replace the parameters with errors (all possible combinations)
def put2files(fName, pars):
    import os
    if os.path.isdir(jobName) == False:
        os.mkdir(jobName)
    inFiles = []
    inFiles.append(replace(fName, 0, 0))
    for i in range(len(pars)):
        inFiles.append(replace(fName, i+1, 1))
        inFiles.append(replace(fName, i+1, 2))
    return inFiles


pars= []
pars.append(('alphaS',   valErr(0.106, 0.002)))
#pars.append(('alphaS',   valErr(0.118, 0.002)))
pars.append(('mCharm',   valErr(1.4, 0.2)))
pars.append(('mBottom',  valErr(4.5, 0.5)))
pars.append(('a0_IR',   valErr(0.5,   0.1)))
pars.append(('ap_IR',   valErr(0.3,   0.6, -0.3)))
pars.append(('b0_IR',   valErr(1.6,   0.4, -1.6)))
pars.append(('ap_IP',   valErr(0.06, 0.19, -0.06)))
pars.append(('b0_IP',   valErr(5.5,   0.7, -2.0)))
pars.append(('Q0',    [sqrt(1.75), sqrt(2.05), sqrt(1.15)] ))
pars.append(('mu',    [1, 2, 0.5]))


#Create files with all variations
inFiles = put2files('H1diff_templ.str', pars)

#run all variations in parallel
procs = []
from subprocess import Popen
for fIn in inFiles:
    procs.append(Popen(['alpos', fIn]))
for p in procs: p.wait()
