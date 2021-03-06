#!/usr/bin/env python
#The script to run over all model variations
#First make a soft link to the datafiles directory in the current one, like `ln -s  ../datafiles`

from math import sqrt

jobName='testPdfNLO'

def valErr(cnt, err1, err2 = -999):
    if err2 == -999:
        return [cnt, cnt+err1, cnt-err1]
    else:
        return [cnt, cnt+err1, cnt+err2]
#def valErr(cnt, err):

#For the template file fName replace the parameters with errors
def replace(fName, sh, shId, sign): #0=cnt, 1=up, 2=dn
    global pars
    with open(fName, "rt") as fin:
        #tag = ''
        tag = 2*shId + sign - 2 if shId != 0 else 0
        #if sign == 1:
        #    tag = 'u'
        #elif sign == 2:
        #    tag = 'd'
        fOutName = jobName +'/steering.str'+ str(tag) # fName
        with open(fOutName, "wt") as fout:
            #print( sh - 1)
            for line in fin:
                for p in pars:
                    if p in sh:
                        line = line.replace('@'+p, str(pars[p][sign]))
                    else:
                        line = line.replace('@'+p, str(pars[p][0]))
                import os
                line = line.replace('@outFile', os.environ['ALPOS_DIR']+'/farm/variants/'+ fOutName+'_dir'+'/'+'out'+'.root')

                fout.write(line)
        print fOutName
        return fOutName


#For the template file fName replace the parameters with errors
def parScan(fName, sh, sign): #0=cnt, 1=up, 2=dn
    global pars
    with open(fName, "rt") as fin:
        #tag = ''
        tag = sign
        #if sign == 1:
        #    tag = 'u'
        #elif sign == 2:
        #    tag = 'd'
        fOutName = jobName +'/steering.str'+ str(tag) # fName
        with open(fOutName, "wt") as fout:
            #print( sh - 1)
            for line in fin:
                for p in pars:
                    if p in sh:
                        line = line.replace('@'+p, str(pars[p][sign]))
                    else:
                        line = line.replace('@'+p, str(pars[p][0]))
                import os
                line = line.replace('@outFile', os.environ['ALPOS_DIR']+'/farm/variants/'+ fOutName+'_dir'+'/'+'out'+'.root')

                fout.write(line)
        print fOutName
        return fOutName







#For the template file fName replace the parameters with errors (all possible combinations)
def put2files(fName, shifts):
    import os
    if os.path.isdir(jobName) == False:
        os.mkdir(jobName)
    if os.path.isdir(jobName+'/logs') == False:
        os.mkdir(jobName+'/logs')
    inFiles = []
    inFiles.append(replace(fName, [], 0, 0))

    if len(shifts) > 1:
        for i in range(len(shifts)):
            inFiles.append(replace(fName, shifts[i], i+1, 1))
            inFiles.append(replace(fName, shifts[i], i+1, 2))
    else:
        sh = shifts[0][0]
        print pars
        for i in range(len(pars[sh])):
            inFiles.append(parScan(fName, shifts[0], i))
        


    return inFiles


pars = {} 
shifts = []

#import sys
#sys.exit(0)

#
##pars['alphaS']  =  valErr(0.118, 0.002)
#pars['alphaS'] =  valErr(0.106, 0.002)
#pars['mCharm']  =  valErr(1.4, 0.2)
#pars['mBottom'] =  valErr(4.5, 0.5)
#pars['a0_IR']   =  valErr(0.5,   0.1)
#pars['ap_IR']   =  valErr(0.3,   0.6, -0.3)
#pars['b0_IR']   =  valErr(1.6,   0.4, -1.6)
#pars['ap_IP']   =  valErr(0.06, 0.19, -0.06)
#pars['b0_IP']   =  valErr(5.5,   0.7, -2.0)
#pars['mu']      =  [1, 2, 0.5]
#pars['Q0']      =  [sqrt(1.75), sqrt(2.05), sqrt(1.15)] 
##pars['Q0']      =  [sqrt(2.5), sqrt(2.5), sqrt(2.5)] 
#
#shifts.append(['alphaS'])
#shifts.append(['mCharm'])
#shifts.append(['mBottom'])
#shifts.append(['a0_IR'])
#shifts.append(['ap_IR'])
#shifts.append(['b0_IR'])
#shifts.append(['ap_IP', 'b0_IP'])
#shifts.append(['mu'])
#shifts.append(['Q0'])


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--test", help="run small number of shifts", action="store_true")
parser.add_argument("-r", "--run", help="run in parallel the produced steerings", action="store_true")
parser.add_argument("steering_template", help="input steering file")
#parser.add_argument("output_directory", help="output folder")
args = parser.parse_args()
#print args.output


with open(args.steering_template, "rt") as fin:
    isIn = False
    for line in fin:
        if  '##Model systematics {{' == line.rstrip():
            isIn = True
        elif '##Model systematics }}' == line.rstrip():
            isIn = False
        if isIn:
            exec(line[1:])

print pars
print shifts



#just few shifts
if args.test:
    shifts = shifts[0:2]

jobName = args.steering_template + '_dir'




#Create files with all variations
inFiles = put2files(args.steering_template, shifts)
#run all variations in parallel
if args.run:
    procs = []
    from subprocess import Popen
    for fIn in inFiles:
        procs.append(Popen(['alpos', fIn]))
    for p in procs: p.wait()
