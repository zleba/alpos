#!/usr/bin/env python


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("super_steering", help="input super-steering file")
args = parser.parse_args()

if  'AExt' in args.super_steering:
   orders   = ['nlo', 'nnlo']
else:
   orders   = ['nloF']

#orders   = ['nloF', 'nlo', 'nnlo']
if 'FitA' in args.super_steering:
   dataSets = ['heraI']
else:
   dataSets = ['heraI', 'heraIjets', 'heraC', 'heraCjets', 'fps3D', 'fps4D']

#dataSets = ['heraI', 'heraC']


for Ord in orders:
    for data in dataSets:
        fName = args.super_steering[:-6] +'_'+ Ord +'_'+ data +'.str'
        with open(args.super_steering, "rt") as fin:
            with open('variants/'+fName, "wt") as fout:
                print fName
                for line in fin:
                    #Data replacement
                    if 'heraI' in data:
                        line = line.replace('$heraI', '')
                        line = line.replace('$heraC', '#')
                        line = line.replace('$fps3D', '#')
                        line = line.replace('$fps4D', '#')
                    elif 'heraC' in data:
                        line = line.replace('$heraI', '#')
                        line = line.replace('$heraC', '')
                        line = line.replace('$fps3D', '#')
                        line = line.replace('$fps4D', '#')
                    else:
                        line = line.replace('$heraI', '#')
                        line = line.replace('$heraC', '#')

                    if 'jets' in data:
                        line = line.replace('$jets', '')
                        line = line.replace('$fps3D', '#')
                        line = line.replace('$fps4D', '#')
                    else:
                        line = line.replace('$jets', '#')

                    if 'fps3D' in data:
                        line = line.replace('$fps3D', '')
                    else:
                        line = line.replace('$fps3D', '#')

                    if 'fps4D' in data:
                        line = line.replace('$fps4D', '')
                    else:
                        line = line.replace('$fps4D', '#')


                    #Order replacement
                    if 'nnlo' == Ord:
                        line = line.replace('$iOrd', '2')
                        line = line.replace('$scheme', 'FONLL-C')
                    else:
                        line = line.replace('$iOrd', '1')
                        line = line.replace('$scheme', 'FONLL-B')

                    if 'nloF' == Ord:
                        line = line.replace('$nf', '3')
                        line = line.replace('$asN ', '#')
                        line = line.replace('$asL ', '')
                    else:
                        line = line.replace('$nf', '0')
                        line = line.replace('$asN ', '')
                        line = line.replace('$asL ', '#')
                    fout.write(line)
        

