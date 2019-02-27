#!/usr/bin/env python


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("super_steering", help="input super-steering file")
args = parser.parse_args()

#orders   = ['nloF', 'nlo', 'nnlo']
dataSets = ['heraI', 'heraIjets', 'heraC', 'heraCjets']

orders   = ['nlo', 'nnlo']
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
                    elif 'heraC' in data:
                        line = line.replace('$heraI', '#')
                        line = line.replace('$heraC', '')
                    if 'jets' in data:
                        line = line.replace('$jets', '')
                    else:
                        line = line.replace('$jets', '#')

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
        
