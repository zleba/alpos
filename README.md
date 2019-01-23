[![Build Status](https://travis-ci.org/zleba/alpos.svg?branch=master)](https://travis-ci.org/zleba/alpos)

# Installation
For a machine with cvmfs access 
1. `git clone git@github.com:zleba/alpos.git`
2. `cd alpos`
3. `. ./setup.sh`
4. `./installQCDnum.sh`
5. `./installFast.sh`
6. `./installAlpos.sh`

To run diffractive fit:
1. `cd alpos`
2. `alpos steering/H1diff.str`

Make sure that the environment was set by `. ./setup.sh` before.  
This will initialize ROOT and LHAPDF, so they don't need to be install.
The `setup.sh` also initialize the new gcc compiler.

In case, that some changes were done in alpos, please recompile like
1. `cd alposBuild`
2. `make`

Notice, that if new files were added, one needs to call `cmake` again, for example using `./installAlpos.sh`.
