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
