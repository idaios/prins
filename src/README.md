# This is my README

PRInS v0.01

Detecting interesting aminoacids through Residual-Interaction Statistics

PRInS is available under the GPLv3.0 and it's free. It works in Linux
systems and it has been tested in Ubuntu 12.04 64bit machines.

Compilation
--------------

To compile PRInS type:
make -f Makefile.gcc

This will create an executable file, named prins. Please copy and
paste this file into your PATH, e.g. $HOME/bin or /usr/bin. Then, you
will be able to execute PRInS from any location in your computer.

Examples 
-------------- 

We provide two main sets of examples. Please see the directory
'data'. There are two directories there, called set1 and set2. The
first, set1, contains about 80 protein structures and it has been used
in Jha et al paper. The second dataset, set2, contains over 1,000
protein structures.

To re-run the data yourself, you can do the following, assuming you
are located within the PRInS main directory:

cd data

# this will run PRInS for set1. Please open and read the run1.sh to
# see how the program has been called.  

./run1.sh

# this will run PRInS for set2. Please open and read the run1.sh to
# see how the program has been called.  

./run2.sh


