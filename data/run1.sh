#!/bin/bash

# set the path only for the script

PATH=../src/:$PATH

# download the data
wget http://pop-gen.eu/wordpress/wp-content/uploads/2014/05/proteinsJhaetal.tar1.gz

# decompress them

tar xvfz proteinsJhaetal.tar1.gz

# run prins

cd ./proteins_used_in_paper/

find `pwd` -iname '*.pdb' > list

cd ../

mkdir scores

cd scores

prins -list ../proteins_used_in_paper/list -oMatrix interMatrix.txt

cd ..; 
