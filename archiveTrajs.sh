#!/bin/bash

# Archives the current trajs and wrf files necessary to recreate the
# run on Quanah. Gets the names/locations of the files, creates links
# to them in the same location, tars them together, and sends it off
# to Quanah via ssh.

# Must be run before the next wrf run starts and after the plots are
# done being made and the long-term forecast has started.

# NOTE: When extracting the files, do so one at a time and carefully
# as there will be multiple files named the same thing. Extract one
# at a time and assign them new names to keep this from being an issue.


# Get the current date
mm=$(date -u +'%m')
dd=$(date -u +'%d')
hh="12"

# Define the members used
mems=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
       21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37
       38 39 40 41 42 )

# Instantiate the array of file names
arr=()


# namelists/DART stuff
arr+=( "/home/bancell/DART/lanai/models/wrf/membersV351HC/mem1/namelist.input" )
arr+=( "/home/bancell/DART/lanai/models/wrf/work/input.nml" )
arr+=( "/home/bancell/DART/lanai/models/wrf/work/prior_inflate_restart" )


# wrf.bdy
for mem in ${mems[@]}
do
  arr+=( "/home/bancell/DART/lanai/models/wrf/membersV351HC/mem${mem}/wrfbdy_d01" )
done


# wrf.input
for mem in ${mems[@]}
do
  arr+=( "/home/bancell/DART/lanai/models/wrf/membersV351HC/mem${mem}/wrfinput_d01" )
  arr+=( "/home/bancell/DART/lanai/models/wrf/membersV351HC/mem${mem}/wrfinput_d02" )
done


# Trajs from the links
for entry in "/home/iarsenea/trajs/trajFiles_links"/*
do
  arr+=( "$entry" )
done


# Plots
arr+=( "/home/iarsenea/trajs/targeting/2020${mm}${dd}${hh}/pdf_2020${mm}${dd}${hh}.pdf" )


# Targets
arr+=( "/home/iarsenea/trajs/targeting/2020${mm}${dd}${hh}/log_2020${mm}${dd}${hh}.txt" )


# r files (values and response box)
arr+=( "/home/bancell/enkfALT/src/targout/Rvalues_2020${mm}${dd}${hh}" )
arr+=( "/home/bancell/enkfALT/src/targout/subsetTARG_2020${mm}${dd}${hh}.txt" )


# Not time-sensitive DART stuff
arr+=( "/storage/RT_ENS/longsave/2020${mm}${dd}${hh}/obs_seq.out" )


# GEFS Files
#arr+=( "/storage/RT_longsave/2020${mm}${dd}06/" )


echo ${arr[@]}


# tar the files together
tar -chvzf /home/iarsenea/trajs/archive/trajBundle_${mm}${dd}${hh}.tar.gz ${arr[@]}
chmod 777 /home/iarsenea/trajs/archive/trajBundle_${mm}${dd}${hh}.tar.gz


# Delete temp files
rm /home/iarsenea/trajs/archive/temp/*


# Send the files in trajFiles to temp
mv /home/iarsenea/trajs/trajFiles/* /home/iarsenea/trajs/archive/temp/.