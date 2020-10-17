#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N trajJob2
#$ -o submit_out_traj2
#$ -e submit_error_traj2
#$ -q SERIAL1
#$ -pe fill 1

# Go to the proper location for the output files
# cd /home/iarsenea/trajs/mFiles

echo Starting traj job for submission 2 at $(date -u)


### Get the proper date to send to the python script ###
#echo $(date -u +'%m')
#echo $(date -u +'%d') 
#echo "12"

mm=$(date -u +'%m')
dd=$(date -u +'%d')
hh="12"



### Start the trajs using the proper environment ###
#time python /home/iarsenea/trajs/trajDriver_v.py 1 $mm $dd $hh 1 2 3 4
/home/iarsenea/miniconda3/envs/trajEnv/bin/python /home/iarsenea/trajs/trajDriver_v.py 22 $mm $dd $hh


### Clean up files ###
#mv /home/iarsenea/trajs/jobFile* /home/iarsenea/trajs/mFiles/.

### Give the proper permissions to the new files
chmod 777 /home/iarsenea/trajs/trajFiles/traj*
