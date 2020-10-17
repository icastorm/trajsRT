#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N archJob1
#$ -o submit_out_arch1
#$ -e submit_error_arch1
#$ -q SERIAL1
#$ -pe fill 1

echo Starting archival job at $(date -u)

# Send current day's data to a tar file
bash /home/iarsenea/trajs/archiveTrajs.sh