#!/bin/bash

# Clear ALL original traj files in /home/iarsenea/trajs/trajFiles.
# Use cautiously.

read -p "Do you wish to permanently empty all files in /trajFiles? (y/n) " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    exit 1
else
	echo "Deleting files..."
	rm /home/iarsenea/trajs/trajFiles/traj_*
fi
