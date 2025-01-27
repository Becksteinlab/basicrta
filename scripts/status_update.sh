#!/usr/bin/bash                                                                 

if [ -z "$1" ]; then 
	echo "Provide cutoff"
	exit 1
fi
cutoff=$1

sacct -n -X --format jobname,jobid -s RUNNING > .running_jobs.txt               
filenames=`awk -v col1=1 -v col2=2 '{printf "basicrta-%s/%s/slurm_%s.out \n", "'$cutoff'", $col1, $col2}' < .running_jobs.txt`

for file in $filenames; do
	if [ -e $file ]; then
		echo $file | xargs tail -n 1 | sed 's@[^a-z A-Z0-9.:%|</]@@;s/\].*/\]/'; echo
	fi;
done
