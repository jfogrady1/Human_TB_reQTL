#!/usr/bin/env bash
set -eu

# Check status of Slurm job

jobid=$(echo $1 | awk '{print $(NF)}')

if [[ "$jobid" == Submitted ]]
then
  echo smk-simple-slurm: Invalid job ID: "$jobid" >&2
  echo smk-simple-slurm: Did you remember to add the flag --parsable to your sbatch call? >&2
  exit 1
fi

output=$(scontrol show jobid -dd "$jobid" | grep -m 1 -oP "JobState=\w*" | cut -c 10-)

if [[ $output =~ ^(COMPLETED).* ]]
then
  echo success
elif [[ $output =~ ^(RUNNING|PENDING|COMPLETING|CONFIGURING|SUSPENDED).* ]]
then
  echo running
else
  echo failed
fi
