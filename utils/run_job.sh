#!/bin/bash
# run_job.sh <name> <command...> : run a job niced under nohup-friendly bookkeeping.
# Writes logs/jobs/<name>.{pid,log,start}, and on exit <name>.done (exit 0) or <name>.fail (exit code).
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
mkdir -p logs/jobs
name=$1; shift
rm -f logs/jobs/$name.done logs/jobs/$name.fail
date +%FT%T > logs/jobs/$name.start
echo $$ > logs/jobs/$name.pid
nice -n 15 "$@" > logs/jobs/$name.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then date +%FT%T > logs/jobs/$name.done; else echo "exit $rc $(date +%FT%T)" > logs/jobs/$name.fail; fi
