#!/bin/bash
# overnight_monitor.sh: every 10 min append a health snapshot to logs/overnight_status.log until all jobs
# listed in logs/jobs/*.start have a .done or .fail flag (or 16 h elapsed).
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
S=logs/overnight_status.log
for i in $(seq 1 96); do
  {
    echo "===== $(date +%FT%T)  $(uptime | sed 's/.*load/load/')  mem_avail=$(free -g | awk '/Mem/{print $7}')G  disk_free=$(df -h ~ | awk 'NR==2{print $4}')"
    open=0
    for st in logs/jobs/*.start; do
      n=$(basename $st .start)
      if [ -f logs/jobs/$n.done ]; then state="DONE $(cat logs/jobs/$n.done)"
      elif [ -f logs/jobs/$n.fail ]; then state="FAIL $(cat logs/jobs/$n.fail)"
      elif kill -0 $(cat logs/jobs/$n.pid) 2>/dev/null; then state="running since $(cat $st)"; open=$((open+1))
      else state="DIED (no flag, pid gone)"; fi
      err=$(grep -c -i -E "traceback|error|killed" logs/jobs/$n.log 2>/dev/null)
      echo "  $n: $state | errors=$err | last: $(tail -c 300 logs/jobs/$n.log 2>/dev/null | tail -1 | cut -c1-120)"
    done
    echo "  our processes: $(ps -u $USER -o pcpu= | awk '$1>5' | wc -l) busy, RSS $(ps -u $USER -o rss= | awk '{s+=$1} END {printf "%.1f", s/1e6}') GB"
  } >> $S
  [ $open -eq 0 ] && [ $i -gt 1 ] && { echo "===== all jobs finished $(date +%FT%T)" >> $S; exit 0; }
  sleep 600
done
