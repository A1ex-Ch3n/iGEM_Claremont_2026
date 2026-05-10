#!/bin/bash
# Usage: bash scripts/watch_job.sh <jobid>
# Monitors a Laguna job — shows status + live log tail every 30 seconds

JOBID=${1:?"Usage: $0 <jobid>"}
LOGDIR="/project/jespinoza_1541/CChen29/iGEM_Claremont_2026/logs"
LOGFILE="$LOGDIR/boltz2_${JOBID}.out"

echo "Watching job $JOBID  (Ctrl+C to stop)"
echo "Log: $LOGFILE"
echo "============================================"

while true; do
    clear
    echo "=== $(date) ==="
    echo ""

    # Job status
    STATUS=$(squeue -j "$JOBID" 2>/dev/null | tail -1)
    if [ -z "$STATUS" ] || echo "$STATUS" | grep -q "JOBID"; then
        echo "Job $JOBID: COMPLETED or not found"
        echo ""
        echo "=== Final log tail ==="
        tail -30 "$LOGFILE" 2>/dev/null || echo "(no log file yet)"
        break
    else
        echo "Queue: $STATUS"
        echo ""
        echo "=== Log tail ==="
        tail -20 "$LOGFILE" 2>/dev/null || echo "(log not yet created — job may still be queuing)"
    fi

    sleep 30
done
