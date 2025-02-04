#!/bin/bash

LOG_FILE="examl_run.log"
echo "==== ExaML Run Log - $(date) ====" > "$LOG_FILE"

for i in {1..4}
do
    echo "[$(date)] Processing test${i}..." | tee -a "$LOG_FILE"
    START_TIME=$(date +%s)

    rm -f test${i}_partitions.binary
    rm -f RAxML_info.test${i}_partitions

    echo "[$(date)] Running parse-examl for test${i}..." | tee -a "$LOG_FILE"
    ./parser/parse-examl -s "./data/test${i}_results/test${i}_TRUE.phy" -q "./data/test${i}_results/test${i}_partitions.txt" -m DNA -n test${i}_partitions >> "$LOG_FILE" 2>&1

    mv test${i}_partitions.binary "./data/test${i}_results/"

    END_TIME=$(date +%s)
    ELAPSED_TIME=$((END_TIME - START_TIME))
    echo "[$(date)] Parsing completed for test${i} in $ELAPSED_TIME seconds." | tee -a "$LOG_FILE"
done

extract_time_from_log() {
    log_file="$1"

    # Extract the "Overall accumulated Time" field
    if [[ -f "$log_file" ]]; then
        exec_time=$(grep -oP 'Overall accumulated Time \(in case of restarts\): \K[\d.]+' "$log_file")
        if [[ -n "$exec_time" ]]; then
            echo "[TIME] ExaML execution took: ${exec_time} seconds." | tee -a "$LOG_FILE"
        else
            echo "[WARNING] Could not extract execution time from log file: $log_file" | tee -a "$LOG_FILE"
        fi
    else
        echo "[WARNING] Log file not found: $log_file" | tee -a "$LOG_FILE"
    fi
}

NUM_CORES=$(nproc)  # Detect available CPU cores
echo "[$(date)] Running ExaML with $NUM_CORES cores..." | tee -a "$LOG_FILE"

for i in {1..4}
do
    echo "[START] Running ExaML for test${i}..." | tee -a "$LOG_FILE"
    
    rm -f ExaML_info.test${i}_output

    START_TIME=$(date +%s)

    # Run ExaML and redirect output to log file
    mpirun.openmpi --use-hwthread-cpus -np "$NUM_CORES" ./examl/examl-AVX \
        -t "./data/test${i}_results/newick_trees.txt" \
        -m GAMMA \
        -s "./data/test${i}_results/test${i}_partitions.binary" \
        -n test${i}_output >> "$LOG_FILE" 2>&1

    END_TIME=$(date +%s)
    ELAPSED_TIME=$((END_TIME - START_TIME))

    # Extract execution time from log
    extract_time_from_log "./ExaML_info.test${i}_output"

    echo "[$(date)] ExaML completed for test${i} in $ELAPSED_TIME seconds." | tee -a "$LOG_FILE"
done

echo "[$(date)] ExaML execution completed." | tee -a "$LOG_FILE"
