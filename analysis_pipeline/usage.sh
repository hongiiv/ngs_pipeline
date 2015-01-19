#!/bin/bash

cluster_name=$1

# log function
LOG_FILE="/SNUH/SNUH_2/scripts/usage/$cluster_name/usage.log"
LOG() {
  echo "[$(date "+%F %T")][$cluster_name] $*" >> ${LOG_FILE}
}

while true; do
    LOG "`qhost`"
    sleep 600
done
