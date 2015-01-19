#!/bin/bash

cluster_name=$2

# Mount slaves now
for ((i=1;i<=$1;i++)); do
  echo ""
  echo "*****************run $i*********************"
  ssh root@bioworker$i "apt-get install glusterfs-server -y" 
done
