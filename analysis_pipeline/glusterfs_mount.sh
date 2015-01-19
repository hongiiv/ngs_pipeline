#!/bin/bash

for ((i=1;i<=$1;i++)); do
  #ssh bioworker$i "apt-get update"
  ssh bioworker$i "apt-get install -y glusterfs-server"

  ssh bioworker$i "mkdir /snu_wgs"
  ssh bioworker$i "mount -t glusterfs bioworker1:/dist_vol /snu_wgs"
  ssh bioworker$i "mountpoint /snu_wgs"
  if [ "$?" == "0" ]; then
    echo "********Mounted glusterfs /snu_wgs on bioworker$i successfully!**********"
    echo ""
  else
    echo "********Failed to mount glusterfs /snu_wgs on bioworker$i!**********"
    echo ""
  fi
done
