#!/bin/bash
# $1 - number of slaves

if [ $# -eq 0 ]; then
  echo "Usage: ./nfs_mount.sh [num_slave]"
  exit 0
fi

# Mount master
echo "Mounting biomaster..."
ifup eth1 dhcp
dhclient eth1

mkdir /SNUH
mount -t nfs 172.16.0.251:/SNUH /SNUH
if mountpoint -q "/SNUH"; then
  echo "Done!"
  echo ""
else
  echo "**********Mounting biomaster failed!!!**********"
  echo ""
  exit 0
fi
 

# Mount slaves now
for ((i=1;i<=$1;i++)); do
  echo "Enable eth1 from bioworker$i..."
  ssh root@bioworker$i "ifup eth1 dhcp"
  ssh root@bioworker$i "dhclient eth1"

  echo "Mounting bioworker$i..."
  ssh root@bioworker$i "mkdir /SNUH"
  ssh root@bioworker$i "mount -t nfs 172.16.0.251:/SNUH /SNUH"

  ssh root@bioworker$i "mountpoint /SNUH"
  if [ "$?" == "0" ]; then
    echo "Done!"
    echo ""
  else
    echo "**********Mounting bioworker$i failed!!!**********"
    echo ""
  fi

done
