#!/bin/bash

(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvdc
(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvde
(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvdf
(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvdg
(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvdh
(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvdi
(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvdj
(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvdk
(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvdl
(echo n; echo p; echo 1; echo ; echo; echo w) | fdisk /dev/xvdm

mkdir /BIO/tmp
apt-get install -y lvm2

vgcreate xvdn /dev/xvdc1 /dev/xvde1 /dev/xvdf1 /dev/xvdg1 /dev/xvdh1 /dev/xvdi1 /dev/xvdj1 /dev/xvdk1 /dev/xvdl1 /dev/xvdm1
if [ "$?" == "0" ]; then
  echo "Created Volume Group xvdn successfully"
  echo ""
else
  echo "Failed to create Volume Group xvdn!"
  exit 1
fi

lvcreate -l 767990 -n xvdn1 xvdn
if [ "$?" == "0" ]; then
  echo "Created logical volume xvdn1 successfully"
  echo ""
else
  echo "Failed to create logical volume xvdn1!"
  exit 1
fi

mkfs.ext3 /dev/xvdn/xvdn1
if [ "$?" == "0" ]; then
  echo "Format xvdn1 successfully"
  echo ""
else
  echo "Failed to format xvdn1!"
  exit 1
fi

mkdir /share
mount /dev/xvdn/xvdn1 /share
mountpoint /share
if [ "$?" == "0" ]; then
  echo "Done!"
  echo ""
else
  echo "Failed to mount /share"
  echo ""
fi
