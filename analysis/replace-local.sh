#!/bin/sh
read -p "This will REPLACE the local data store -- are you sure? " -r
if [[ $REPLY =~ ^[Yy]$ ]]; then
  aws s3 sync s3://phd-biosurv/ biosurv/ --delete --sse --storage-class REDUCED_REDUNDANCY $@
fi
