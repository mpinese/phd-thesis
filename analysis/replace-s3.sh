#!/bin/sh
aws s3 sync biosurv/ s3://phd-biosurv/ --delete --sse --storage-class REDUCED_REDUNDANCY $@
