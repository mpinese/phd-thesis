#!/bin/sh
aws s3 sync s3://phd-biosurv/ biosurv/ --sse --storage-class REDUCED_REDUNDANCY $@
