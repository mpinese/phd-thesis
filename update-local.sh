#!/bin/sh
aws s3 sync s3://marpin-phd/ . --sse --storage-class REDUCED_REDUNDANCY $@
