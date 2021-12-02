#!/bin/bash
echo "Check CMSSW environment"
echo "$CMSSW_BASE"
echo "Will run:"
echo "$@"
eval "$@"
