#!/bin/bash

#bash killall.sh
rm -rf startPos*
nohup bash multcasa.sh > multout.out &
