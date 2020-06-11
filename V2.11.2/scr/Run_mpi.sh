#!/bin/sh

if [ $# -eq 1 ]
then
  cp $1 defruns
fi
mpirun ./coherens
