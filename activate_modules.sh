#!/bin/bash
# file to activate modules 


export mkProfile=/u/sw && source $mkProfile/etc/profile && module load gcc-glibc/5 && module load trilinos && module load muparser && module avail
#module list
