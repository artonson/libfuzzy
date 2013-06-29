#!/bin/bash

#MR_USER=tmp ./yamr_mc_fos_quality.py \
#-o output_probability.txt \
#-i 100000 \
#-j 1000 \
#-m probability \
#--alpha-num=100 \
#--beta-num=100 \
#2> ERR.TXT

MR_USER=tmp ./yamr_mc_fos_quality.py \
-r \
-o output_possibility.txt \
-i 1000000 \
-j 1000 \
-m possibility \
--alpha-num=100 \
--beta-num=100 

