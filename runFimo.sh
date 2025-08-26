#!/usr/bin/bash

module load meme

fas=$1
fimo --thresh 0.0001 --no-qvalue --text --verbosity 1 human.mouse.nonredund.final.meme $fas

