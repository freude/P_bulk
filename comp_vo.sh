#!/bin/bash

WDIR=$(pwd)

matlab -nojvm -nodisplay -nosplash -r "run('$WDIR/main_script0.m');exit;"

