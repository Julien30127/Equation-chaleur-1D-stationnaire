#######################################
# rjulien@Ark30.mk
# Default options for rjulien@Ark30 computer
#######################################
CC=gcc
LIBSLOCAL=-L/usr/lib -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include
OPTCLOCAL=-fPIC -march=native