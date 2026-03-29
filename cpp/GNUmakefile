## path to the 'Chombo/lib' directory
CHOMBO_HOME := /home/zkeheng/Chombo-3.2/lib

## defines the Chombo variables
include $(CHOMBO_HOME)/mk/Make.defs

## this is the name of the file with 'main()'
ebase := main

## Chombo libraries needed by this program
LibNames := AMRElliptic AMRTools BoxTools

## Other libraries needed by this program (using -L and -l options)
XTRALIBFLAGS :=

## ’all’ is the default target and ’all-test’ is defined in Make.rules
all: all-test

include $(CHOMBO_HOME)/mk/Make.rules
