# --------------------------------------------------------------
# GNUmakefile
# --------------------------------------------------------------

name := Darwin4.0
G4TARGET := $(name)
G4EXLIB := true

G4DEBUG := 0

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)
ROOTGLIBS       = $(shell root-config --glibs)

EXTRALIBS +=$(ROOTLIBS)
EXTRALIBS +=$(ROOTGLIBS)
EXTRALIBS += $(CLHEPLIBS)


CPPFLAGS += $(ROOTCFLAGS) #-g -g3 -ggdb
#CPPFLAGS += $(ROOTCFLAGS) -g -g3 -ggdb

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
