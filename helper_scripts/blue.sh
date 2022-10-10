#!/bin/bash

PDB=$1

/software/rosetta/latest/bin/make_blueprint.hdf5.linuxgccrelease -s $PDB   -blue $PDB.bp  #-mute all

