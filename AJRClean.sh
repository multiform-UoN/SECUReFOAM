#!/bin/bash

cd src/darcyBoundaryConditions/
wclean
cd ../densityModels
wclean
cd ../viscosityModels
wclean
cd ../../applications/solvers/rhoDarcyFoam/
wclean
# cd ../AMRrhoDarcyFoam
# wclean
