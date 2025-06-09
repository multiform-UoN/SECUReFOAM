#!/bin/bash

cd src/darcyBoundaryConditions/
wmake
cd ../densityModels
wmake
cd ../viscosityModels
wmake
cd ../../applications/solvers/rhoDarcyFoam/
wmake
#cd ../AMRrhoDarcyFoam/
#wmake
