#!/bin/bash

cd src/darcyBoundaryConditions/
wmake
cd ../densityModels
wmake
cd ../viscosityModels
wmake
cd ../../applications/solvers/AJRrhoDarcyFoam/
wmake
#cd ../AMRrhoDarcyFoam/
#wmake
