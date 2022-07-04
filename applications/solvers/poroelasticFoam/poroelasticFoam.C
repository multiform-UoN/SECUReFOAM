/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    poroelasticFoam

Description
    Solves the equations of poroelasticity (coupled Darcy and elasticity).
    A total lagrangian formulation is employed for the stress.

Contributors
    Federico Municchi, Nottingham (2020)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"

    pimpleControl pimple(mesh);

    #include "readPoroelasticProperties.H"
    #include "createFields.H"
    #include "readStressControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while(runTime.run())
    {

        runTime++;
        Info<< "Time: " << runTime.value() << nl << endl;

        while(pimple.loop())
        {
            #include "solveDarcy.H"

            while(pimple.correct())
            {
                #include "solveElasticity.H"
            }

            #include "updateStress.H"

        }


        runTime.write();

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
