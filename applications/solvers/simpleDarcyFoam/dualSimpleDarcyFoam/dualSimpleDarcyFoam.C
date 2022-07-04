/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | foam-extend: Open Source CFD
\\    /   O peration     | Version:     4.0
\\  /    A nd           | Web:         http://www.foam-extend.org
\\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
This file is part of foam-extend.

foam-extend is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

foam-extend is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
dualSimpleDarcyFoam

Description
Solves the Darcy equation for matrix and fractures

Developers
Federico Municchi
Matteo Icardi

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "reverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

  simpleControl simple(mesh);

  #include "readGravitationalAcceleration.H"
  #include "createFields.H"


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nStarting time loop\n" << endl;

  while (simple.loop(runTime))
  {

    Info << "\nTime = " << runTime.timeName() << nl << endl;

    #include "pEqn.H"

    // - Reconstruct velocity fields
    U = fvc::reconstruct(phi);
    U_fr = fvc::reconstruct(phi_fr);

    runTime.write();

    Info << "Mean velocity  = "
        << gSum(U.primitiveField()*mesh.V())/gSum(mesh.V())
        << nl << endl;

    Info << "Mean fracture velocity  = "
        << gSum(U_fr.primitiveField()*mesh.V())/gSum(mesh.V())
        << nl << endl;

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;


  }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
