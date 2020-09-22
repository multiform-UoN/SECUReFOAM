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
    rhoDarcyFoam

Description
    Solves the Darcy equation with variable density

Developers
    Federico Municchi
    Matteo Icardi

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

  pimpleControl pimple(mesh);

  #include "readGravitationalAcceleration.H"
  #include "createFields.H"
  #include "initContinuityErrs.H"




  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



  while(runTime.loop())
  {

    Info << "\nTime = " << runTime.timeName() << nl << endl;

    while (pimple.loop())
    {
      #include "updateRho.H"
      #include "updatePhiG.H"

      {
        #include "pEqn.H"
      }

      // - Reconstruct velocity fields
      U = fvc::reconstruct(phi + phiG);

      {
        #include "cEqn.H"
      }
    }

    #include "continuityErrs.H"

    Info << "Mean velocity  = "
        << gSum(U.primitiveField()*mesh.V())/gSum(mesh.V())
        << nl << endl;

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;

  }

  Info<< "End\n" << endl;

  return 0;
}

// ************************************************************************* //
