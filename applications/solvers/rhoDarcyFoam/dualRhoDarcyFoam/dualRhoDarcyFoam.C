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
dualRhoDarcyFoam

Description
Solves the Darcy equation with variable density for matrix and fractures

Developers
Federico Municchi
Matteo Icardi

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "densityModel.H"
#include "viscosityModel.H"
#include "reverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

  #include "setRootCase.H"
  #include "createTime.H"
  #include "createTimeControls.H"
  #include "createMesh.H"

  pimpleControl pimple(mesh);

  #include "readGravitationalAcceleration.H"
  #include "createFields.H"
  #include "initContinuityErrs.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  #include "CourantNo.H"
  #include "CourantNoG_fr.H"
  #include "setInitialDeltaT.H"

  while(runTime.loop())
  {
    #include "readTimeControls.H"

    Info << "\nTime = " << runTime.timeName() << nl << endl;

    while (pimple.loop())
    {
      //- Concentration equation
      {
        #include "cEqn.H"
      }

      //- Update density/viscosity and, if needed, PhiG and Mf
      #include "updatePhiG.H"
      #include "updateMf.H"

      //- Pressure equation
      {
        if (hydrostaticPressure)
        {
          volScalarField&  pp(p_rgh);
          volScalarField&  pp_fr(p_rgh_fr);
          #include "pEqn.H"
          p = p_rgh + rho*gh;
          p_fr = p_rgh_fr + rho_fr*gh;
          p.correctBoundaryConditions();
          p_fr.correctBoundaryConditions();
        }
        else
        {
          volScalarField&   pp(p);
          volScalarField&   pp_fr(p_fr);
          #include "pEqn.H"
          p_rgh = p - rho*gh;
          p_rgh_fr = p_fr - rho_fr*gh;
        }
      }

      // - Reconstruct velocity fields
      U = fvc::reconstruct(phi);
      U_fr = fvc::reconstruct(phi_fr);
      #include "updateD.H"
    }

    #include "continuityErrs.H"
    #include "CourantNo.H"
    #include "CourantNoG_fr.H"
    #include "setDeltaT.H"
    #include "writeInfo.H"

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;

  }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
