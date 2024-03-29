/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    simpleDarcyFoam

Description
    Steady-state solver for incompressible flow in porous medium with
    tensor formulation for the permeability.

Developers
    Federico Municchi, University of Nottingham (2019)

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
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {

            fvScalarMatrix pEqn
            (
                -  fvm::laplacian(Mf,p)
                +  fvc::div(phiG)
                ==
                dims*fvOptions(p)

            );

            fvOptions.constrain(pEqn);
            pEqn.relax();
            pEqn.solve();
            fvOptions.correct(p);

            if (simple.finalNonOrthogonalIter())
            {
              phi = pEqn.flux() + phiG;
              U = fvc::reconstruct(phi);
              // Info << "Mass conservation "
              //     << fvc::domainIntegrate(fvc::div(phi))/gSum(mesh.V())
              //     << endl;
              // Info << "Mean vel from eq flux "
              //     << gSum(U.primitiveField()*mesh.V())/gSum(mesh.V())
              //     << nl << endl;
              // U = -(K & (fvc::grad(p) - rho*g) / mu);
              // phi = fvc::flux(U);
              // Info << "Mass conservation "
              //     << fvc::domainIntegrate(fvc::div(U))/gSum(mesh.V())
              //     << endl;
              Info << "Mean velocity  = "
                  << gSum(U.primitiveField()*mesh.V())/gSum(mesh.V())
                  << nl << endl;
            }

        }

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
