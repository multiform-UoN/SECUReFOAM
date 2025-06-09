/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | 
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
Application
    rhoDarcyFoam

Description
    Solves the Darcy equation with variable density/viscosity

Developers
    Federico Municchi
    Matteo Icardi
    Juan Hidalgo
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
      Info<< "A0\n" << endl;

      //- Update density/viscosity and, if needed, PhiG and Mf
      #include "updatePhiG.H"
      Info<< "A1\n" << endl;
      #include "updateMf.H"
      Info<< "A2\n" << endl;

      //- Pressure equation
      if (hydrostaticPressure)
      {
        volScalarField&  pp(p_rgh);
        #include "pEqn.H"
        Info<< "A33\n" << endl;
        p = p_rgh + rho*gh;
      }
      else
      {
        volScalarField&   pp(p);
        #include "pEqn.H"
        Info<< "A3\n" << endl;
        p_rgh = p - rho*gh;
      }
    
      // - Reconstruct velocity fields and recompute dispersion
      U = fvc::reconstruct(phi);
      Info<< "A4\n" << endl;
      #include "updateD.H"
      Info<< "A5\n" << endl;
    }

    Info<< "A\n" << endl;
    
    #include "postProcessing.H"
    Info<< "B\n" << endl;
    #include "continuityErrs.H"
    Info<< "C\n" << endl;
    #include "CourantNo.H"
    Info<< "D\n" << endl;
    #include "setDeltaT.H"
    Info<< "E\n" << endl;

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
