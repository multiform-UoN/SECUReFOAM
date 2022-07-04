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

      //- Update density/viscosity and, if needed, PhiG and Mf
      #include "updatePhiG.H"
      #include "updateMf.H"

      //- Pressure equation
      if (hydrostaticPressure)
      {
        volScalarField&  pp(p_rgh);
        #include "pEqn.H"
        p = p_rgh + rho*gh;
      }
      else
      {
        volScalarField&   pp(p);
        #include "pEqn.H"
        p_rgh = p - rho*gh;
      }
    
      // - Reconstruct velocity fields and recompute dispersion
      U = fvc::reconstruct(phi);
      #include "updateD.H"
    }
    
    #include "postProcessing.H"
    #include "continuityErrs.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"

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
