/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | 
  \\    /   O peration     | 
   \\  /    A nd           | 
    \\/     M anipulation  | 
-------------------------------------------------------------------------------
Application
    setRandomField

Description
    Create Gaussian Random Fields using a spectral method
    Ref: Heße et al, Env Model Soft 55 (2014)

Developers
    Matteo Icardi

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"

using namespace Foam;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  #include "setRootCase.H"

  #include "createTime.H"
  #include "createMesh.H"

  #include "setupGRF.H"

  Info << endl << "Reading field " << varName << endl;

  if (fieldHeader.typeHeaderOk<volScalarField>(true))
  {
    Info << "  Field is a scalar" << endl;
    volScalarField f(fieldHeader, mesh);
    const scalar f0(1.0);
    #include "createVolFields.H"
    #include "computeRandomField.H"
  }
  else if (fieldHeader.typeHeaderOk<volVectorField>(true))
  {
    Info << "  Field is a vector" << endl;
    volVectorField f(fieldHeader, mesh);
    const vector f0(vector::one);
    #include "createVolFields.H"
    #include "computeRandomField.H"
  }
  else if (fieldHeader.typeHeaderOk<volTensorField>(true))
  {
    Info << "  Field is a tensor" << endl;
    volTensorField f(fieldHeader, mesh);
    const tensor f0(tensor::I);
    #include "createVolFields.H"
    #include "computeRandomField.H"
  }
  else if (fieldHeader.typeHeaderOk<volSymmTensorField>(true))
  {
    Info << "  Field is a symmetric tensor" << endl;
    volSymmTensorField f(fieldHeader, mesh);
    const symmTensor f0(symmTensor::I);
    #include "createVolFields.H"
    #include "computeRandomField.H"
  }
  else if (fieldHeader.typeHeaderOk<surfaceScalarField>(true))
  {
    Info << "  Field is a surface scalar" << endl;
    surfaceScalarField f(fieldHeader, mesh);
    const scalar f0(1.0);
    #include "createSurfaceFields.H"
    #include "computeRandomField.H"
  }
  else if (fieldHeader.typeHeaderOk<surfaceVectorField>(true))
  {
    Info << "  Field is a surface vector" << endl;
    surfaceVectorField f(fieldHeader, mesh);
    const vector f0(vector::one);
    #include "createSurfaceFields.H"
    #include "computeRandomField.H"
  }
  else if (fieldHeader.typeHeaderOk<surfaceTensorField>(true))
  {
    Info << "  Field is a surface tensor" << endl;
    surfaceTensorField f(fieldHeader, mesh);
    const tensor f0(tensor::I);
    #include "createSurfaceFields.H"
    #include "computeRandomField.H"
  }
  else if (fieldHeader.typeHeaderOk<surfaceSymmTensorField>(true))
  {
    Info << "  Field is a surface symmetric tensor" << endl;
    surfaceSymmTensorField f(fieldHeader, mesh);
    const symmTensor f0(symmTensor::I);
    #include "createSurfaceFields.H"
    #include "computeRandomField.H"
  }
  else
  {
    FatalError
    << "  There is no field " << varName << endl
    << exit(FatalError);
  }

  Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s"
      << nl << endl;


  Info<< "End\n" << endl;

  return 0;

}


// ************************************************************************* //
