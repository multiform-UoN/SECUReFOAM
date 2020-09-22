/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | foam-extend: Open Source CFD
\\    /   O peration     | Version:     4.0
\\  /    A nd           | Web:         http://www.foam-extend.org
\\/     M anipulation  | Matteo Icardi, 2019
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
setRandomField

Description
Sets a random field
Ref: Heße et al, Env Model Soft 55 (2014)

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

  Info<< "  Reading field " << varName << endl << endl;

  if (fieldHeader.typeHeaderOk<volScalarField>(true))
  {
    Info << "  Field is a scalar" << endl;
    volScalarField f(fieldHeader, mesh);
    const scalar f0(1.0);
    #include "computeRandomField.H"
  }
  else if (fieldHeader.typeHeaderOk<volVectorField>(true))
  {
    Info << "  Field is a vector" << endl;
    volVectorField f(fieldHeader, mesh);
    const vector f0(vector::one);
    #include "computeRandomField.H"
  }
  else if (fieldHeader.typeHeaderOk<volTensorField>(true))
  {
    Info << "  Field is a tensor" << endl;
    volTensorField f(fieldHeader, mesh);
    const tensor f0(tensor::I);
    #include "computeRandomField.H"
  }
  else
  {
    FatalError
    << "  There is no field " << varName << endl
    << exit(FatalError);
  }

  Info<< "End\n" << endl;

  return 0;
}


// ************************************************************************* //
