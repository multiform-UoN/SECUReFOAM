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

\*---------------------------------------------------------------------------*/

#include "densityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace densityModels
  {
    defineTypeNameAndDebug(densityModel, 0);
    defineRunTimeSelectionTable(densityModel, dictionary);
  }
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

densityModels::densityModel::densityModel
(
    const dictionary& dict,
    const volScalarField& C,
    volScalarField& rho,
    volScalarField& drhodC
)
    :
    dict_(dict),
    C_(C),
    rho_(rho),
    drhodC_(drhodC)
{}

/*------------------------  Member functions----------------------------------*/


tmp<scalarField> densityModels::densityModel::rho()
{
  return rho(C_);
}

tmp<scalarField> densityModels::densityModel::rho(const volScalarField& C)
{
  return rho(C.field());
}

tmp<scalarField> densityModels::densityModel::drhodC()
{
  return drhodC(C_);
}

tmp<scalarField> densityModels::densityModel::drhodC(const volScalarField& C)
{
  return drhodC(C.field());
}


void densityModels::densityModel::correctRho()
{
  rho_.field() = rho();
  forAll (rho_.boundaryField(), patchI)
    {
      rho_.boundaryFieldRef()[patchI] = rho(C_.boundaryField()[patchI]);
    }
}


void densityModels::densityModel::correctDrhodC()
{
  drhodC_.field() = drhodC();
  forAll (drhodC_.boundaryField(), patchI)
    {
      drhodC_.boundaryFieldRef()[patchI] = drhodC(C_.boundaryField()[patchI]);
    } 
}


void
densityModels::densityModel::write(Ostream& os) const
{
	dict_.write(os);
}

// ************************************************************************* //
