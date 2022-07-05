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

#include "viscosityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace viscosityModels
  {
    defineTypeNameAndDebug(viscosityModel, 0);
    defineRunTimeSelectionTable(viscosityModel, dictionary);
  }
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

viscosityModels::viscosityModel::viscosityModel
(
    const dictionary& dict,
    const volScalarField& C,
    volScalarField& mu
)
    :
    dict_(dict),
    C_(C),
    mu_(mu)
{}

/*------------------------  Member functions----------------------------------*/

tmp<scalarField> viscosityModels::viscosityModel::mu()
{
  return mu(C_);
}

tmp<scalarField> viscosityModels::viscosityModel::mu(const volScalarField& C)
{
  return mu(C.field());
}


void viscosityModels::viscosityModel::correctMu()
{
  mu_.field() = mu();
  forAll (mu_.boundaryField(), patchI)
    {
      mu_.boundaryFieldRef()[patchI] = mu(C_.boundaryField()[patchI]);
    }
}


void
viscosityModels::viscosityModel::write(Ostream& os) const
{
	dict_.write(os);
}

// ************************************************************************* //
