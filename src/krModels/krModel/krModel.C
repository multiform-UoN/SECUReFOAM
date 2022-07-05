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

#include "krModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace krModels
  {
    defineTypeNameAndDebug(krModel, 0);
    defineRunTimeSelectionTable(krModel, dictionary);
  }
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

krModels::krModel::krModel
(
    const dictionary& dict,
    const volScalarField& Sw,
    volScalarField& kr
)
    :
    dict_(dict),
    Sw_(Sw),
    kr_(kr)
{}

/*------------------------  Member functions----------------------------------*/

tmp<scalarField> krModels::krModel::kr()
{
  return kr(Sw_);
}

tmp<scalarField> krModels::krModel::kr(const volScalarField& Sw)
{
  return kr(Sw.field());
}


void krModels::krModel::correctKr()
{
  kr_.field() = kr();
  forAll (kr_.boundaryField(), patchI)
    {
      kr_.boundaryFieldRef()[patchI] = kr(Sw_.boundaryField()[patchI]);
    }
}

void
krModels::krModel::write(Ostream& os) const
{
	dict_.write(os);
}

// ************************************************************************* //
