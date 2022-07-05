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

#include "constantkr.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace krModels
{
defineTypeNameAndDebug(constantKr, 0);

addToRunTimeSelectionTable
(
    krModel,
    constantKr,
    dictionary
);
}
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

krModels::constantKr::constantKr
(
    const dictionary& dict,
    const volScalarField& Sw,
    volScalarField& kr
)
    :
    krModel(dict, Sw, kr),
    kr0_(readScalar(dict.lookup("kr0")))
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //


tmp<scalarField> krModels::constantKr::kr(const scalarField& Sw)
{
  tmp<scalarField> tf
  (
      new scalarField(Sw.size(),kr0_)
  );

  return tf;
}


void
krModels::constantKr::write(Ostream& os) const
{
	krModel::write(os);
  writeEntry<scalar>(os, "kr0", kr0_);
}



// ************************************************************************* //


// ************************************************************************* //
