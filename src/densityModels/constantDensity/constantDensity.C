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

#include "constantDensity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace densityModels
{
defineTypeNameAndDebug(constantDensity, 0);

addToRunTimeSelectionTable
(
    densityModel,
    constantDensity,
    dictionary
);
}
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

densityModels::constantDensity::constantDensity
(
    const dictionary& dict,
    const volScalarField& C,
    volScalarField& rho,
    volScalarField& drhodC
)
    :
    densityModel(dict, C, rho, drhodC),
    rho0_(readScalar(dict.lookup("rho0")))
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //


tmp<scalarField> densityModels::constantDensity::rho(const scalarField& C)
{
  tmp<scalarField> tf
  (
      new scalarField(C.size(),rho0_)
  );

  return tf;
}

tmp<scalarField> densityModels::constantDensity::drhodC(const scalarField& C)
{
  tmp<scalarField> tf
  (
      new scalarField(C.size(),scalar(0))
  );

  return tf;
}



void
densityModels::constantDensity::write(Ostream& os) const
{
	densityModel::write(os);
  writeEntry<scalar>(os, "rho0", rho0_);
}





// ************************************************************************* //


// ************************************************************************* //
