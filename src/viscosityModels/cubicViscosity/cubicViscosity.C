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

#include "cubicViscosity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
defineTypeNameAndDebug(cubicViscosity, 0);

addToRunTimeSelectionTable
(
    viscosityModel,
    cubicViscosity,
    dictionary
);
}
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

viscosityModels::cubicViscosity::cubicViscosity
(
    const dictionary& dict,
    const volScalarField& C,
    volScalarField& muField
)
    :
    viscosityModel(dict, C, muField),
    mu0_(readScalar(dict.lookup("mu0"))),
    a1_(readScalar(dict.lookup("a1"))),
    a2_(readScalar(dict.lookup("a2"))),
    a3_(readScalar(dict.lookup("a3")))
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //

tmp<scalarField> viscosityModels::cubicViscosity::mu(const scalarField& C)
{
  tmp<scalarField> tf
  (
      new scalarField(C.size(),mu0_)
  );

  scalarField& f = tf.ref();

  f += (a1_ + (a2_  + a3_*C)*C)*C;

  return tf;
}


void
viscosityModels::cubicViscosity::write(Ostream& os) const
{
	viscosityModel::write(os);
  writeEntry<scalar>(os, "mu0", mu0_);
  writeEntry<scalar>(os, "a1", a1_);
  writeEntry<scalar>(os, "a2", a2_);
  writeEntry<scalar>(os, "a3", a3_);
}

// ************************************************************************* //
