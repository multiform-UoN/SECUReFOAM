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

#include "tabulatedViscosity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
defineTypeNameAndDebug(tabulatedViscosity, 0);

addToRunTimeSelectionTable
(
    viscosityModel,
    tabulatedViscosity,
    dictionary
);
}
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

viscosityModels::tabulatedViscosity::tabulatedViscosity
(
    const dictionary& dict,
    const volScalarField& C,
    volScalarField& muField
)
    :
    viscosityModel(dict, C, muField),
    mu_("mu", dict.subDict("values"))
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //


tmp<scalarField> viscosityModels::tabulatedViscosity::mu(const scalarField& C)
{
  tmp<scalarField> tf
  (
      new scalarField(C.size(), 0.0)
  );

  scalarField& f = tf.ref();
  
  forAll (C, i)
  {
    f[i] = mu_.f(0, C[i]);
  }
  return tf;
}


void
viscosityModels::tabulatedViscosity::write(Ostream& os) const
{
	viscosityModel::write(os);
  dictionary dict("mu");
  dict.add("values", mu_.values());
  os  << indent << dict.dictName() << dict;
}


// ************************************************************************* //
