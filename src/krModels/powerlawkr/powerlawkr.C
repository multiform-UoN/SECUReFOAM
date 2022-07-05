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

#include "powerlawkr.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace krModels
{
defineTypeNameAndDebug(powerlawKr, 0);

addToRunTimeSelectionTable
(
    krModel,
    powerlawKr,
    dictionary
);
}
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

krModels::powerlawKr::powerlawKr
(
    const dictionary& dict,
    const volScalarField& Sw,
    volScalarField& kr
)
    :
    krModel(dict, Sw, kr),
    a_(readScalar(dict.lookup("a"))),
    Srw_(readScalar(dict.lookup("Srw")))
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //


tmp<scalarField> krModels::powerlawKr::kr(const scalarField& Sw)
{
  tmp<scalarField> tf
  (
      new scalarField(Sw.size(),Srw_)
  );

  scalarField& f = tf.ref();

  f = Foam::pow((Sw - Srw_)/(1. - Srw_), a_);

  return tf;
}


void
krModels::powerlawKr::write(Ostream& os) const
{
	krModel::write(os);
  writeEntry<scalar>(os, "Srw", Srw_);
  writeEntry<scalar>(os, "a", a_);
}


// ************************************************************************* //
