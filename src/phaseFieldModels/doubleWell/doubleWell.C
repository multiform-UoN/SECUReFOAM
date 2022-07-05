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

#include "doubleWell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFieldModels
{
defineTypeNameAndDebug(doubleWell, 0);

addToRunTimeSelectionTable
(
    phaseFieldModel,
    doubleWell,
    dictionary
);
}
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

phaseFieldModels::doubleWell::doubleWell
(
    const dictionary& dict,
    const volScalarField& C,
    volScalarField& mu,
    volScalarField& kappa
)
    :
    phaseFieldModel(dict, C, mu, kappa),
    roots_(dict.lookupOrDefault<vector>("roots",vector(0.0,0.5,1.0))),
    k_(readScalar(dict.lookup("k"))),
    eps_(readScalar(dict.lookup("eps")))
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //

tmp<scalarField> phaseFieldModels::doubleWell::psi(const scalarField& C)
{
  tmp<scalarField> tf
  (
      new scalarField(C.size(), 0.)
  );

  scalarField& f = tf.ref();

  f = k_*(C-roots_[0])*(C-roots_[1])*(C-roots_[2]);

  return tf;
}


tmp<scalarField> phaseFieldModels::doubleWell::mu(const scalarField& C)
{
  tmp<scalarField> tf
  (
      new scalarField(C.size(), 0.)
  );

  scalarField& f = tf.ref();

  f = k_*
      (
        (C-roots_[0])*(C-roots_[1])
        +
        (C-roots_[0])*(C-roots_[2])
        +
        (C-roots_[1])*(C-roots_[2])
      );

  return tf;
}


tmp<scalarField> phaseFieldModels::doubleWell::dpsidC(const scalarField& C)
{
  tmp<scalarField> tf
  (
      new scalarField(C.size(), 0.)
  );

  scalarField& f = tf.ref();

  f = k_*
      (
        (C-roots_[0])*(C-roots_[1])
        +
        (C-roots_[0])*(C-roots_[2])
        +
        (C-roots_[1])*(C-roots_[2])
      );

  return tf;
}


tmp<scalarField> phaseFieldModels::doubleWell::kappa(const scalarField& C)
{
  return kappa_.field();  // does nothing, return initialised value
}



void
phaseFieldModels::doubleWell::write(Ostream& os) const
{
	phaseFieldModel::write(os);
  writeEntry<vector>(os, "roots", roots_);
  writeEntry<scalar>(os, "k", k_);
  writeEntry<scalar>(os, "eps", eps_);
}

// ************************************************************************* //
