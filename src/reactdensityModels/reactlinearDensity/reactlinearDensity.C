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

#include "reactlinearDensity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactdensityModels
{
defineTypeNameAndDebug(reactlinearDensity, 0);

addToRunTimeSelectionTable
(
    reactdensityModel,
    reactlinearDensity,
    dictionary
);
}
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactdensityModels::reactlinearDensity::reactlinearDensity
(
    const dictionary& dict,
    const volScalarField& cA,
    const volScalarField& cB,
    const volScalarField& cC,
    volScalarField& rho,
    volScalarField& drhodC
)
    :
  reactdensityModel(dict, cA, cB, cC, rho, drhodC),
    rho0_(readScalar(dict.lookup("rho0"))),
    Ra_(readScalar(dict.lookup("Ra"))),
    Rb_(readScalar(dict.lookup("Rb"))),
    Rc_(readScalar(dict.lookup("Rc")))
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //


tmp<scalarField> reactdensityModels::reactlinearDensity::rho(const scalarField& cA, const scalarField& cB, const scalarField& cC)
{
  tmp<scalarField> tf
  (
      new scalarField(cA.size(),rho0_)
  );

  scalarField& f = tf.ref();

  f +=  Ra_*cA + Rb_*cB + Rc_*cC;

  return tf;
}


tmp<scalarField> reactdensityModels::reactlinearDensity::drhodC(const scalarField& cA)
{
  tmp<scalarField> tf
  (
      new scalarField(cA.size(),Ra_)
  );

  return tf;
}


void
reactdensityModels::reactlinearDensity::write(Ostream& os) const
{
	reactdensityModel::write(os);
  writeEntry<scalar>(os, "rho0", rho0_);
  writeEntry<scalar>(os, "Ra", Ra_);
  writeEntry<scalar>(os, "Rb", Rb_);
  writeEntry<scalar>(os, "Rc", Rc_);
}



// ************************************************************************* //
