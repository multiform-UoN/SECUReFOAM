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

#include "leverettPsi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace psiModels
{
defineTypeNameAndDebug(leverettPsi, 0);

addToRunTimeSelectionTable
(
    psiModel,
    leverettPsi,
    dictionary
);
}
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

psiModels::leverettPsi::leverettPsi
(
    const dictionary& dict,
    const volScalarField& Sw,
    volScalarField& psi,
    volScalarField& dpsidSw,
    volScalarField& kappa
)
    :
    psiModel(dict, Sw, psi, dpsidSw, kappa),
    alpha_(readScalar(dict.lookup("alpha"))),
    beta_(readScalar(dict.lookup("beta"))),
    ve_(readScalar(dict.lookup("ve"))),
    hcap_(readScalar(dict.lookup("hcap"))),
    delta_(readScalar(dict.lookup("delta")))
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //

tmp<scalarField> psiModels::leverettPsi::psi(const scalarField& Sw)
{
  tmp<scalarField> tf
  (
      new scalarField(Sw.size(), 0.)
  );

  scalarField& f = tf.ref();

  f = (hcap_/Foam::pow(Sw, 1./alpha_))*
      (1. - Foam::exp(beta_*(Sw - ve_)))*
      (1. + beta_*(alpha_/(alpha_ - 1.))*Sw);

 return tf;
}


tmp<scalarField> psiModels::leverettPsi::dpsidSw(const scalarField& Sw)
{
  tmp<scalarField> tf
  (
      new scalarField(Sw.size(), 0.)
  );

  scalarField& f(tf.ref());

  f =
    psi(Sw) * (-1./(alpha_*Sw))
    -
    (hcap_/Foam::pow(Sw, 1/alpha_))*(beta_*Foam::exp(beta_*(Sw - ve_)))
    *
    (1. + beta_*(alpha_/(alpha_ - 1.))*Sw)
    +
    (hcap_/Foam::pow(Sw, 1/alpha_))*(1. - Foam::exp(beta_*(Sw - ve_)))
    *
    (beta_*(alpha_/(alpha_ - 1.)));

  return tf;
}


tmp<scalarField> psiModels::leverettPsi::kappa(const scalarField& Sw)
{
  tmp<scalarField> tf
  (
      new scalarField(Sw.size(), 0.)
  );

    scalarField& f = tf.ref();

    f = hcap_*
        Foam::pow(delta_, 2)*(alpha_/(alpha_ - 1.))
        *
        Foam::pow(Sw, (alpha_ - 1.)/alpha_)
        *
        (1. - Foam::exp(beta_*(Sw - ve_)));

  
  return tf;
}


void
psiModels::leverettPsi::write(Ostream& os) const
{
	psiModel::write(os);
  writeEntry<scalar>(os, "alpha", alpha_);
  writeEntry<scalar>(os, "beta", beta_);
  writeEntry<scalar>(os, "ve", ve_);
  writeEntry<scalar>(os, "hcap", hcap_);
  writeEntry<scalar>(os, "delta", delta_);
}


// ************************************************************************* //
