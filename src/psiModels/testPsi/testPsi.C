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

#include "testPsi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace psiModels
{
defineTypeNameAndDebug(testPsi, 0);

addToRunTimeSelectionTable
(
    psiModel,
    testPsi,
    dictionary
);
}
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

psiModels::testPsi::testPsi
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
    delta_(readScalar(dict.lookup("delta"))),
    tau_(readScalar(dict.lookup("tau")))
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //


tmp<scalarField> psiModels::testPsi::psi(const scalarField& Sw)
{
    tmp<scalarField> tf
  (
      new scalarField(Sw.size(), 0.)
  );

  scalarField& f = tf.ref();

  f = hcap_*
       (
          Foam::pow(Sw, -1/alpha_)*(1. - exp(beta_*(Sw - ve_)))
          -
          (alpha_/(alpha_ - 1.))*Foam::pow(Sw, (alpha_ - 1.)/alpha_)
          *
          (beta_*exp(beta_*(Sw - ve_)))
       );
  
 return tf;
 }


tmp<scalarField> psiModels::testPsi::dpsidSw(const scalarField& Sw)
{
  tmp<scalarField> tf
  (
      new scalarField(Sw.size(), 0.)
  );

  scalarField& f = tf.ref();

  f = hcap_*
      (
	      (-1./alpha_)*Foam::pow(Sw, -1./alpha_ -1.)*(1. - exp(beta_*(Sw - ve_)))
        -
        2.*Foam::pow(Sw, -1./alpha_)*beta_*exp(beta_*(Sw - ve_))
        -
        (alpha_/(alpha_ -1.))
        *
        Foam::pow(Sw,(alpha_ - 1.)/alpha_)
	      *
        (beta_*beta_*exp(beta_*(Sw - ve_)))
	    );

  return tf;
}


tmp<scalarField> psiModels::testPsi::kappa(const scalarField& Sw)
{
  tmp<scalarField> tf
  (
      new scalarField(Sw.size(), tau_)
  );

  scalarField& f = tf.ref();
  
  return tf;
}

void
psiModels::testPsi::write(Ostream& os) const
{
	psiModel::write(os);
  writeEntry<scalar>(os, "alpha", alpha_);
  writeEntry<scalar>(os, "beta", beta_);
  writeEntry<scalar>(os, "ve", ve_);
  writeEntry<scalar>(os, "hcap", hcap_);
  writeEntry<scalar>(os, "delta", delta_);
  writeEntry<scalar>(os, "tau", tau_);
}


// ************************************************************************* //
