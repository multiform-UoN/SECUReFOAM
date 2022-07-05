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

#include "phaseFieldModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace phaseFieldModels
  {
    defineTypeNameAndDebug(phaseFieldModel, 0);
    defineRunTimeSelectionTable(phaseFieldModel, dictionary);
  }
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

phaseFieldModels::phaseFieldModel::phaseFieldModel
(
    const dictionary& dict,
    const volScalarField& C,
    volScalarField& mu,
    volScalarField& kappa
)
    :
    dict_(dict),
    C_(C),
    mu_(mu),
    kappa_(kappa)
{}

/*------------------------  Member functions----------------------------------*/

tmp<scalarField> phaseFieldModels::phaseFieldModel::psi()
{
  return psi(C_);
}

tmp<scalarField> phaseFieldModels::phaseFieldModel::psi(const volScalarField& C)
{
  return psi(C.field());
}

tmp<scalarField> phaseFieldModels::phaseFieldModel::dpsidC()
{
  return dpsidC(C_);
}

tmp<scalarField> phaseFieldModels::phaseFieldModel::dpsidC(const volScalarField& C)
{
  return dpsidC(C.field());
}

tmp<scalarField> phaseFieldModels::phaseFieldModel::mu()
{
  return mu(C_);
}

tmp<scalarField> phaseFieldModels::phaseFieldModel::mu(const volScalarField& C)
{
  return mu(C.field());
}

tmp<scalarField> phaseFieldModels::phaseFieldModel::kappa()
{
  return kappa(C_);
}

tmp<scalarField> phaseFieldModels::phaseFieldModel::kappa(const volScalarField& C)
{
  return kappa(C.field());
}

tmp<scalarField> phaseFieldModels::phaseFieldModel::freeEnergy()
{
  return freeEnergy(C_);
}

tmp<scalarField> phaseFieldModels::phaseFieldModel::freeEnergy(const volScalarField& C)
{
  return freeEnergy(C.field());
}


void phaseFieldModels::phaseFieldModel::correctMu()
{
  mu_.field() = mu();
  forAll (mu_.boundaryField(), patchI)
    {
      mu_.boundaryFieldRef()[patchI] = mu(C_.boundaryField()[patchI]);
    }
}


void phaseFieldModels::phaseFieldModel::correctKappa()
{
  kappa_.field() = kappa();
  forAll (kappa_.boundaryField(), patchI)
    {
      kappa_.boundaryFieldRef()[patchI] = kappa(C_.boundaryField()[patchI]);
    } 
}


void
phaseFieldModels::phaseFieldModel::write(Ostream& os) const
{
	dict_.write(os);
}

// ************************************************************************* //
