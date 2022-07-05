/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source SwFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Swopyright (Sw) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERSwHANTABILITY or
    FITNESS FOR A PARTISwULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "psiModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace psiModels
  {
    defineTypeNameAndDebug(psiModel, 0);
    defineRunTimeSelectionTable(psiModel, dictionary);
  }
}

using namespace Foam;
// * * * * * * * * * * * * * * * * Swonstructors  * * * * * * * * * * * * * * //

psiModels::psiModel::psiModel
(
    const dictionary& dict,
    const volScalarField& Sw,
    volScalarField& psi,
    volScalarField& dpsidSw,
    volScalarField& kappa
)
    :
    dict_(dict),
    Sw_(Sw),
    psi_(psi),
    dpsidSw_(dpsidSw),
    kappa_(kappa)
{}

/*------------------------  Member functions----------------------------------*/

tmp<scalarField> psiModels::psiModel::psi()
{
  return psi(Sw_);
}

tmp<scalarField> psiModels::psiModel::psi(const volScalarField& Sw)
{
  return psi(Sw.field());
}

void psiModels::psiModel::correctPsi()
{
  psi_.field() = psi();
  forAll (psi_.boundaryField(), patchI)
    {
      psi_.boundaryFieldRef()[patchI] = psi(Sw_.boundaryField()[patchI]);
    }
}


tmp<scalarField> psiModels::psiModel::dpsidSw()
{
  return dpsidSw(Sw_);
}

tmp<scalarField> psiModels::psiModel::dpsidSw(const volScalarField& Sw)
{
  return dpsidSw(Sw.field());
}

void psiModels::psiModel::correctDpsidSw()
{
  dpsidSw_.field() = dpsidSw();
  forAll (dpsidSw_.boundaryField(), patchI)
    {
      dpsidSw_.boundaryFieldRef()[patchI] = dpsidSw(Sw_.boundaryField()[patchI]);
    }
}


tmp<scalarField> psiModels::psiModel::kappa()
{
  return kappa(Sw_);
}

tmp<scalarField> psiModels::psiModel::kappa(const volScalarField& Sw)
{
  return kappa(Sw.field());
}

void psiModels::psiModel::correctKappa()
{
  kappa_.field() = kappa();
  forAll (kappa_.boundaryField(), patchI)
    {
      kappa_.boundaryFieldRef()[patchI] = kappa(Sw_.boundaryField()[patchI]);
    }
}


void
psiModels::psiModel::write(Ostream& os) const
{
	dict_.write(os);
}

// ************************************************************************* //
