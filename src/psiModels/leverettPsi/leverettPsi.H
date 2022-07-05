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

Class
    Foam::psiModels::leverettPsi

Description

       TO DO

ourceFiles
    leverettPsi.C

Developer
    Juan J. Hidalgo, Matteo Icardi

\*---------------------------------------------------------------------------*/

#ifndef leverettPsi_H
#define leverettPsi_H

#include "psiModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace psiModels
{

/*---------------------------------------------------------------------------*\
                  Class leverettPsi Declaration
\*---------------------------------------------------------------------------*/

class leverettPsi
:
    public psiModel
{
private:

    //- model coefficients
  scalar alpha_, beta_, ve_, hcap_, delta_;



public:

    //- Runtime type information
    TypeName("leverettPsi");

    // Constructors

    //- Construct from components
    leverettPsi
    (
        const dictionary& dict,
        const volScalarField& Sw,
        volScalarField& psi,
        volScalarField& dpsidSw,
        volScalarField& kappa
    );

    //- Destructor
    ~leverettPsi()
        {}

    //- Member Functions

        //- Return psi, dpsidSw, kappa
        virtual tmp<scalarField> psi(const scalarField& Sw);
        virtual tmp<scalarField> dpsidSw(const scalarField& Sw);
        virtual tmp<scalarField> kappa(const scalarField& Sw);


    //- Write
    virtual void write(Ostream& o) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace psiModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
