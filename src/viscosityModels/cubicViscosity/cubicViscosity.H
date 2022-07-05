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
    Foam::viscosityModels::cubicViscosity

Description
    Cubic law viscosity model based on

        mu = mu0 + a1*c + a2*c^2 + a3*c^3

SourceFiles
    cubicViscosity.C

Developer
    Juan J. Hidalgo

\*---------------------------------------------------------------------------*/

#ifndef cubicViscosity_H
#define cubicViscosity_H

#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                  Class cubicViscosity Declaration
\*---------------------------------------------------------------------------*/

class cubicViscosity
:
    public viscosityModel
{
private:

    //- model coefficients
  scalar mu0_,a1_,a2_,a3_;


public:

    //- Runtime type information
    TypeName("cubicViscosity");

    // Constructors

    //- Construct from components
    cubicViscosity
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& muField
    );

    //- Destructor
    ~cubicViscosity()
        {}

    //- Member Functions

        //- Return viscosity
        virtual tmp<scalarField> mu(const scalarField& C);

        //- Write
        virtual void write(Ostream& o) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
