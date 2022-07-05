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
    Foam::densityModels::linearDensity

Description
    Linear density model based on rho0 and beta coefficients

        rho =rho0 + beta*C

SourceFiles
    linearDensity.C

Developer
    Juan J. Hidalgo, Matteo Icardi

\*---------------------------------------------------------------------------*/

#ifndef linearDensity_H
#define linearDensity_H

#include "densityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace densityModels
{

/*---------------------------------------------------------------------------*\
                  Class linearDensity Declaration
\*---------------------------------------------------------------------------*/

class linearDensity
:
    public densityModel
{
private:

    //- model coefficients
    scalar rho0_,beta_;



public:

    //- Runtime type information
    TypeName("linearDensity");

    // Constructors

    //- Construct from components
    linearDensity
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& rho,
        volScalarField& drhodC
    );

    //- Destructor
    ~linearDensity()
        {}

    //- Member Functions

        //- Return density
        virtual tmp<scalarField> rho(const scalarField& C);

        //- Return density derivative
        virtual tmp<scalarField> drhodC(const scalarField& C);

        //- Write
        virtual void write(Ostream& o) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace densityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
