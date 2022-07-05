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
    Foam::viscosityModels::constantViscosity

Description
    Constant viscosity model based on mu0.

        mu = mu0

SourceFiles
    constantViscosity.C

Developer
    Juan J. Hidalgo

\*---------------------------------------------------------------------------*/

#ifndef constantViscosity_H
#define constantViscosity_H

#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                  Class constantViscosity Declaration
\*---------------------------------------------------------------------------*/

class constantViscosity
:
    public viscosityModel
{
private:

    //- model coefficients
    scalar mu0_;

        
public:

    //- Runtime type information
    TypeName("constantViscosity");

    // Constructors

    //- Construct from components
    constantViscosity
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& muField
    );

    //- Destructor
    ~constantViscosity()
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
