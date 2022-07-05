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
    Foam::viscosityModels::tabulatedViscosity

Description
    Viscosity model based on tabulated data.

        mu = mu(c)

Based on thermophysicalModels/specie/transport/tabulated/tabulatedTransport.

SourceFiles
    tabulatedViscosity.C

Developer
    Juan J. Hidalgo

\*---------------------------------------------------------------------------*/

#ifndef tabulatedViscosity_H
#define tabulatedViscosity_H

#include "viscosityModel.H"
#include "nonUniformTableThermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                  Class tabulatedViscosity Declaration
\*---------------------------------------------------------------------------*/

class tabulatedViscosity
:
    public viscosityModel
{
  typedef thermophysicalFunctions::nonUniformTable nonUniformTable;

    // Private Data
private:

        //- Dynamic viscosity table
        nonUniformTable mu_;


public:

    //- Runtime type information
    TypeName("tabulatedViscosity");

    // Constructors

    //- Construct from components
    tabulatedViscosity
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& muField
    );

    //- Destructor
    ~tabulatedViscosity()
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
