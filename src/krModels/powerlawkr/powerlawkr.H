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
    Foam::krModels::powerlawkr

Description

        kr =Se^a
        Se = (Sw - Srw)/(1.- Srw)

SourceFiles
    powerlawkr.C

Developer
    Juan J. Hidalgo, Matteo Icardi

\*---------------------------------------------------------------------------*/

#ifndef powerlawKr_H
#define powerlawKr_H

#include "krModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace krModels
{

/*---------------------------------------------------------------------------*\
                  Class powelawKr Declaration
\*---------------------------------------------------------------------------*/

class powerlawKr
:
    public krModel
{
private:

    //- model coefficients
    scalar a_, Srw_;



public:

    //- Runtime type information
    TypeName("powerlawKr");

    // Constructors

    //- Construct from components
    powerlawKr
    (
        const dictionary& dict,
        const volScalarField& Sw,
        volScalarField& kr
    );

    //- Destructor
    ~powerlawKr()
        {}

    //- Member Functions

    //- Return kr
    tmp<scalarField> kr(const scalarField& Sw);

    void write(Ostream& os) const;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace krModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
