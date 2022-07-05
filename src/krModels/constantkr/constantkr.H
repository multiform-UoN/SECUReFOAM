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
    Foam::krModels::constantKr

Description
    Constant kr model based on kr0.

        kr = kr0

SourceFiles
    constantKr.C

Developer
    Juan J. Hidalgo, Matteo Icardi

\*---------------------------------------------------------------------------*/

#ifndef constantKr_H
#define constantKr_H

#include "krModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace krModels
{

/*---------------------------------------------------------------------------*\
                  Class constantKr Declaration
\*---------------------------------------------------------------------------*/

class constantKr
:
    public krModel
{
private:

    //- model coefficients
    scalar kr0_;


public:

    //- Runtime type information
    TypeName("constantKr");

    // Constructors

    //- Construct from components
    constantKr
    (
        const dictionary& dict,
        const volScalarField& Sw,
        volScalarField& kr
    );

    //- Destructor
    ~constantKr()
        {}

    //- Member Functions

    //- True unless trivial constant law
    bool needsUpdate()
        {
            return 0;
        }
    
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
