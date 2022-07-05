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
    Foam::phaseFieldModels::doubleWell

Description

       TO DO

SourceFiles
    doubleWell.C

Developer
    Juan J. Hidalgo, Matteo Icardi

\*---------------------------------------------------------------------------*/

#ifndef doubleWell_H
#define doubleWell_H

#include "phaseFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFieldModels
{

/*---------------------------------------------------------------------------*\
                  Class doubleWell Declaration
\*---------------------------------------------------------------------------*/

class doubleWell
:
    public phaseFieldModel
{
private:

    //- roots of the double well derivative
    Foam::vector roots_;

    // constant
    scalar k_;

    // epsilon
    scalar eps_;


public:

    //- Runtime type information
    TypeName("doubleWell");

    // Constructors

    //- Construct from components
    doubleWell
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& mu,
        volScalarField& kappa
    );

    //- Destructor
    ~doubleWell()
        {}

    //- Member Functions

        //- Return values
        virtual tmp<scalarField> psi(const scalarField& C);
        virtual tmp<scalarField> dpsidC(const scalarField& C);
        virtual tmp<scalarField> mu(const scalarField& C);
        virtual tmp<scalarField> kappa(const scalarField& C);
        virtual tmp<scalarField> freeEnergy(const scalarField& C);

        //- Write
        virtual void write(Ostream& o) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseFieldModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
