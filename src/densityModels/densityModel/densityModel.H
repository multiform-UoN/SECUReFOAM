/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Juan Hidalgo, Matteo Icardi
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

Namespace
    Foam::densityModel

Description
    A namespace for various density implementations.

Class
    Foam::densityModel

Description
    An abstract base class for density models.

SourceFiles
    densityModel.C
    densityModelNew.C

Acknowledgements

  Based on the porusModels from https://github.com/phorgue/porousMultiphaseFoam/tree/openfoam-v8/libraries/porousModels by P. Horgue.
\*---------------------------------------------------------------------------*/

#ifndef densityModel_H
#define densityModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace densityModels
{
/*---------------------------------------------------------------------------*\
                   Class densityModel Declaration
\*---------------------------------------------------------------------------*/

class densityModel
{
// Protected data
protected:

    //- dictionary
    dictionary dict_;

    //- Reference to the concentration field
    const volScalarField& C_;

    //- Reference to the density field
    volScalarField& rho_;
    volScalarField& drhodC_;

    //- Disallow copy construct
    densityModel(const densityModel&);

    //- Disallow default bitwise assignment
    void operator=(const densityModel&);


public:

    //- Runtime type information
    TypeName("densityModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        densityModel,
        dictionary,
        (
            const dictionary& dict,
            const volScalarField& C,
            volScalarField& rho,
            volScalarField& drhodC
        ),
        (dict, C, rho, drhodC)
    );


    // Selectors

    //- Return a reference to the selected dispersion model
    static autoPtr<densityModel> New
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& rho,
        volScalarField& drhodC
    );


    // Constructors

    //- Construct from components
    densityModel
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& rho,
        volScalarField& drhodC
    );

    //- Destructor
    virtual ~densityModel()
        {}


    // Member Functions

    //- True unless trivial constant law
    virtual bool needsUpdate()
        {
            return 1;
        }

    //- Return the phase transport properties dictionary
    const dictionary& dict() const
        {
            return dict_;
        }


    //- Return density
    tmp<scalarField> rho();
    tmp<scalarField> rho(const volScalarField& C);
    virtual tmp<scalarField> rho(const scalarField& C) = 0;

    //- Return density derivative
    tmp<scalarField> drhodC();
    tmp<scalarField> drhodC(const volScalarField& C);
    virtual tmp<scalarField> drhodC(const scalarField& C) = 0;

    //- correct density including BC
    void correctRho();
    void correctDrhodC();

    //- Write
    virtual void write(Ostream& o) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace densityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
