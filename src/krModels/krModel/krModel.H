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
    Foam::krModel

Description
    A namespace for various kr implementations.

Class
    Foam::krModel

Description
    An abstract base class for kr models.

SourceFiles
    krModel.C
    krModelNew.C

Acknowledgements

  Based on the porusModels from https://github.com/phorgue/porousMultiphaseFoam/tree/openfoam-v8/libraries/porousModels by P. Horgue.
\*---------------------------------------------------------------------------*/

#ifndef krModel_H
#define krModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace krModels
{
/*---------------------------------------------------------------------------*\
                   Class krModel Declaration
\*---------------------------------------------------------------------------*/

class krModel
{
// Protected data
protected:

    //- dictionary
    dictionary dict_;

    //- Reference to the concentration field
    const volScalarField& Sw_;

    //- Reference to the kr field
    volScalarField& kr_;

    //- Disallow copy construct
    krModel(const krModel&);

    //- Disallow default bitwise assignment
    void operator=(const krModel&);


public:

    //- Runtime type information
    TypeName("krModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        krModel,
        dictionary,
        (
            const dictionary& dict,
            const volScalarField& Sw,
            volScalarField& kr
        ),
        (dict, Sw, kr)
    );


    // Selectors

    //- Return a reference to the selected kr model
    static autoPtr<krModel> New
    (
        const dictionary& dict,
        const volScalarField& Sw,
        volScalarField& kr
    );


    // Constructors

    //- Construct from components
    krModel
    (
        const dictionary& dict,
        const volScalarField& Sw,
        volScalarField& kr
    );

    //- Destructor
    virtual ~krModel()
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


    //- Return kr
    tmp<scalarField> kr();
    tmp<scalarField> kr(const volScalarField& Sw);
    virtual tmp<scalarField> kr(const scalarField& Sw) = 0;


    //- correct kr including BC
    void correctKr();

    //- Write
    virtual void write(Ostream& o) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace krModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
