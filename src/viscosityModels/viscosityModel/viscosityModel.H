/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Namespace
    Foam::viscosityModel

Description
    A namespace for various viscosity implementations.

Class
    Foam::viscosityModel

Description
    An abstract base class for viscosity models.

SourceFiles
    viscosityModel.C
    viscosityModelNew.C
\*---------------------------------------------------------------------------*/

#ifndef viscosityModel_H
#define viscosityModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
/*---------------------------------------------------------------------------*\
                   Class viscosityModel Declaration
\*---------------------------------------------------------------------------*/

class viscosityModel
{
// Protected data
protected:

    //- dictionary
    dictionary dict_;

    //- Reference to the concentration field
    const volScalarField& C_;

    //- Reference to the viscosity field
    volScalarField& mu_;

    //- Disallow copy construct
    viscosityModel(const viscosityModel&);

    //- Disallow default bitwise assignment
    void operator=(const viscosityModel&);


public:

    //- Runtime type information
    TypeName("viscosityModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        viscosityModel,
        dictionary,
        (
            const dictionary& dict,
            const volScalarField& C,
            volScalarField& muField
        ),
        (dict, C, muField)
    );


    // Selectors

    //- Return a reference to the selected dispersion model
    static autoPtr<viscosityModel> New
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& muField
    );


    // Constructors

    //- Construct from components
    viscosityModel
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& muField
    );

    //- Destructor
    virtual ~viscosityModel()
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


    //- Return viscosity
    tmp<scalarField> mu();
    tmp<scalarField> mu(const volScalarField& C);
    virtual tmp<scalarField> mu(const scalarField& C) = 0;

    //- correct viscosity including BC
    void correctMu();


    //- Write
    virtual void write(Ostream& o) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
