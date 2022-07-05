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
    Foam::phaseFieldModel

Description
    A namespace for various psi implementations.

Class
    Foam::phaseFieldModel

Description
    An abstract base class for psi models.

SourceFiles
    phaseFieldModel.C
    phaseFieldModelNew.C

Acknowledgements

  Based on the porusModels from https://github.com/phorgue/porousMultiphaseFoam/tree/openfoam-v8/libraries/porousModels by P. Horgue.
\*---------------------------------------------------------------------------*/

#ifndef phaseFieldModel_H
#define phaseFieldModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFieldModels
{
/*---------------------------------------------------------------------------*\
                   Class phaseFieldModel Declaration
\*---------------------------------------------------------------------------*/

class phaseFieldModel
{
// Protected data
protected:

    //- dictionary
    dictionary dict_;

    //- Reference to the concentration field
    const volScalarField& C_;

    //- Reference to the chemical potential
    volScalarField& mu_;

    //- Reference to the mobility
    volScalarField& kappa_;

    //- Disallow copy construct
    phaseFieldModel(const phaseFieldModel&);

    //- Disallow default bitwise assignment
    void operator=(const phaseFieldModel&);


public:

    //- Runtime type information
    TypeName("phaseFieldModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseFieldModel,
        dictionary,
        (
            const dictionary& dict,
            const volScalarField& C,
            volScalarField& mu,
            volScalarField& kappa
        ),
        (dict, C, mu, kappa)
    );


    // Selectors

    //- Return a reference to the selected psi model
    static autoPtr<phaseFieldModel> New
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& mu,
        volScalarField& kappa
    );


    // Constructors

    //- Construct from components
    phaseFieldModel
    (
        const dictionary& dict,
        const volScalarField& C,
        volScalarField& mu,
        volScalarField& kappa
    );

    //- Destructor
    virtual ~phaseFieldModel()
        {}


    // Member Functions

    //- Return the phase transport properties dictionary
    const dictionary& dict() const
        {
            return dict_;
        }


    //- Return values
    tmp<scalarField> psi();
    tmp<scalarField> psi(const volScalarField& C);
    virtual tmp<scalarField> psi(const scalarField& C) = 0;
    tmp<scalarField> dpsidC();
    tmp<scalarField> dpsidC(const volScalarField& C);
    virtual tmp<scalarField> dpsidC(const scalarField& C) = 0;
    tmp<scalarField> mu();
    tmp<scalarField> mu(const volScalarField& C);
    virtual tmp<scalarField> mu(const scalarField& C) = 0;
    tmp<scalarField> kappa();
    tmp<scalarField> kappa(const volScalarField& C);
    virtual tmp<scalarField> kappa(const scalarField& C) = 0;
    tmp<scalarField> freeEnergy();
    tmp<scalarField> freeEnergy(const volScalarField& C);
    virtual tmp<scalarField> freeEnergy(const scalarField& C) = 0;

    //- correct psi including BC
    void correctMu();
    void correctKappa();

    //- Write
    virtual void write(Ostream& o) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseFieldModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
