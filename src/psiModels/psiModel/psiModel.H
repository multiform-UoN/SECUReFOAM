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
    Foam::psiModel

Description
    A namespace for various psi implementations.

Class
    Foam::psiModel

Description
    An abstract base class for psi models.

SourceFiles
    psiModel.C
    psiModelNew.C

Acknowledgements

  Based on the porusModels from https://github.com/phorgue/porousMultiphaseFoam/tree/openfoam-v8/libraries/porousModels by P. Horgue.
\*---------------------------------------------------------------------------*/

#ifndef psiModel_H
#define psiModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace psiModels
{
/*---------------------------------------------------------------------------*\
                   Class psiModel Declaration
\*---------------------------------------------------------------------------*/

class psiModel
{
// Protected data
protected:

    //- dictionary
    dictionary dict_;

    //- Reference to the concentration field
    const volScalarField& Sw_;

    //- Reference to the psi, kappa fields
    volScalarField& psi_;
    volScalarField& dpsidSw_;
    volScalarField& kappa_;

    //- Disallow copy construct
    psiModel(const psiModel&);

    //- Disallow default bitwise assignment
    void operator=(const psiModel&);


public:

    //- Runtime type information
    TypeName("psiModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        psiModel,
        dictionary,
        (
            const dictionary& dict,
            const volScalarField& Sw,
            volScalarField& psi,
            volScalarField& dpsidSw,            
            volScalarField& kappa
        ),
        (dict, Sw, psi, dpsidSw, kappa)
    );


    // Selectors

    //- Return a reference to the selected psi model
    static autoPtr<psiModel> New
    (
        const dictionary& dict,
        const volScalarField& Sw,
        volScalarField& psi,
        volScalarField& dpsidSw,
        volScalarField& kappa
    );


    // Constructors

    //- Construct from components
    psiModel
    (
        const dictionary& dict,
        const volScalarField& Sw,
        volScalarField& psi,
        volScalarField& dpsidSw,
        volScalarField& kappa
    );

    //- Destructor
    virtual ~psiModel()
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


    //- Return psi, dpsidSw, kappa
    tmp<scalarField> psi();
    tmp<scalarField> psi(const volScalarField& Sw);
    virtual tmp<scalarField> psi(const scalarField& Sw) = 0;
    tmp<scalarField> dpsidSw();
    tmp<scalarField> dpsidSw(const volScalarField& Sw);
    virtual tmp<scalarField> dpsidSw(const scalarField& Sw) = 0;
    tmp<scalarField> kappa();
    tmp<scalarField> kappa(const volScalarField& Sw);
    virtual tmp<scalarField> kappa(const scalarField& Sw) = 0;

    //- correct psi including BC
    void correctPsi();
    void correctDpsidSw();
    void correctKappa();

    //- Write
    virtual void write(Ostream& o) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace psiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
