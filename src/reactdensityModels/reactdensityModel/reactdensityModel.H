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
    Foam::reactdensityModel

Description
    A namespace for various density implementations.

Class
    Foam::reactdensityModel

Description
    An abstract base class for density models in reactive transport problems.

SourceFiles
    reactdensityModel.C
    reactdensityModelNew.C

Acknowledgements

  Based on the porusModels from https://github.com/phorgue/porousMultiphaseFoam/tree/openfoam-v8/libraries/porousModels by P. Horgue.
\*---------------------------------------------------------------------------*/

#ifndef reactdensityModel_H
#define reactdensityModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactdensityModels
{
/*---------------------------------------------------------------------------*\
                   Class reactdensityModel Declaration
\*---------------------------------------------------------------------------*/

class reactdensityModel
{
// Protected data
protected:

    //- dictionary
    dictionary dict_;

    //- Reference to the concentration fields
    const volScalarField& cA_;
    const volScalarField& cB_;
    const volScalarField& cC_;
  
    //- Reference to the density field
    volScalarField& rho_;
    volScalarField& drhodC_;

    //- Disallow copy construct
    reactdensityModel(const reactdensityModel&);

    //- Disallow default bitwise assignment
    void operator=(const reactdensityModel&);


public:

    //- Runtime type information
    TypeName("reactdensityModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        reactdensityModel,
        dictionary,
        (
            const dictionary& dict,
            const volScalarField& cA,
	    const volScalarField& cB,
	    const volScalarField& cC,
            volScalarField& rho,
            volScalarField& drhodC
        ),
        (dict, cA, cB, cC, rho, drhodC)
    );


    // Selectors

    //- Return a reference to the selected dispersion model
    static autoPtr<reactdensityModel> New
    (
        const dictionary& dict,
        const volScalarField& cA,
        const volScalarField& cB,
	const volScalarField& cC,
        volScalarField& rho,
        volScalarField& drhodC
    );


    // Constructors

    //- Construct from components
    reactdensityModel
    (
        const dictionary& dict,
        const volScalarField& cA,
	const volScalarField& cB,	
        const volScalarField& cC,
        volScalarField& rho,
        volScalarField& drhodC
    );

    //- Destructor
    virtual ~reactdensityModel()
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
    tmp<scalarField> rho(const volScalarField& cA,const volScalarField& cB,const volScalarField& cC );
    virtual tmp<scalarField> rho(const scalarField& cA, const scalarField& cB, const scalarField& cC) = 0;

    //- Return density derivative
    tmp<scalarField> drhodC();
    tmp<scalarField> drhodC(const volScalarField& cA);
    virtual tmp<scalarField> drhodC(const scalarField& cA) = 0;

    //- correct density including BC
    void correctRho();
    void correctDrhodC();

    //- Write
    virtual void write(Ostream& o) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactdensityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
