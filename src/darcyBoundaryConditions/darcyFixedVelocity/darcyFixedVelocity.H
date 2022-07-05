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
    Foam::darcyFixedVelocity

Description
    Set the pressure gradient boundary condition for a given Darcy total velocity
    Permeability is anisotropic (K volTensorField)

SourceFiles
    darcyFixedVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef darcyFixedVelocity_H
#define darcyFixedVelocity_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class darcyGradPressureFvPatch Declaration
    \*---------------------------------------------------------------------------*/

class darcyFixedVelocity
:
    public fixedGradientFvPatchScalarField
{
    // Private data

    word MfName_;
    scalar velocity_;

public:

    //- Runtime type information
    TypeName("darcyFixedVelocity");


    // Constructors

    //- Construct from patch and internal field
    darcyFixedVelocity
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    darcyFixedVelocity
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given
    //  darcyFixedVelocity onto a new patch
    darcyFixedVelocity
    (
        const darcyFixedVelocity&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    darcyFixedVelocity
    (
        const darcyFixedVelocity&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
                (
                    new darcyFixedVelocity(*this)
                );
        }

    //- Construct as copy setting internal field reference
    darcyFixedVelocity
    (
        const darcyFixedVelocity&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
        {
            return tmp<fvPatchScalarField>
                (
                    new darcyFixedVelocity(*this, iF)
                );
        }


    // Member functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
