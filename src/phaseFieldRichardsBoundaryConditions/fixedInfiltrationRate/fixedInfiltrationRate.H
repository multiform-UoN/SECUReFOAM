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
    Foam::fixedInfiltrationRate

Description
    Set the saturation gradient boundary condition for a given infiltration rate

SourceFiles
    fixedInfiltrationRate.C

\*---------------------------------------------------------------------------*/

#ifndef fixedInfiltrationRate_H
#define fixedInfiltrationRate_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class fixedInfiltrationRateFvPatch Declaration
    \*---------------------------------------------------------------------------*/

class fixedInfiltrationRate
:
    public fixedGradientFvPatchScalarField
{
    // Private data

    word MfName_;
    word dpsidSwName_;
    word ThetaName_;
    scalar infilRate_;

public:

    //- Runtime type information
    TypeName("fixedInfiltrationRate");


    // Constructors

    //- Construct from patch and internal field
    fixedInfiltrationRate
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    fixedInfiltrationRate
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given
    //  fixedInfiltrationRate onto a new patch
    fixedInfiltrationRate
    (
        const fixedInfiltrationRate&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    fixedInfiltrationRate
    (
        const fixedInfiltrationRate&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
                (
                    new fixedInfiltrationRate(*this)
                );
        }

    //- Construct as copy setting internal field reference
    fixedInfiltrationRate
    (
        const fixedInfiltrationRate&,
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
                    new fixedInfiltrationRate(*this, iF)
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
