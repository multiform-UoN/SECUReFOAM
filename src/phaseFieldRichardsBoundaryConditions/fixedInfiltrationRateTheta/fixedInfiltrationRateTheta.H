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
    Foam::fixedInfiltrationRateTheta

Description
    Set the saturation gradient boundary condition for a given infiltration rate

SourceFiles
    fixedInfiltrationRateTheta.C

\*---------------------------------------------------------------------------*/

#ifndef fixedInfiltrationRateTheta_H
#define fixedInfiltrationRateTheta_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class fixedInfiltrationRateThetaFvPatch Declaration
    \*---------------------------------------------------------------------------*/

class fixedInfiltrationRateTheta
:
    public fixedGradientFvPatchScalarField
{
    // Private data

    word MfName_;
    word dpsidSwName_;
    word SwName_;
    scalar infilRate_;

public:

    //- Runtime type information
    TypeName("fixedInfiltrationRateTheta");


    // Constructors

    //- Construct from patch and internal field
    fixedInfiltrationRateTheta
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    fixedInfiltrationRateTheta
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given
    //  fixedInfiltrationRateTheta onto a new patch
    fixedInfiltrationRateTheta
    (
        const fixedInfiltrationRateTheta&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    fixedInfiltrationRateTheta
    (
        const fixedInfiltrationRateTheta&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
                (
                    new fixedInfiltrationRateTheta(*this)
                );
        }

    //- Construct as copy setting internal field reference
    fixedInfiltrationRateTheta
    (
        const fixedInfiltrationRateTheta&,
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
                    new fixedInfiltrationRateTheta(*this, iF)
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
