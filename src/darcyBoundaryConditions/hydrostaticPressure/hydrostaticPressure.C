/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hydrostaticPressure.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hydrostaticPressure::
hydrostaticPressure
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}



Foam::hydrostaticPressure::
hydrostaticPressure
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}



Foam::hydrostaticPressure::
hydrostaticPressure
(
    const hydrostaticPressure& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}



Foam::hydrostaticPressure::
hydrostaticPressure
(
    const hydrostaticPressure& ptf
)
:
    fixedValueFvPatchScalarField(ptf)
{}



Foam::hydrostaticPressure::
hydrostaticPressure
(
    const hydrostaticPressure& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::hydrostaticPressure::
updateCoeffs()
{
    if (fixedValueFvPatchScalarField::updated())
    {
        return;
    }

    fixedValueFvPatchScalarField::updateCoeffs();

    const scalarField& rhop = this->patch().template
        lookupPatchField<volScalarField, scalar>
        (
            "rho"
        );

    const scalarField& gh = this->patch().template
        lookupPatchField<volScalarField, scalar>
        (
            "rho"
        );

    this->operator==
    (
        *this - rhop*gh
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    hydrostaticPressure
);
}


// ************************************************************************* //
