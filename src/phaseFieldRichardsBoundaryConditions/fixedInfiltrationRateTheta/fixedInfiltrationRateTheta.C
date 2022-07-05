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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fixedInfiltrationRateTheta.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedInfiltrationRateTheta::fixedInfiltrationRateTheta
(
    const fvPatch& Sw,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(Sw, iF),
    MfName_("Mf"),
    dpsidSwName_("dpsidSw"),
    SwName_("Sw"),
    infilRate_(scalar(0))
{}

Foam::fixedInfiltrationRateTheta::fixedInfiltrationRateTheta
(
    const fvPatch& Sw,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedGradientFvPatchScalarField(Sw, iF),
    MfName_(dict.lookupOrDefault<word>("Mf", "Mf")),
    dpsidSwName_(dict.lookupOrDefault<word>("dpsidSw", "dpsidSw")),
    SwName_(dict.lookupOrDefault<word>("Sw", "Sw")),
    infilRate_(readScalar(dict.lookup("infilRate")))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}

Foam::fixedInfiltrationRateTheta::fixedInfiltrationRateTheta
(
    const fixedInfiltrationRateTheta& ptf,
    const fvPatch& Sw,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedGradientFvPatchScalarField(ptf, Sw, iF, mapper),
    MfName_(ptf.MfName_),
    dpsidSwName_(ptf.dpsidSwName_),
    SwName_(ptf.SwName_),
    infilRate_(ptf.infilRate_)
{}

Foam::fixedInfiltrationRateTheta::fixedInfiltrationRateTheta
(
    const fixedInfiltrationRateTheta& ptf
)
    :
    fixedGradientFvPatchScalarField(ptf),
    MfName_(ptf.MfName_),
    dpsidSwName_(ptf.dpsidSwName_),
    SwName_(ptf.SwName_),
    infilRate_(ptf.infilRate_)
{}

Foam::fixedInfiltrationRateTheta::fixedInfiltrationRateTheta
(
    const fixedInfiltrationRateTheta& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(ptf, iF),
    MfName_(ptf.MfName_),
    dpsidSwName_(ptf.dpsidSwName_),
    SwName_(ptf.SwName_),
    infilRate_(ptf.infilRate_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedInfiltrationRateTheta::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //fvc::interpolate(dpsidSw)
    if(db().foundObject<volScalarField>("theta"))
    {

        const fvsPatchField<symmTensor>& Mf=
            patch().lookupPatchField<surfaceSymmTensorField, tensor>(MfName_);

        const fvPatchField<scalar>& dpsidSw=
            patch().lookupPatchField<volScalarField, scalar>(dpsidSwName_);

        const fvPatchField<scalar>& Sw =
            patch().lookupPatchField<volScalarField, scalar>(SwName_);

        const uniformDimensionedVectorField& g
            (
                    db().lookupObject<IOobject>("g")
            );

            gradient() =
                    -
                    infilRate_ * ( patch().nf()  & ( inv(Mf) &  patch().nf() ) )
                    - 
                    (
                        (  g.value() & patch().nf() ) +  ( dpsidSw * Sw.snGrad()  )
                    );

	}

    else
    {
        FatalErrorInFunction
            <<"This BC can only be applied to theta"
            << exit(FatalError);
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::fixedInfiltrationRateTheta::write(Ostream& os) const
{
    writeEntryIfDifferent<word>(os, "Mf", "Mf", MfName_);
    writeEntryIfDifferent<word>(os, "dpsidSw", "dpsidSw", dpsidSwName_);
    writeEntryIfDifferent<word>(os, "Sw", "Sw", SwName_);
    writeEntry<scalar>(os, "infilRate", infilRate_);
    fixedGradientFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    fixedInfiltrationRateTheta
);
}


// ************************************************************************* //
