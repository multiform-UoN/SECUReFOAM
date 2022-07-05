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
#include "fixedInfiltrationRate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedInfiltrationRate::fixedInfiltrationRate
(
    const fvPatch& Sw,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(Sw, iF),
    MfName_("Mf"),
    dpsidSwName_("dpsidSw"),
    ThetaName_("theta"),
    infilRate_(scalar(0))
{}

Foam::fixedInfiltrationRate::fixedInfiltrationRate
(
    const fvPatch& Sw,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedGradientFvPatchScalarField(Sw, iF),
    MfName_(dict.lookupOrDefault<word>("Mf", "Mf")),
    dpsidSwName_(dict.lookupOrDefault<word>("dpsidSw", "dpsidSw")),
    ThetaName_(dict.lookupOrDefault<word>("theta", "theta")),
    infilRate_(readScalar(dict.lookup("infilRate")))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}

Foam::fixedInfiltrationRate::fixedInfiltrationRate
(
    const fixedInfiltrationRate& ptf,
    const fvPatch& Sw,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedGradientFvPatchScalarField(ptf, Sw, iF, mapper),
    MfName_(ptf.MfName_),
    dpsidSwName_(ptf.dpsidSwName_),
    ThetaName_(ptf.ThetaName_),
    infilRate_(ptf.infilRate_)
{}

Foam::fixedInfiltrationRate::fixedInfiltrationRate
(
    const fixedInfiltrationRate& ptf
)
    :
    fixedGradientFvPatchScalarField(ptf),
    MfName_(ptf.MfName_),
    dpsidSwName_(ptf.dpsidSwName_),
    ThetaName_(ptf.ThetaName_),
    infilRate_(ptf.infilRate_)
{}

Foam::fixedInfiltrationRate::fixedInfiltrationRate
(
    const fixedInfiltrationRate& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(ptf, iF),
    MfName_(ptf.MfName_),
    dpsidSwName_(ptf.dpsidSwName_),
    ThetaName_(ptf.ThetaName_),
    infilRate_(ptf.infilRate_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedInfiltrationRate::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if(db().foundObject<volScalarField>("Sw"))
    {

        const fvsPatchField<symmTensor>& Mf=
            patch().lookupPatchField<surfaceSymmTensorField, tensor>(MfName_);

        const fvPatchField<scalar>& dpsidSw=
            patch().lookupPatchField<volScalarField, scalar>(dpsidSwName_);

        const fvPatchField<scalar>& theta =
            patch().lookupPatchField<volScalarField, scalar>(ThetaName_);

        const uniformDimensionedVectorField& g
            (
                    db().lookupObject<IOobject>("g")
            );

            gradient() =
                    (
                    -
                    infilRate_ * ( patch().nf()  & ( inv(Mf) &  patch().nf() ) )
                    - 
                    
                        (  g.value() & patch().nf() ) +  (  theta.snGrad()  )
                    )
                    /dpsidSw;

	}

    else
    {
        FatalErrorInFunction
            <<"This BC can only be applied to Sw"
            << exit(FatalError);
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::fixedInfiltrationRate::write(Ostream& os) const
{
    writeEntryIfDifferent<word>(os, "Mf", "Mf", MfName_);
    writeEntryIfDifferent<word>(os, "dpsidSw", "dpsidSw", dpsidSwName_);
    writeEntryIfDifferent<word>(os, "theta", "theta", ThetaName_);
    writeEntry<scalar>(os, "infilRate", infilRate_);
    fixedGradientFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    fixedInfiltrationRate
);
}


// ************************************************************************* //
