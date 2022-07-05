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
#include "darcyGradPressureAniso.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_("Mf"),
    phiName_("phi")
{}

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_(dict.lookupOrDefault<word>("Mf", "Mf")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const darcyGradPressureAniso& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_)
{}

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const darcyGradPressureAniso& ptf
)
    :
    fixedGradientFvPatchScalarField(ptf),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_)
{}

Foam::darcyGradPressureAniso::darcyGradPressureAniso
(
    const darcyGradPressureAniso& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(ptf, iF),
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyGradPressureAniso::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<symmTensor>& Mf=
        patch().lookupPatchField<surfaceSymmTensorField, tensor>(MfName_);

    const fvPatchField<vector>& U=
        patch().lookupPatchField<volVectorField, vector>("U");


    if(db().foundObject<volScalarField>("p_rgh"))
    {

        const fvsPatchField<scalar>& ghf=
                patch().lookupPatchField<surfaceScalarField, scalar>("ghf");

        const fvPatchField<scalar>& rho =
                patch().lookupPatchField<volScalarField, scalar>("rho");

        gradient() =
        - ( (inv(Mf) & U) & patch().nf())
        - ( rho.snGrad() * ghf );
    }
    else if(db().foundObject<volScalarField>("p"))
    {
        const uniformDimensionedVectorField& g
        (
            db().lookupObject<IOobject>("g")
        );

        dimensionedScalar rho
        (
            db().lookupObject<IOdictionary>
            (
                "transportProperties"
            ).lookup("rho")
        );

        gradient() =
            - ( (inv(Mf) & U)& patch().nf())
            + ( rho.value() * g.value() & patch().nf() );
    }
    else
    {
        FatalErrorInFunction
            <<"This BC can only be applied to p or p_rgh"
            << exit(FatalError);
    }


    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::darcyGradPressureAniso::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "Mf", "Mf", MfName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    darcyGradPressureAniso
);
}


// ************************************************************************* //
