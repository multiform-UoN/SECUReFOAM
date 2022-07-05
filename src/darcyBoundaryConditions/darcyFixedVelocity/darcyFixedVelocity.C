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
#include "darcyFixedVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyFixedVelocity::darcyFixedVelocity
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_("Mf"),
    velocity_(scalar(0))
{}

Foam::darcyFixedVelocity::darcyFixedVelocity
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    MfName_(dict.lookupOrDefault<word>("Mf", "Mf")),
    velocity_(readScalar(dict.lookup("velocity")))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}

Foam::darcyFixedVelocity::darcyFixedVelocity
(
    const darcyFixedVelocity& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    MfName_(ptf.MfName_),
    velocity_(ptf.velocity_)
{}

Foam::darcyFixedVelocity::darcyFixedVelocity
(
    const darcyFixedVelocity& ptf
)
    :
    fixedGradientFvPatchScalarField(ptf),
    MfName_(ptf.MfName_),
    velocity_(ptf.velocity_)
{}

Foam::darcyFixedVelocity::darcyFixedVelocity
(
    const darcyFixedVelocity& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(ptf, iF),
    MfName_(ptf.MfName_),
    velocity_(ptf.velocity_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyFixedVelocity::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<symmTensor>& Mf=
        patch().lookupPatchField<surfaceSymmTensorField, tensor>(MfName_);

    if(db().foundObject<volScalarField>("p_rgh"))
    {
        const fvPatchField<scalar>& rho =
                patch().lookupPatchField<volScalarField, scalar>("rho");

        const bool hydrostaticPressure
        (
        db().lookupObject<IOdictionary>
                (
                    "transportProperties"
                ).lookupOrDefault<bool>("hydrostaticPressure",true)
        );

        if (hydrostaticPressure)
        {
            // const bool dynamicHydrostaticPressure
            //     (
            //     db().lookupObject<IOdictionary>
            //             (
            //                 "transportProperties"
            //             ).lookupOrDefault<bool>("dynamicHydrostaticPressure",true)
            //     );
            // if (dynamicHydrostaticPressure)
            {
                const fvsPatchField<scalar>& ghf=
                        patch().lookupPatchField<surfaceScalarField, scalar>("ghf");

                gradient() =
                - ( (inv(Mf) & (-velocity_ * patch().nf())) & patch().nf() )
                - ( rho.snGrad() * ghf );
            }
            // else
            {
                // TODO
            }
        }
        else
        {
            const uniformDimensionedVectorField& g
            (
                db().lookupObject<IOobject>("g")
            );

            gradient() =
                - ( (inv(Mf) & (-velocity_ * patch().nf())) & patch().nf() )
                + ( rho * g.value() & patch().nf() );
        }
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
            - ( (inv(Mf) & (-velocity_ * patch().nf())) & patch().nf() )
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


void Foam::darcyFixedVelocity::write(Ostream& os) const
{
    writeEntry(os, "type", type());
    if (patchType().size())
       {
           writeEntry(os, "patchType", patchType());
       }
    writeEntryIfDifferent<word>(os, "Mf", "Mf", MfName_);
    writeEntry<scalar>(os, "velocity", velocity_);
    writeEntry(os,"value",*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    darcyFixedVelocity
);
}


// ************************************************************************* //
