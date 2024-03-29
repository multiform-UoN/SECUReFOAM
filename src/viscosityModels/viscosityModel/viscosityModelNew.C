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

#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::viscosityModels::viscosityModel>
Foam::viscosityModels::viscosityModel::New
(
    const dictionary& dict,
    const volScalarField& C,
    volScalarField& mu
)
{
    const word modelType(dict.lookup("type"));

    Info << endl << "Selecting viscosity model => " << modelType << "\n" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
            (
                "viscosityModel::New(const word&, "
                " const dictionary& transportProperties,"
                " const volScalarField& C"
            )   << "Unknown viscosityModel type "
                << modelType << nl << nl
                << "Valid viscosityModels are : " << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
    }

    return autoPtr<viscosityModel>
        (cstrIter()(dict, C, mu));
}


// ************************************************************************* //
