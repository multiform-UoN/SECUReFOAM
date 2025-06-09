//- Viscosity
Info << endl << "Creating viscosity mu" << endl;
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    transportProperties.lookupOrDefault<dimensionedScalar>
    (
        "mu", dimensionedScalar(dimensionSet(1,-1,-1,0,0),scalar(1))
    )
);

autoPtr<viscosityModels::viscosityModel> muModel
        (
          viscosityModels::viscosityModel::New
          (
            transportProperties.subDict("viscosityModel"),
            c,
            mu
          )
        );

#include "updateMu.H"

bool updateMf = muModel->needsUpdate();
//- Densities
Info << endl << "Creating density rho" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    transportProperties.lookupOrDefault<dimensionedScalar>
    (
        "rho", dimensionedScalar(dimensionSet(1,-3,0,0,0),scalar(1))
    )
);

volScalarField drhodc
(
    IOobject
    (
        "drhodc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
      rho.dimensions()/c.dimensions(),
      scalar(0)
    )
);

autoPtr<densityModels::densityModel> rhoModel
        (
          densityModels::densityModel::New
          (
            transportProperties.subDict("densityModel"),
            c,
            rho,
            drhodc
          )
        );

#include "updateRho.H"

bool updatePhiG = rhoModel->needsUpdate();
