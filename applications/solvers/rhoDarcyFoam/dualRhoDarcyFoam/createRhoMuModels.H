// Viscosity
Info << endl << "Creating viscosity mu_fr" << endl;
volScalarField mu_fr
(
    IOobject
    (
        "mu_fr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    transportProperties.lookupOrDefault<dimensionedScalar>
    (
        "mu_fr", dimensionedScalar(dimensionSet(1,-1,-1,0,0),scalar(1))
    )
);

autoPtr<viscosityModels::viscosityModel> muModel_fr
        (
          viscosityModels::viscosityModel::New
          (
            transportProperties.subDict("viscosityModel"),
            c_fr,
            mu_fr
          )
        );

#include "updateMu.H"

//- Densities
Info << endl << "Creating density rho_fr" << endl;
volScalarField rho_fr
(
    IOobject
    (
        "rho_fr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    transportProperties.lookupOrDefault<dimensionedScalar>
    (
        "rho_fr", dimensionedScalar(dimensionSet(1,-3,0,0,0),scalar(1))
    )
);

volScalarField drhodc_fr
(
    IOobject
    (
        "drhodc_fr",
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

autoPtr<densityModels::densityModel> rhoModel_fr
        (
          densityModels::densityModel::New
          (
            transportProperties.subDict("densityModel"),
            c_fr,
            rho_fr,
            drhodc_fr
          )
        );

#include "updateRho.H"