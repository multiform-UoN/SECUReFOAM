Info << endl << "Reading field c" << endl;
volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << endl << "Reading transportProperties" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

#include "createRhoMuModels.H"

//- Read reference h if present
Info << endl << "Reading hRef" << endl;
const uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimLength, 0)
);

Info << endl << "Calculating field g.h\n" << endl;
const dimensionedScalar ghRef(- mag(g)*hRef);
const volScalarField gh("gh", (g & mesh.C()) - ghRef);
const surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

// - choose if forcing incompressibility (remove storativity and drho/dt terms)
const bool incompressible
(
  transportProperties.lookupOrDefault<bool>("incompressible",true)
);
Info<< "Incompressibility is " << incompressible << nl << endl;

// - choose if using standard weak cmopressibility (time derivative of density)
const bool standardCompressibility
(
  transportProperties.lookupOrDefault<bool>("standardCompressibility",true)
);
Info << endl << "Standard compressibility is " << standardCompressibility << endl;

// - choose if enabling hydrostatic contribution to the pressure
const bool hydrostaticPressure
(
  transportProperties.lookupOrDefault<bool>("hydrostaticPressure",true)
);
Info << endl << "hydrostaticPressure is " << hydrostaticPressure << endl;

// // - choose if enabling dynamic hydrostatic pressure (standard p_rgh in OF)
// const bool constantHydrostaticPressure
// (
//   transportProperties.lookupOrDefault<bool>("dynamicHydrostaticPressure",true)
// );
// Info<< "dynamicHydrostaticPressure is " << dynamicHydrostaticPressure << nl << endl;

Info << endl << "Reading field p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << endl << "Reading field p_rgh" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
// p_rgh = p - hydrostaticPressure*rho*gh;

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(p.name());

Info << endl << "Reading porosity beta" << endl;
volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
      transportProperties.lookupOrDefault<dimensionedScalar>
      (
        "beta", dimensionedScalar(dimless,scalar(1.))
      )
    )
);


Info << endl << "Reading storativity Ss" << endl;
volScalarField Ss
(
 IOobject
 (
  "Ss",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  transportProperties.lookupOrDefault<dimensionedScalar>
  (
   "Storativity", dimensionedScalar(dimless/dimPressure,scalar(0.))
  )
 )
);

//- Read permeability
// ORIGINAL
/*Info << endl << "Reading permeability field K" << endl;
const volSymmTensorField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
      "K",
      dimArea,
      transportProperties.lookupOrDefault<scalar>("K",scalar(1.))
      *
      symmTensor::I
    )
);*/
// EUGENIO
Info << endl << "Reading permeability field K" << endl;
const volSymmTensorField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceSymmTensorField Mf
(
    IOobject
    (
        "Mf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    inv
    (
        reverseLinear<SymmTensor<scalar>>
            (
                mesh
            ).interpolate(inv(K/mu))
    )
);

surfaceScalarField MfIso("MfIso",tr(Mf)/scalar(3));
surfaceSymmTensorField MfAniso("MfAniso",dev(Mf));

Info << endl << "Calculating face flux field phi" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (
        (
          Mf
          & (
                fvc::interpolate
                (
                    fvc::grad(p_rgh)
                )
            )
        )
    )
   & mesh.Sf()
);

//- Create density-driven flux
surfaceScalarField phiG
(
    "phiG",
    phi
);

Info << endl << "Calculating field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::reconstruct(phi)
);

// - Dimensional unit to use fvOptions for the pressure
dimensionedScalar dims(dimTime*dimTime/dimArea,1.0);
dimensionedScalar USMALL(dimLength/dimTime,SMALL);

//- Diffusion
const dimensionedScalar  Dm(transportProperties.lookup("Dm"));

//- Diffusion tensoor
volSymmTensorField  D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    Dm * symmTensor::I
);

//- Dispersion model
const bool dispersionModel
(
  transportProperties.lookupOrDefault<bool>("dispersionModel",true)
);
Info << endl << "dispersionModel is " << dispersionModel;

//- Transversal dispersion length scale
volScalarField alphaT
(
    IOobject
    (
        "alphaT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
      transportProperties.lookupOrDefault<dimensionedScalar>
      (
        "alphaT", dimensionedScalar(dimLength, scalar(0))
      )
    )
);

//- Longitudinal dispersion length scale
volScalarField alphaL
(
    IOobject
    (
        "alphaL",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
      transportProperties.lookupOrDefault<dimensionedScalar>
      (
        "alphaL", dimensionedScalar(dimLength, scalar(0))
      )
    )
);

surfaceSymmTensorField Df
(
    IOobject
    (
        "Df",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    inv
    (
        reverseLinear<SymmTensor<scalar>>
            (
                mesh
            ).interpolate(inv(D))
    )  // harmonic average of diffusivity
);

#include "updateD.H"

Info << endl << "Mean dispersion tensor is "
     << (fvc::domainIntegrate(D)/gSum(mesh.V())).value() << endl;

#include "createFvOptions.H"
