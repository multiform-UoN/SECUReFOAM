Info << endl << "Reading field cA" << endl;
volScalarField cA
(
    IOobject
    (
        "cA",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << endl << "Reading field cB" << endl;
volScalarField cB
(
    IOobject
    (
     "cB",
     runTime.timeName(),
     mesh,
     IOobject::MUST_READ,
     IOobject::AUTO_WRITE
     ),
    mesh
 );

Info << endl << "Reading field cC" << endl;
volScalarField cC
(
    IOobject
    (
     "cC",
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


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


Info << endl << "Reading Damköhler number" << endl;
volScalarField Dam
(
    IOobject
    (
     "Dam",
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
      "Dam", dimensionedScalar(dimVol/dimTime/dimMass,scalar(1.))
      )
     )
 );


Info << endl << "Reading Da" << endl;
const dimensionedScalar  Da(transportProperties.lookup("Da"));
//- Diffusion tensor
// volSymmTensorField  Da
// (
//     IOobject
//     (
//         "Da",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     deltaa * symmTensor::I
// );
// Info << endl << "Reading deltab=Db/Da" << endl;
Info << endl << "Reading Db" << endl;
const dimensionedScalar  Db(transportProperties.lookup("Da"));
//- Diffusion tensor
// volSymmTensorField  Db
// (
//     IOobject
//     (
//         "Db",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     deltab * symmTensor::I
// );
// Info << endl << "Reading deltac=Dc/Da" << endl;

Info << endl << "Reading Dc" << endl;
const dimensionedScalar  Dc(transportProperties.lookup("Dc"));
//- Diffusion tensor
// volSymmTensorField  Dc
// (
//     IOobject
//     (
//         "Dc",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     deltac * symmTensor::I
// );

//- Read permeability
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
    mesh,
    dimensionedSymmTensor
    (
      "K",
      dimArea,
      transportProperties.lookupOrDefault<scalar>("K",scalar(1.))
      *
      symmTensor::I
    )
);

surfaceSymmTensorField Mf
(
    IOobject
    (
        "Mf",
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
        IOobject::AUTO_WRITE
    ),
    (
        (
          Mf
          & (
                fvc::interpolate
                (
                    fvc::grad(p)
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
dimensionedScalar dims(dimTime*dimTime*dimLength/dimMass,1.0);
//dimensionedScalar dims(dimTime*dimTime/dimArea,1.0);
dimensionedScalar USMALL(dimLength/dimTime,SMALL);

#include "createFvOptions.H"
