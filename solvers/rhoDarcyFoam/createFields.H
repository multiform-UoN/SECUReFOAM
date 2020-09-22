Info << nl << "Reading field c\n" << endl;
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


Info << nl << "Reading field p\n" << endl;
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


Info << "Reading transportProperties\n" << endl;
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


//- These could become volScalarFields
const dimensionedScalar  rho1(transportProperties.lookup("rho1"));
const dimensionedScalar  rho2(transportProperties.lookup("rho2"));

const dimensionedScalar  mu(transportProperties.lookup("mu"));
const dimensionedScalar  Dm(transportProperties.lookup("Dm"));

const dimensionedScalar D(Dm);
// const dimensionedTensor D0((tensor::I));
//volTensorField D("D", D0*(mag(U)*sqrt(mag(K)) + Dm) );

// - choose if enabling hydrostatic contribution to the pressure
const bool hydrostaticPressure
(
  transportProperties.lookupOrDefault<bool>("hydrostaticPressure",true)
);
Info<< "hydrostaticPressure is " << hydrostaticPressure << nl << endl;

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
    rho2*c + rho1*(scalar(1.0)-c)
);


//- Read reference h if present
Info<< "\nReading hRef\n" << endl;
uniformDimensionedScalarField hRef
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

Info<< "Calculating field g.h\n" << endl;
dimensionedScalar ghRef(- mag(g)*hRef);
const volScalarField gh("gh", (g & mesh.C()) - ghRef);
const surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

Info<< "Reading field p_rgh\n" << endl;
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
p_rgh = p - hydrostaticPressure*rho*gh;

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

mesh.setFluxRequired(p_rgh.name());

Info << nl << "Reading porosity beta\n" << endl;
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
    dimensionedScalar( "beta0", dimless, scalar(1.) )
);

//- Read permeability
Info << "Reading field K\n" << endl;
volTensorField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


//- Create mobility surface tensor field
const surfaceTensorField Mf("Mf",fvc::interpolate(K,"K")/mu);
const surfaceScalarField MfIso("Mf",tr(Mf)/scalar(3));
const surfaceTensorField MfAniso("Mf",dev(Mf));

Info<< "Calculating face flux field phi\n" << endl;
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

Info << nl << "Calculating field U\n" << endl;
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

// #include "phiG.H"

#include "createFvOptions.H"

// - Dimensional unit to use fvOptions for the pressure
dimensionedScalar dims(Mf.dimensions()*dimTime/dimArea,1.0);
