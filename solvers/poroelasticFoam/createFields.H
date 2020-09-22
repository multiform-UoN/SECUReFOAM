Info<< "Reading field p\n" << endl;
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

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading permeability field K\n" << endl;
volSymmTensorField K
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

//- gradient of displacement
volTensorField gradD(fvc::grad(D));

//- Jacobian
volScalarField J(scalar(1.0) + fvc::div(D));

Info<< "Create stress field Sigma\n" << endl;
volSymmTensorField Sigma
(
    IOobject
    (
        "Sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mu*symm(gradD + T(gradD)) + lambda*tr(gradD)*I
);

//- Create mobility surface tensor field
surfaceSymmTensorField Mf("Mf",fvc::interpolate(K,"K")/(nu*rhof));

//- Create gravitational flux
surfaceScalarField phiG("phiG",(rhof * (Mf & g)) & mesh.Sf());

#include "createFvOptions.H"

// - Dimensional unit to use fvOptions for the pressure
dimensionedScalar dims(Mf.dimensions()*dimTime/dimArea,1.0);
