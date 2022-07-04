//- Read pressure and velocity
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

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

//- Read fluid properties
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


const dimensionedScalar rho(transportProperties.lookup("rho"));
const dimensionedScalar  mu(transportProperties.lookup("mu"));;

//- Read permeability
Info << "Reading field K\n" << endl;
const volSymmTensorField K
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
const surfaceSymmTensorField Mf
(
    IOobject
    (
        "Mf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    // (fvc::interpolate((K)))
    inv
    (
        reverseLinear<SymmTensor<scalar>>
            (
                mesh
            ).interpolate(inv(K/mu))
    )  // harmonic average of permeabilities
);

//- Create gravitational flux
// NB: this is different than phiG in rhoDarcy as we don't have p_rgh here
surfaceScalarField phiG("phiG",(rho * (Mf & g)) & mesh.Sf());


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
    fvc::reconstruct(phiG)
);

#include "createPhi.H"

// - Dimensional unit to use fvOptions for the pressure
dimensionedScalar dims(Mf.dimensions()*dimTime/dimArea,1.0);

const dimensionedScalar tau0(transportProperties.lookup("tau"));

// - fracture/matrix transfer rate
volScalarField tau
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    tau0
);

Info<< "Reading field p_fr\n" << endl;
volScalarField p_fr
(
    IOobject
    (
        "p_fr",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

mesh.setFluxRequired(p_fr.name());

Info << "Reading field K_fr\n" << endl;
const volSymmTensorField K_fr
(
    IOobject
    (
        "K_fr",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
      "K_fr",
      dimArea,
      transportProperties.lookupOrDefault<scalar>("K_fr",scalar(1.))
      *
      symmTensor::I
    )
);

const surfaceSymmTensorField Mf_fr
(
    IOobject
    (
        "Mf_fr",
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
            ).interpolate(inv(K_fr/mu))
    )  // harmonic average of permeabilities
);

surfaceScalarField phi_fr
(
    IOobject
    (
        "phi_fr",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (
        (
            Mf_fr
            & (
                fvc::interpolate
                (
                    fvc::grad(p_fr)
                )
            )
        )
    )
    & mesh.Sf()
);

volVectorField U_fr
(
    IOobject
    (
        "U_fr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::reconstruct(phi_fr)
);


#include "createFvOptions.H"
