#include "../createFields.H"

Info << endl << "Reading field c_fr" << endl;
volScalarField c_fr
(
    IOobject
    (
        "c_fr",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createRhoMuModels.H"

Info << endl << "Reading field p_fr" << endl;
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

Info << endl << "Reading field p_rgh_fr" << endl;
volScalarField p_rgh_fr
(
    IOobject
    (
        "p_rgh_fr",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
// p_rgh = p - hydrostaticPressure*rho*gh;

mesh.setFluxRequired(p_rgh_fr.name());
mesh.setFluxRequired(p_fr.name());

Info << endl << "Reading porosity beta_fr" << endl;
volScalarField beta_fr
(
    IOobject
    (
        "beta_fr",
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
        "beta_fr", dimensionedScalar(dimless,scalar(1.))
      )
    )
);


Info << endl << "Reading storativity Ss_fr" << endl;
volScalarField Ss_fr
(
 IOobject
 (
  "Ss_fr",
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
   "Storativity_fr", dimensionedScalar(dimless/dimPressure,scalar(0.))
  )
 )
);

//- Read permeability
Info << endl << "Reading permeability field K_fr" << endl;
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

surfaceSymmTensorField Mf_fr
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
            ).interpolate(inv(K_fr/mu_fr))
    )
);

surfaceScalarField MfIso_fr("MfIso_fr",tr(Mf)/scalar(3));
surfaceSymmTensorField MfAniso_fr("MfAniso_fr",dev(Mf));

Info << endl << "Calculating face flux field phi_fr" << endl;
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
                    fvc::grad(p_rgh_fr)
                )
            )
        )
    )
   & mesh.Sf()
);

//- Create density-driven flux
surfaceScalarField phiG_fr
(
    "phiG_fr",
    phi_fr
);

Info << endl << "Calculating field U_fr" << endl;
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

//- Diffusion tensoor
volSymmTensorField  D_fr
(
    IOobject
    (
        "D_fr",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    Dm * symmTensor::I
);

//- Transversal dispersion length scale
volScalarField alphaT_fr
(
    IOobject
    (
        "alphaT_fr",
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
        "alphaT_fr", dimensionedScalar(dimLength, scalar(0))
      )
    )
);

//- Longitudinal dispersion length scale
volScalarField alphaL_fr
(
    IOobject
    (
        "alphaL_fr",
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
        "alphaL_fr", dimensionedScalar(dimLength, scalar(0))
      )
    )
);

surfaceSymmTensorField Df_fr
(
    IOobject
    (
        "Df_fr",
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
            ).interpolate(inv(D_fr))
    )  // harmonic average of diffusivity
);

#include "updateD.H"

Info << endl << "Mean dispersion tensor is "
     << (fvc::domainIntegrate(D_fr)/gSum(mesh.V())).value() << endl;

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


//- Store old mass to compute mass rate
scalar oldMassTot
(
	gSum((rho.primitiveField() + rho_fr.primitiveField())*mesh.V())
);

scalar oldMassFracture
(
	gSum(rho_fr.primitiveField()*mesh.V())
);
