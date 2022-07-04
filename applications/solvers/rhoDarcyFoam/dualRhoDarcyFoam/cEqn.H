#include "CourantNo.H"
#include "CourantNoG_fr.H"

const volScalarField dp
(
  Foam::pos(p_fr-p)
  *
  (p_fr-p)
);

const volScalarField dpm
(
  (p_fr-p)
  -
  dp
);

// // const volTensorField Deff
// const dimensionedTensor Deff
// (
//   "Deff",
//   D0 //+ K/mu*gh*(rho2-rho1)
// );
//
// // const volTensorField Deff_fr
// const dimensionedTensor Deff_fr
// (
//   "Deff_fr",
//   D0 //+ K_fr/mu*gh*(rho2-rho1)
// );


while(pimple.correct())
{

  while(pimple.correctNonOrthogonal())
  {

    fvScalarMatrix cEqn
    (
        fvm::ddt(beta,c)
      + fvm::div(phi,c)
      - fvm::laplacian(Df,c,"laplacian(D,c)")
      ==
        tau*dp*c_fr*beta_fr/rho_fr
      + fvm::Sp(tau*dpm*beta/rho,c)
      + fvOptions(c)
    );

    fvOptions.constrain(cEqn);
    cEqn.relax();
    cEqn.solve();
    fvOptions.correct(c);
  }
  // //- Bound c
  // c = (mag(c) + c)/scalar(2.0);
  // c.correctBoundaryConditions();


  while(pimple.correctNonOrthogonal())
  {
    fvScalarMatrix cEqnfr
    (
        fvm::ddt(beta_fr,c_fr)
      + fvm::div(phi_fr,c_fr,"div(phi,c)")
      - fvm::laplacian(Df_fr,c_fr,"laplacian(D,c)")
      ==
      - tau*dpm*c*beta/rho
      + fvm::Sp(-tau*dp*beta_fr/rho_fr,c_fr)
      + fvOptions(c_fr)
  );

    fvOptions.constrain(cEqnfr);
    cEqnfr.relax();
    cEqnfr.solve();
    fvOptions.correct(c_fr);
  }
  // //- Bound c
  // c_fr = (mag(c_fr) + c_fr)/scalar(2.0);
  // c_fr.correctBoundaryConditions();

  // #include "updatePhiG.H"

}
