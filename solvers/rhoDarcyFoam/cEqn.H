  #include "CourantNoG.H"

  while(pimple.correctNonOrthogonal())
  {
    fvScalarMatrix cEqn
    (
        fvm::ddt(beta,c)
      + fvm::div(phi,c)
      + fvm::div(phiG,c)
      - fvm::laplacian(D,c,"laplacian(D,c)")
      ==
        fvOptions(c)
    );

    fvOptions.constrain(cEqn);
    cEqn.relax();
    cEqn.solve();
    fvOptions.correct(c);
  }
  // //- Bound c
  // c = (mag(c) + c)/scalar(2.0);
  // c.correctBoundaryConditions();
