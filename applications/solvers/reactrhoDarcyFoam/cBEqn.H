//while(pimple.correctNonOrthogonal())
//  {
    fvScalarMatrix cBEqn
    (
        fvm::ddt(1., cB)
        + fvm::div(phi, cB)
        - fvm::laplacian(Db, cB, "laplacian(Db, cB)")
        + Dam*cA*cB
        ==
           fvOptions(cB)
    );

    fvOptions.constrain(cBEqn);
    cBEqn.relax();
    cBEqn.solve();
    fvOptions.correct(cB);
  //  }


