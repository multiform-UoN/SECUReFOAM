//while(pimple.correctNonOrthogonal())
//  {
    fvScalarMatrix cCEqn
    (
        fvm::ddt(1., cC)
      + fvm::div(phi, cC)
      - fvm::laplacian(Dc, cC, "laplacian(Dc, cC)")
      - Dam*cA*cB

        ==
           fvOptions(cC)
    );

    fvOptions.constrain(cCEqn);
    cCEqn.relax();
    cCEqn.solve();
    fvOptions.correct(cC);
//  }

