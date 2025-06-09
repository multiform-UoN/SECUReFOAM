//- update concentration in the immobile regions
multiRate->correct();

while(pimple.correctNonOrthogonal())
{
    fvScalarMatrix cEqn
    (
        fvm::ddt(beta,c)
      + multiRate->source()
      + fvm::div(phi, c)
      - fvm::laplacian(Df,c,"laplacian(D,c)")
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
