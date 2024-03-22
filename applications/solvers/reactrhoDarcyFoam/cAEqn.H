//while(pimple.correctNonOrthogonal())
//  {
    fvScalarMatrix cAEqn
    (
        fvm::ddt(1., cA)
	+ fvm::div(phi, cA)
	- fvm::laplacian(Da, cA,"laplacian(Da, cA)")
	+ Dam*cA*cB
        ==
           fvOptions(cA)
    );

    fvOptions.constrain(cAEqn);
    cAEqn.relax();
    cAEqn.solve();
    fvOptions.correct(cA);
    //  }


