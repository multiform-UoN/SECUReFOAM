while(pimple.correctNonOrthogonal())
{
  fvScalarMatrix pEqn
  (
    -  fvm::laplacian(Mf,p_rgh)
    +  fvc::div(phiG)
    ==
       dims*fvOptions(p_rgh)
  );

  fvOptions.constrain(pEqn);
  pEqn.relax();
  pEqn.solve();
  fvOptions.correct(p_rgh);

  if (pimple.finalNonOrthogonalIter())
  {
    phi = pEqn.flux();// + phiG;
  }

}

//- Update pressure w.o. hydrostatic contributions
p = p_rgh + hydrostaticPressure*rho*gh;
