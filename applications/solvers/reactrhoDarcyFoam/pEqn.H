while(pimple.correctNonOrthogonal())
{
  fvScalarMatrix pEqn
  (
    -  fvm::laplacian(Mf,p)
    -  fvc::div(phiG)
   ==
    dims*fvOptions(p)
  );


  fvOptions.constrain(pEqn);
  pEqn.relax();
  pEqn.setReference(pRefCell, pRefValue);
  pEqn.solve();
  fvOptions.correct(p);

  if (pimple.finalNonOrthogonalIter())
  {
    phi = pEqn.flux() - phiG;
  }

}

