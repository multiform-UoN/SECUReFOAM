while (pimple.correctNonOrthogonal())
{

	fvScalarMatrix pEqn
	(
		fvm::ddt(rhof*stor,p) -  fvm::laplacian(Mf,p)
		==
		- fvc::div(phiG)
		- fvc::ddt(alpha,J)
    + dims*fvOptions(p)


	);

	fvOptions.constrain(pEqn);
  pEqn.relax();
  pEqn.solve();
  fvOptions.correct(p);

	if (pimple.finalNonOrthogonalIter())
  {
		U = fvc::reconstruct(pEqn.flux() + phiG);
	}
}
