int iCorr = 0;
scalar initialResidual = 0;

do
{
	//- Create stress flux
	surfaceVectorField phiStress
	(
		fvc::interpolate
		(
			(mu*T(gradD) + lambda*tr(gradD)*I),
			"interpolate(stress)"
		)&mesh.Sf()
	);

	while(pimple.correctNonOrthogonal())
	{
		fvVectorMatrix DEqn
		(
			fvm::d2dt2(rhos,D)
		  - fvm::laplacian(mu,D)
		  ==
			fvc::div(phiStress)
		  + rhos*g
		  - (alpha*fvc::grad(p))
		);

		DEqn.relax();
		initialResidual = DEqn.solve().max().initialResidual();
	}

	gradD = fvc::grad(D);
	J = scalar(1.0) + fvc::div(D);


} while (initialResidual > convergenceTolerance && ++iCorr < nCorr);
