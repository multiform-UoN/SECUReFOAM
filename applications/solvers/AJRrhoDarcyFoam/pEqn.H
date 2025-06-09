while(pimple.correctNonOrthogonal())
{
  fvScalarMatrix pEqn
  (
    -  fvm::laplacian(Mf,pp)
    -  fvc::div(phiG)
    //==
    //(dims/rho)*fvOptions(pp)
  );

  if (!incompressible)
  {
    pEqn += fvm::ddt(Ss,pp);

    if (standardCompressibility)
    {
      pEqn
      +=
      (drhodc/rho)
      *
      (
        beta*fvc::ddt(c)
        +
        (U&fvc::grad(c))
      );
    }
    else
    {
      fvScalarMatrix cS(fvOptions(c));

      pEqn
      +=
      drhodc/(rho-drhodc*c)
      *
      (
        fvc::div
              (
                fvc::interpolate
                (
                  (D & fvc::grad(c))
                )
                &
                mesh.Sf()
              )
        +
        (cS.A()*c - cS.H()) // sources of c 
      );
    }
    if (hydrostaticPressure)
    {
      pEqn
      +=
      Ss*drhodc*fvc::ddt(c)*gh;
    }
  }


  fvOptions.constrain(pEqn);
  pEqn.relax();
  pEqn.setReference(pRefCell, pRefValue);
  pEqn.solve();
  fvOptions.correct(pp);

  if (pimple.finalNonOrthogonalIter())
  {
    phi = pEqn.flux() - phiG;
  }

}
