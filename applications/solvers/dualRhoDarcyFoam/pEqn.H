while(pimple.correct())
{
  while(pimple.correctNonOrthogonal())
  {
    fvScalarMatrix p_frEqn // same as below but without tau*p_rgh
    (
      -  fvm::laplacian(Mf_fr,p_rgh_fr)
      +  fvc::div(phiG_fr)
      ==
         fvm::Sp(-tau,p_rgh_fr)
      +  dims*fvOptions(p_rgh_fr)
    );

    // A*rgh_fr = H + tau*rgh
    // rgh_fr = HbyA + fvm::Sp(tau*InvA,rgh)

    const volScalarField InvA("InvA",scalar(1)/p_frEqn.A());
    const volScalarField HbyA("HbyA",p_frEqn.H()*InvA);

    fvScalarMatrix pEqn
    (
      -  fvm::laplacian(Mf,p_rgh)
      +  fvc::div(phiG)
      ==
         fvm::Sp(-tau,p_rgh)
      // +  tau*p_rgh_fr // this is replaced by the line below
      +  tau*(HbyA + fvm::Sp(tau*InvA,p_rgh))
      +  dims*fvOptions(p_rgh)
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

  while(pimple.correctNonOrthogonal())
  {
    fvScalarMatrix pEqn // same as below but without tau*p_rgh
    (
      -  fvm::laplacian(Mf,p_rgh)
      +  fvc::div(phiG)
      ==
         fvm::Sp(-tau,p_rgh)
      +  dims*fvOptions(p_rgh)
    );

    // A*rgh = H + tau*rgh_fr
    // rgh = HbyA + fvm::Sp(tau*InvA,rgh_fr)

    const volScalarField InvA("InvA",scalar(1)/pEqn.A());
    const volScalarField HbyA("HbyA",pEqn.H()*InvA);

    fvScalarMatrix p_frEqn
    (
      -  fvm::laplacian(Mf_fr,p_rgh_fr)
      +  fvc::div(phiG_fr)
      ==
         fvm::Sp(-tau,p_rgh_fr)
      // +  tau*p_rgh
      +  tau*(HbyA + fvm::Sp(tau*InvA,p_rgh_fr))
      +  dims*fvOptions(p_rgh_fr)
    );

    fvOptions.constrain(p_frEqn);
    p_frEqn.relax();
    p_frEqn.solve();
    fvOptions.correct(p_rgh_fr);

    if (pimple.finalNonOrthogonalIter())
    {
      phi_fr = p_frEqn.flux();// + phiG_fr;
    }

  }

}

//- Update pressure w.o. hydrostatic contributions
p = p_rgh + hydrostaticPressure*rho*gh;
p_fr = p_rgh_fr + hydrostaticPressure*rho*gh;
