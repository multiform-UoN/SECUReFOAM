while(pimple.correct())
{
  while(pimple.correctNonOrthogonal())
  {
    #include "updatePEqn.H"

    fvOptions.constrain(pEqn);
    pEqn.relax();
    pEqn.setReference(pRefCell, pRefValue);

    //- Solve for pp with A-H splitting
    const volScalarField InvA("InvA",scalar(1)/p_frEqn.A());
    const volScalarField HbyA("HbyA",p_frEqn.H()*InvA);
    solve(pEqn==tau/rho*(HbyA + fvm::Sp(tau/rho_fr*InvA,pp)));

    // //- Solve for pp uncoupled
    // solve(pEqn==tau/rho*pp_fr);

    fvOptions.correct(pp);

    if (pimple.finalNonOrthogonalIter())
    {
      phi = pEqn.flux() - phiG;
    }

  }

  while(pimple.correctNonOrthogonal())
  {
    #include "updatePEqn.H"

    fvOptions.constrain(p_frEqn);
    p_frEqn.relax();
    
    //- Solve for pp_fr with A-H splitting
    const volScalarField InvA("InvA",scalar(1)/pEqn.A());
    const volScalarField HbyA("HbyA",pEqn.H()*InvA);
    solve(p_frEqn==tau/rho_fr*(HbyA + fvm::Sp(tau/rho*InvA,pp_fr)));

    // //- Solve for pp_fr uncoupled
    // solve(p_frEqn==tau/rho_fr*pp);

    fvOptions.correct(pp_fr);

    if (pimple.finalNonOrthogonalIter())
    {
      phi_fr = p_frEqn.flux() - phiG_fr;
    }

  }

}
