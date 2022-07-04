while(simple.correctNonOrthogonal())
{
  #include "updatePEqn.H"

  fvOptions.constrain(pEqn);
  pEqn.relax();
  pEqn.setReference(pRefCell, pRefValue);

  //- Solve for pp with A-H splitting
  const volScalarField InvA("InvA",scalar(1)/p_frEqn.A());
  const volScalarField HbyA("HbyA",p_frEqn.H()*InvA);
  solve(pEqn==tau*(HbyA + fvm::Sp(tau*InvA,p)));

  // //- Solve for p uncoupled
  // solve(pEqn==tau*p_fr);

  fvOptions.correct(p);


  if (simple.finalNonOrthogonalIter())
  {
    phi = pEqn.flux();
  }

}

while(simple.correctNonOrthogonal())
{
  #include "updatePEqn.H"

  fvOptions.constrain(p_frEqn);
  p_frEqn.relax();
  
  //- Solve for p_fr with A-H splitting
  const volScalarField InvA("InvA",scalar(1)/pEqn.A());
  const volScalarField HbyA("HbyA",pEqn.H()*InvA);
  solve(p_frEqn==tau*(HbyA + fvm::Sp(tau*InvA,p_fr)));

  // //- Solve for p_fr uncoupled
  // solve(p_frEqn==tau*p);

  fvOptions.correct(p_fr);

  if (simple.finalNonOrthogonalIter())
  {
    phi_fr = p_frEqn.flux();
  }

}
