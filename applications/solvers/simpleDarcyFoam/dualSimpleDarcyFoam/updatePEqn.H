fvScalarMatrix pEqn
(
  -  fvm::laplacian(Mf,p)
  ==
      fvm::Sp(-tau,p)
  +  dims*fvOptions(p)
);

fvScalarMatrix p_frEqn
(
  -  fvm::laplacian(Mf_fr,p_fr)
  ==
      fvm::Sp(-tau,p_fr)
  +  dims*fvOptions(p_fr)
);
