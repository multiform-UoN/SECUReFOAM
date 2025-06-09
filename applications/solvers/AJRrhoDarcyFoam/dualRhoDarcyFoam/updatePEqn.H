//- TODO add weak compressibility

fvScalarMatrix p_frEqn // same as below but without tau*pp
(
  -  fvm::laplacian(Mf_fr,pp_fr)
  -  fvc::div(phiG_fr)
  ==
      fvm::Sp(-tau/rho_fr,pp_fr)
  +  (dims/rho_fr)*fvOptions(pp_fr)
  +  (hydrostaticPressure)*tau/rho_fr*gh*(rho-rho_fr)
);

fvScalarMatrix pEqn // same as below but without tau*pp_fr
(
  -  fvm::laplacian(Mf,pp)
  -  fvc::div(phiG)
  ==
      fvm::Sp(-tau/rho,pp)
  +  (dims/rho)*fvOptions(pp)
  +  (hydrostaticPressure)*tau/rho*gh*(rho_fr-rho)
);