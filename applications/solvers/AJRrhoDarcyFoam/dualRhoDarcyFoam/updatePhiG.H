//- Update density
#include "updateRho.H"

//- Update density-driven flux
if (updatePhiG)
{
  if (hydrostaticPressure)
  {
    // - Spherical (isotropic) part
    phiG    =  MfIso * fvc::snGrad(rho) * ghf * mesh.magSf();
    phiG_fr =  MfIso_fr * fvc::snGrad(rho_fr) * ghf * mesh.magSf();

    // - deviatoric (anisotropic) part
    phiG    +=
              ghf
              *
              (
                MfAniso
                &
                fvc::interpolate(fvc::grad(rho))
              )
              &
              mesh.Sf();
    phiG_fr +=
              ghf
              *
              (
                MfAniso_fr
                &
                fvc::interpolate(fvc::grad(rho_fr))
              )
              &
              mesh.Sf();
  }
  else
  {
    phiG = ( -fvc::interpolate(rho) * (Mf & g) ) & mesh.Sf();
    phiG_fr = ( -fvc::interpolate(rho_fr) * (Mf_fr & g) ) & mesh.Sf();
  }
}
