//- Update density
#include "updateRho.H"

//- Update density-driven flux
if (updatePhiG)
{
  if (hydrostaticPressure)
  {
    // - Spherical (isotropic) part
    phiG    =  MfIso * fvc::snGrad(c) * ghf * mesh.magSf();

    // - deviatoric (anisotropic) part
    phiG    +=
              ghf
              *
              (
                MfAniso
                &
                fvc::interpolate(fvc::grad(c))
              )
              &
              mesh.Sf();

    phiG *= fvc::interpolate(drhodc);
  }
  else
  {
    phiG = ( -fvc::interpolate(rho) * (Mf & g) ) & mesh.Sf();
  }
}
