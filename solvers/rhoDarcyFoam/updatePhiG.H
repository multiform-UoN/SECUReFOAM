//- Update density-driven flux
#include "updateRho.H"

if (hydrostaticPressure)
{
  // - Spherical (isotropic) part
  phiG    = - MfIso * fvc::snGrad(rho) * ghf * mesh.magSf();

  // - deviatoric (anisotropic) part
  phiG    +=
            ghf
            *
            (
              - MfAniso
              &
              fvc::interpolate(fvc::grad(rho))
            )
            &
            mesh.Sf();
}
else
{
  phiG = ( fvc::interpolate(rho) * (Mf & g) ) & mesh.Sf();
}
