//- Update density
#include "updateRho.H"
if (updatePhiG)
  {
   phiG = ( -fvc::interpolate(rho) * (Mf & g) ) & mesh.Sf();
  }
 
