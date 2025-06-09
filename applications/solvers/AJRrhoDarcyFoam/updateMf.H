//- Update viscosity
#include "updateMu.H"

//- Update mobilities
if (updateMf)
{
  // standard interpolation
  // Mf = fvc::interpolate(K/mu);
  // two-point flux approximation based on harmonic
  Mf = inv(reverseLinear<SymmTensor<scalar>>(mesh).interpolate(inv(K/mu)));
  MfIso = tr(Mf)/scalar(3);
  MfAniso = dev(Mf);
}
