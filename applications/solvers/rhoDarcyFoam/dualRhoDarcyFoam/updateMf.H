//- Update viscosity
#include "updateMu.H"

//- Update mobilities
if (updateMf)
{
  Mf = inv(reverseLinear<SymmTensor<scalar>>(mesh).interpolate(inv(K/mu)));
  MfIso = tr(Mf)/scalar(3);
  MfAniso = dev(Mf);
  Mf_fr = inv(reverseLinear<SymmTensor<scalar>>(mesh).interpolate(inv(K_fr/mu_fr)));
  MfIso_fr = tr(Mf_fr)/scalar(3);
  MfAniso_fr = dev(Mf_fr);
}
