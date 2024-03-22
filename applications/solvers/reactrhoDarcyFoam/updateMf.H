//- Update viscosity
#include "updateMu.H"

//- Update mobilities
if (updateMf)
{
  Mf = inv(reverseLinear<SymmTensor<scalar>>(mesh).interpolate(inv(K/mu)));
}
