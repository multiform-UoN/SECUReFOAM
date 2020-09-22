#include "GRF.H"
K1 = K;

// Reading thresholds
scalarField TT0(2,mu);
TT0[0] -= sigma;
TT0[1] += sigma;
const scalarField TT
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("thresholds",TT0)
);

// Reading value
scalarField KK0(3,mu);
KK0[1] *= 10;
KK0[2] /= 10;
const scalarField KK
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("values",KK0)
);

if (KK.size()!=TT.size()+1)
{
  FatalError
  << "truncated GRF needs number of values" << nl
  << "equal to the number of thresholds + 1 " << nl
  << exit(FatalError);
}

// Assume we're in the last case
K = KK[TT.size()];

forAll(K, i) // loop through cells
  {
    forAll(TT,j)  // loop through thresholds
    {
      if (K1[i]<TT[j]) // we're in the j-th case
      {
        K[i] = KK[j];
        break;
      }
    }
  }
