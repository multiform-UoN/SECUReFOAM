// Calculate percentiles and thresholds
#include "erfinv.H"

// Reading thresholds
const scalarField TT
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("thresholds",TT01)
);
Info << endl << "Thresholds " << TT << endl;

// Reading value
scalarField KK0(TT.size()+1);
forAll(KK0,Ki)
{
  KK0[Ki] = Ki; 
}
const scalarField KK
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("values",KK0)
);
Info << endl << "Values " << KK << endl;

if (KK.size()!=TT.size()+1)
{
  FatalError
  << "truncated GRF needs number of values" << nl
  << "equal to the number of thresholds + 1 " << nl
  << exit(FatalError);
}

#include "GRF.H"
K1 = K;

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
