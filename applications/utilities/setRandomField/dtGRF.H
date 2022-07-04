//- Read percentiles and convert them to thresholds
#include "erfinv.H"

//- Read thresholds (overwrite percentiles if found)
const scalarField TT
(
 setRandomFieldDict.lookupOrDefault<scalarField> ("thresholds",TT01)
);
const scalarField TT2
(
 setRandomFieldDict.lookupOrDefault<scalarField> ("thresholds2",TT02)
);
Info << endl << "Thresholds " << TT << "     " << TT2 << endl;

// Read values as a vector
scalarField KK0((TT.size()+1)*(TT2.size()+1));
forAll(KK0,Ki)
{
  KK0[Ki] = Ki;
}
const scalarField KK
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("values",KK0)
);
Info << endl << "Values " << KK << endl;

if (KK.size()!=(TT.size()+1)*(TT2.size()+1))
{
  FatalError
  << "bitruncated GRF needs number of values" << nl
  << "equal to  2*(thresholds + 1) " << nl
  << exit(FatalError);
}

//- create two GRFs
#include "GRF.H"
K2 = K1 - K2;
K1 = K;

// rescale to match mean and variance (brute-force)
{
  const scalar meank(gSum(K1*V)/vol); // spatial average
  const scalar vark(gSum(pow(K1,2)*V)/vol); // spatial variance
  K1 = (K1-meank)/Foam::sqrt((-meank*meank+vark));
  K1 = (K1*sigma + mu);
}
{
  const scalar meank(gSum(K2*V)/vol); // spatial average
  const scalar vark(gSum(pow(K2,2)*V)/vol); // spatial variance
  K2 = (K2-meank)/Foam::sqrt((-meank*meank+vark));
  K2 = (K2*sigma + mu);
}

// Creating thresholds
// we start by assuming is the last column and last row
K = KK[KK.size()-1];

scalarField VF(KK.size(), scalar(0)); // facies volumes

forAll(K, i) // loop over cells
{
  forAll(TT,j) // loop over rows
  {
    if (K1[i]<TT[j]) // it's the j-th row
    {
      // assume it's the last column
      K[i] = KK[TT2.size()*(TT.size()+1)+j];
      VF[TT2.size()*(TT.size()+1)+j] += V[i];
      forAll(TT2,k) // loop through columns
      {
        if (K2[i]<TT2[k]) // it's the k-th column
        {
          K[i] = KK[k*(TT.size()+1)+j];
          VF[k*(TT.size()+1)+j] += V[i];
          VF[TT2.size()*(TT.size()+1)+j] -= V[i];
          break;
        }
      }
      break;
    }
  }
  if (K1[i]>=TT[TT.size()-1]) // it's the last row
  {
    // assume it's the last column
    K[i] = KK[(TT2.size()+1)*(TT.size()+1)-1];
    VF[(TT2.size()+1)*(TT.size()+1)-1] += V[i];
    forAll(TT2,k) // loop through columns
    {
      if (K2[i]<TT2[k]) // it's the k-th column
      {
        K[i] = KK[(k+1)*(TT.size()+1)-1];
        VF[(k+1)*(TT.size()+1)-1] += V[i];
        VF[(TT2.size()+1)*(TT.size()+1)-1] -= V[i];
        break;
      }
    }
  }
}

reduce(VF,sumOp<tmp<scalarField>>());

Info << endl << "truncation rule applied"
<< endl << "field values: " << KK
<< endl << "volumes: " << VF/vol << " " << gSum(VF/vol)
<< endl << "thresholds on GRF1: " << TT
<< endl << "thresholds on GRF2: " << TT2 << endl;
