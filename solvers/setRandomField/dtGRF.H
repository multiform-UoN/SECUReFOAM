#include "GRF.H"
K2 = K1 - K2;
K1 = K;

// Read threasholds
scalarField TT01(1,mu);
scalarField TT02(1,mu);
const scalarField TT
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("thresholds",TT01)
);
const scalarField TT2
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("thresholds2",TT02)
);


// Read values as a vector
scalarField KK0(4,mu);
KK0[1] *= 10;
KK0[2] /= 10;
const scalarField KK
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("values",KK0)
);


if (KK.size()!=(TT.size()+1)*(TT2.size()+1))
{
  FatalError
  << "bitruncated GRF needs number of values" << nl
  << "equal to  2*(thresholds + 1) " << nl
  << exit(FatalError);
}

// // rescale to match mean and variance (theoretical one)
// const scalar renorm
// (
//   Foam::sqrt
//   (
//     Foam::pow(scalar(2),dim)
//     /
//     (n*aasum)
//   )
// );
// K1 *= renorm;
// K2 *= renorm;


// rescale to match mean and variance (brute-force)
{
  const scalar meank(gSum(K1*mesh.V())/vol); // spatial average
  const scalar vark(gSum(pow(K1,2)*mesh.V())/vol); // spatial variance
  K1 = (K1-meank)/Foam::sqrt((-meank*meank+vark));
}
{
  const scalar meank(gSum(K2*mesh.V())/vol); // spatial average
  const scalar vark(gSum(pow(K2,2)*mesh.V())/vol); // spatial variance
  K2 = (K2-meank)/Foam::sqrt((-meank*meank+vark));
}

K1 = (K1*sigma + mu);
K2 = (K2*sigma + mu);

// Creating thresholds
// we start by assuming is the last column and last row
K = KK[KK.size()-1];

scalarField V(KK.size(), scalar(0)); // facies volumes

forAll(K, i) // loop over cells
{
  forAll(TT,j) // loop over rows
  {
    if (K1[i]<TT[j]) // it's the j-th row
    {
      // assume it's the last column
      K[i] = KK[TT2.size()*(TT.size()+1)+j];
      V[TT2.size()*(TT.size()+1)+j] += mesh.V()[i];
      forAll(TT2,k) // loop through columns
      {
        if (K2[i]<TT2[k]) // it's the k-th column
        {
          K[i] = KK[k*(TT.size()+1)+j];
          V[k*(TT.size()+1)+j] += mesh.V()[i];
          V[TT2.size()*(TT.size()+1)+j] -= mesh.V()[i];
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
    V[(TT2.size()+1)*(TT.size()+1)-1] += mesh.V()[i];
    forAll(TT2,k) // loop through columns
    {
      if (K2[i]<TT2[k]) // it's the k-th column
      {
        K[i] = KK[(k+1)*(TT.size()+1)-1];
        V[(k+1)*(TT.size()+1)-1] += mesh.V()[i];
        V[(TT2.size()+1)*(TT.size()+1)-1] -= mesh.V()[i];
        break;
      }
    }
  }
}

Info << "field values: " << KK
<< endl << "volumes: " << V/vol
<< endl << "thresholds on GRF1: " << TT
<< endl << "thresholds on GRF2: " << TT2 << endl;
