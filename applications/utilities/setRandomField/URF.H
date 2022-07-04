//- Totally random field
forAll(K,cellI)
{
   K[cellI] = (scalar(rand())/scalar(RAND_MAX));
}

#include "printStats.H"

// rescale to match mean and variance (brute-force)
{
  const scalar meank(gSum(K*V)/vol); // spatial average
  const scalar vark(gSum(pow(K,2)*V)/vol); // spatial variance
  K = (K-meank)/Foam::sqrt((-meank*meank+vark));
  K = (sigma*K+mu);
}
