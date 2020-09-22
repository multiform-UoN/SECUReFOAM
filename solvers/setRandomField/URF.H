dimensionedScalar mu
(
  setRandomFieldDict.lookup("Kmean")
);

dimensionedScalar sigma
(
  setRandomFieldDict.lookup("Ksigma")
);

//- Totally random field
forAll(K,cellI)
{
  // K[cellI] = (scalar(rand())/scalar(RAND_MAX));
}

#include "printStats.H"

K = (K*sigma.value()) + mu.value();
