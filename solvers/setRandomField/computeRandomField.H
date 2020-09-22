Info << endl <<  "  Start computing Random Field" << endl;

if (grftype=="gaussian")
{
  #include "GRF.H"
}
else if (grftype=="logNormal")
{
  #include "lGRF.H"
}
else if (grftype=="truncated")
{
  #include "tGRF.H"
}
else if (grftype=="bitruncated")
{
  #include "dtGRF.H"
}
else if (grftype=="uniform")
{
  #include "URF.H"
}
else
{
  FatalError
  << "I cannot recognise type " << grftype << nl
  << "Avaliable values are: logNormal, gaussian, truncated, bitruncated, uniform " << nl
  << exit(FatalError);
}

#include "printStats.H"

#include "printLogStats.H"

Info << "Writing the Random Field" << endl;

//  Final scaling
if (setFinalStats)
{
  const scalar meank(gSum(K*mesh.V())/vol); // spatial average
  const scalar vark(gSum(pow(K,2)*mesh.V())/vol); // spatial variance
  K = (K-meank)/Foam::sqrt((-meank*meank+vark));
  K = (sigma*K+mu);
}

//  Compute and write the field f
if (reset)
{
  forAll(f,icell)
  {
    f[icell] = K[icell]*f0;
  }
}
else if (additive)
{
  forAll(f,icell)
  {
    f[icell] = f[icell] + K[icell]*f0;
  }
}
else
{
  forAll(f,icell)
  {
    f[icell] = K[icell]*f[icell];
  }
}
f.correctBoundaryConditions();
f.write();
