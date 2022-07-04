//- Allocate main variables
scalar ff(scalar(0));
scalar aa(scalar(0));
scalarField K(f.size(), scalar(0));
scalarField K1(f.size(), scalar(0));
scalarField K2(f.size(), scalar(0));
scalarField* ptrK(&K);

Info << endl <<  "Start computing Random Field" << endl;

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

Info << endl << "Writing the Random Field" << endl;

//  Final scaling
if (setFinalStats)
{
  #include "setStats.H"
}

#include "printStats.H"
#include "printLogStats.H"

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
//f.correctBoundaryConditions();
f.write();
