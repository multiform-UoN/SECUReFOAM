IOdictionary setRandomFieldDict
(
  IOobject
  (
    "setRandomFieldDict",
    runTime.system(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
  )
);

// Set random seed
clock_t time_= std::clock();
srand (time_);
label seed(0);

if(setRandomFieldDict.found("seed"))
{
    Info << "Reading seed from dictionary" << endl;
    seed = readLabel(setRandomFieldDict.lookup("seed"));
}
else
{
    seed = rand();
}

Info << "Seed is : " << seed << endl;

Random rng(seed);
const scalar pi(M_PI);


// - Reading the field
const word varName = setRandomFieldDict.lookupOrDefault<word>("field","K");

IOobject fieldHeader
(
  varName,
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ
);

// -  Random field type
const word grftype = setRandomFieldDict.lookupOrDefault<word>("type","logNormal");

// -  Correlation type
const word corrtype = setRandomFieldDict.lookupOrDefault<word>("correlation","exponential");

const scalar nu // for Matern only
(
  setRandomFieldDict.lookupOrDefault<scalar>("nu",scalar(1))
);
Info << "nu " << nu << endl;

Info << endl << "Creating a " << grftype <<
  " random field \n  with "<<  corrtype <<
  " correlation for field " << varName <<
  endl;

// -  Correlation length, tolerance (to skip frequencies)
const vector Lcorr
(
  setRandomFieldDict.lookupOrDefault<vector>("Lcorr",vector(1,1,1))
);
Info << "Lcorr " << Lcorr << endl;

const scalar corrTol
(
  setRandomFieldDict.lookupOrDefault<scalar>("corrTol",scalar(1e-8))
);
Info << "corrTol " << corrTol << endl;

// -  Mean and variance and option to force it at the end
const scalar mu
(
  setRandomFieldDict.lookupOrDefault<scalar>("mu",scalar(0))
);
Info << "mu " << mu << endl;

const scalar sigma
(
  setRandomFieldDict.lookupOrDefault<scalar>("sigma",scalar(1))
);
Info << "sigma " << sigma << endl;

const bool setFinalStats
(
  setRandomFieldDict.lookupOrDefault<bool>("setFinalStats",false)
);
Info << "setFinalStats " << setFinalStats << endl;

// - Number of frequencies in all directions
label nfreq
(
  setRandomFieldDict.lookupOrDefault<label>("nfreq",20)
);
label nfreqy
(
  setRandomFieldDict.lookupOrDefault<label>("nfreqy",nfreq)
);
label nfreqz
(
  setRandomFieldDict.lookupOrDefault<label>("nfreqz",nfreq)
);
Info << "nfreq " << nfreq << " " << nfreqy << " " << nfreqz << endl;

// - Skip frequencies to speed up for large domains and small Lcorr
label skipfreq
(
  setRandomFieldDict.lookupOrDefault<label>("skipfreq",0)
);
label skipfreqy
(
  setRandomFieldDict.lookupOrDefault<label>("skipfreqy",skipfreq)
);
label skipfreqz
(
  setRandomFieldDict.lookupOrDefault<label>("skipfreqz",skipfreq)
);
Info << "skipfreq " << skipfreq << " " << skipfreqy << " " << skipfreqz << endl;

// - Overwrite the initial field or use it to add/multiply
const bool reset
(
  setRandomFieldDict.lookupOrDefault<bool>("reset",true)
);
Info << "reset " << reset << endl;

const bool additive
(
  setRandomFieldDict.lookupOrDefault<bool>("additive",false)
);
Info << "additive " << additive << endl;

// - Print stats
const bool printStats
(
  setRandomFieldDict.lookupOrDefault<bool>("printStats",true)
);

const bool printLogStats
(
  setRandomFieldDict.lookupOrDefault<bool>("printLogStats",true)
);

// Set dimensions (X cannot be removed)
const bool disableY
(
  setRandomFieldDict.lookupOrDefault<bool>("disableY",false)
);
Info << "disableY " << disableY << endl;

const bool disableZ
(
  setRandomFieldDict.lookupOrDefault<bool>("disableZ",false)
);
Info << "disableZ " << disableZ << endl;

const label dim(3 - disableY - disableZ );
Info << "Number of dimensions: " << dim << endl;

// - Domain size (cell-based)
vector meshSpan
      (
        Foam::gMax( (mesh.C().component(0))() ) - Foam::gMin( (mesh.C().component(0))() ),
        Foam::gMax( (mesh.C().component(1))() ) - Foam::gMin( (mesh.C().component(1))() )*(!disableY)+disableY,
        Foam::gMax( (mesh.C().component(2))() ) - Foam::gMin( (mesh.C().component(2))() )*(!disableZ)+disableZ
      );

// - Periodicity (for periodic cells)
const bool periodic
(
  setRandomFieldDict.lookupOrDefault<bool>("periodic",true)
);
Info << "periodic " << periodic << endl;

// - for periodic faces
const bool periodicFaces
(
  setRandomFieldDict.lookupOrDefault<bool>("facePeriodicity",false)
);
Info << "facePeriodicity " << periodicFaces << endl;

// - force period equal to domain size
const bool fullLengthPeriodicity
(
  setRandomFieldDict.lookupOrDefault<bool>("fullPeriodicity",true)
);
Info << "fullPeriodicity " << fullLengthPeriodicity << endl;

// - Domain size (face-based)
if (periodicFaces)
{
  meshSpan = mesh.bounds().span();
}
Info << "meshSpan " << meshSpan << endl;

// Maximum wave lengths wrt correlation
scalar maxLengthByL
(
  setRandomFieldDict.lookupOrDefault<scalar>
    (
      "maxLengthByL",
      scalar(4)
    )
);

// Minimum wave lengths (sets the number of frequency used)
scalar minLengthByL
(
  setRandomFieldDict.lookupOrDefault<scalar>
    (
      "minLengthByL",
      scalar(1./4)
    )
);

// Minimum wave lengths (sets the number of frequency used)
scalar minLengthByCellSize
(
  setRandomFieldDict.lookupOrDefault<scalar>
    (
      "minLengthByCellSize",
      scalar(4)
    )
);

// Maximum wave lengths (overwrites the maxLengthByL if found)
vector maxLength
(
  setRandomFieldDict.lookupOrDefault<vector>
    (
      "maxLength",
      Lcorr*maxLengthByL
    )
);

if (fullLengthPeriodicity)
{
  maxLength = meshSpan;
}

// - set maximum wave length such that is a divisor of domain size
if (periodic)
{
  maxLength[0] = meshSpan[0]/ceil(meshSpan[0]/(maxLength[0]+SMALL));
  maxLength[1] = meshSpan[1]/ceil(meshSpan[1]/(maxLength[1]+SMALL));
  maxLength[2] = meshSpan[2]/ceil(meshSpan[2]/(maxLength[2]+SMALL));
}
Info << "maxLength " << maxLength << endl;

// - correlation lenghts
const scalar dx=Lcorr.x();
const scalar dy=Lcorr.y();
const scalar dz=Lcorr.z();

// - Set frequency vectors
nfreqy = nfreqy*(!disableY) + (disableY);
nfreqz = nfreqz*(!disableZ) + (disableZ);

#include "computeCellSize.H"
const vector minLength(max(Lcorr*minLengthByL,cellSize*minLengthByCellSize));
Info << "minLength " << minLength << endl;

// - Check that enough frequencies are present
if (nfreq/maxLength.x()*minLength.x()<scalar(1))
{
  nfreq = ceil(maxLength.x()/(minLength.x()));
  Info << "  WARNING: Not enough frequencies in x "
       << endl << "       Changed to " <<  nfreq << endl;
}
if ((!disableY) && nfreqy/maxLength.y()*minLength.y()<scalar(1))
{
  nfreqy = ceil(maxLength.y()/(minLength.y()));
  Info << "  WARNING: Not enough frequencies in y "
       << endl << "       Changed to " <<  nfreqy << endl;
}
if ((!disableZ) && nfreqz/maxLength.z()*minLength.z()<scalar(1))
{
  nfreqz = ceil(maxLength.z()/(minLength.z()));
  Info << "  WARNING: Not enough frequencies in z "
       << endl << "       Changed to " <<  nfreqz << endl;
}

// - Compute frequencies
scalarField fx(nfreq,scalar(1));
scalarField fy(nfreqy,scalar(1));
scalarField fz(nfreqz,scalar(1));

forAll(fx,ai)
{
  fx[ai] = scalar(ai)/maxLength.x()*dx;
}
forAll(fy,ai)
{
  fy[ai] = scalar(ai)/maxLength.y()*dy*(!disableY);
}
forAll(fz,ai)
{
  fz[ai] = scalar(ai)/maxLength.z()*dz*(!disableZ);
}

Info << endl << "Generating and storing random numbers" << endl;
scalarField rngField(nfreq*nfreqy*nfreqz*8,scalar(0));
if (Pstream::master())
{
    forAll(rngField,rngI)
    {
      rngField[rngI] = rng.scalarNormal();
    }
}
Pstream::scatter(rngField);
