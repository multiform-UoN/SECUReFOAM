// Set random seed
clock_t time_= std::clock();
srand (time_);

label rndSeed(rand());
Info << "Seed is : " << rndSeed << endl;

Random rng(rndSeed+Pstream::myProcNo());

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

const word grftype = setRandomFieldDict.lookupOrDefault<word>("type","logNormal");
const word corrtype = setRandomFieldDict.lookupOrDefault<word>("correlation","exponential");
const word varName = setRandomFieldDict.lookupOrDefault<word>("field","K");
Info << "Creating a " << grftype <<
  " random field \n  with "<<  corrtype <<
  " correlation for field " << varName <<
  endl << endl;


IOobject fieldHeader
(
  varName,
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ
);


const vector Lcorr
(
  setRandomFieldDict.lookupOrDefault<vector>("Lcorr",vector(1,1,1))
);
Info << "Lcorr " << Lcorr << endl;

vector maxLength
(
  setRandomFieldDict.lookupOrDefault<vector>("maxLength",vector(10,10,10))
);
Info << "maxLength " << maxLength << endl;

const scalar corrTol
(
  setRandomFieldDict.lookupOrDefault<scalar>("corrTol",scalar(1e-8))
);
Info << "corrTol " << corrTol << endl;

const scalar nu
(
  setRandomFieldDict.lookupOrDefault<scalar>("nu",scalar(1))
);
Info << "nu " << nu << endl;

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

const label nfreq
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

const bool printStats
(
  setRandomFieldDict.lookupOrDefault<bool>("printStats",true)
);

const bool setFinalStats
(
  setRandomFieldDict.lookupOrDefault<bool>("setFinalStats",true)
);
Info << "setFinalStats " << setFinalStats << endl;

const bool printLogStats
(
  setRandomFieldDict.lookupOrDefault<bool>("printLogStats",true)
);

const bool periodic
(
  setRandomFieldDict.lookupOrDefault<bool>("periodic",true)
);
Info << "periodic " << periodic << endl;

// dimensions (2 or 3)
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

//const tensorField Keye(mesh.nCells(), tensor::I);
const scalar pi(M_PI);

// correlation lenghts
const scalar dx=Lcorr.x();
const scalar dy=Lcorr.y();
const scalar dz=Lcorr.z();

// Set frequency vectors
nfreqy = nfreqy*(!disableY) + (disableY);
nfreqz = nfreqz*(!disableZ) + (disableZ);
scalarField fx(nfreq,scalar(1));
scalarField fy(nfreqy,scalar(1));
scalarField fz(nfreqz,scalar(1));

const vector meshSpan(mesh.bounds().span());

maxLength = cmptMultiply(maxLength,Lcorr*(!periodic));

// Maximum wave length (if simmetric it can't be more than domain)
const scalar maxlx( Foam::max( maxLength.x(), meshSpan.x() ) );
const scalar maxly( Foam::max( maxLength.y(), meshSpan.y() ) );
const scalar maxlz( Foam::max( maxLength.z(), meshSpan.z() ) );

forAll(fx,ai)
{
  fx[ai] = scalar(ai)/maxlx*dx;
}
forAll(fy,ai)
{
  fy[ai] = scalar(ai)/maxly*dy*(!disableY);
}
forAll(fz,ai)
{
  fz[ai] = scalar(ai)/maxlz*dz*(!disableZ);
}

Info << "  WARNING: The highest frequency solved compared to the correlation length is "
     << endl << "  " <<  fx[nfreq-1] << " " << fy[nfreqy-1] << " " << fz[nfreqz-1]
     << "  This should be more at least 1 or 2!" << endl
     << "  Check also that the mesh spacing is smaller than the correlation length!"
     << endl << endl;

scalar ff(scalar(0));
scalar aa(scalar(0));
scalarField K(mesh.nCells(), scalar(0));
scalarField K1(mesh.nCells(), scalar(0));
scalarField K2(mesh.nCells(), scalar(0));
scalarField* ptrK(&K);

scalar vol(gSum(mesh.V()));
