scalar aasum(0.0);

forAll(fx,kx)
{
  forAll(fy,ky)
  {
    forAll(fz,kz)
    {
      #include "corrLength.H"

      if (aa<corrTol)
      {
        break;
      }

      K1 += aa * rng.scalarNormal()*
        cos(scalar(2)*pi*(mesh.C().primitiveField()&vector(fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)));
      K2 += aa * rng.scalarNormal()*
        sin(scalar(2)*pi*(mesh.C().primitiveField()&vector(fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)));

      K1 += aa * rng.scalarNormal()*
        cos(scalar(2)*pi*(mesh.C().primitiveField()&vector(-fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)));
      K2 += aa * rng.scalarNormal()*
        sin(scalar(2)*pi*(mesh.C().primitiveField()&vector(-fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)));

      K1 += aa * rng.scalarNormal()*
        cos(scalar(2)*pi*(mesh.C().primitiveField()&vector(fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)));
      K2 += aa * rng.scalarNormal()*
        sin(scalar(2)*pi*(mesh.C().primitiveField()&vector(fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)));

      K1 += aa * rng.scalarNormal()*
        cos(scalar(2)*pi*(mesh.C().primitiveField()&vector(fx[kx]/dx,fy[ky]/dy,-fz[kz]/dz)));
      K2 += aa * rng.scalarNormal()*
        sin(scalar(2)*pi*(mesh.C().primitiveField()&vector(fx[kx]/dx,fy[ky]/dy,-fz[kz]/dz)));

      K1 += aa * rng.scalarNormal()*
        cos(scalar(2)*pi*(mesh.C().primitiveField()&vector(-fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)));
      K2 += aa * rng.scalarNormal()*
        sin(scalar(2)*pi*(mesh.C().primitiveField()&vector(-fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)));

      K1 += aa * rng.scalarNormal()*
        cos(scalar(2)*pi*(mesh.C().primitiveField()&vector(-fx[kx]/dx,fy[ky]/dy,-fz[kz]/dz)));
      K2 += aa * rng.scalarNormal()*
        sin(scalar(2)*pi*(mesh.C().primitiveField()&vector(-fx[kx]/dx,fy[ky]/dy,-fz[kz]/dz)));

      K1 += aa * rng.scalarNormal()*
        cos(scalar(2)*pi*(mesh.C().primitiveField()&vector(-fx[kx]/dx,-fy[ky]/dy,-fz[kz]/dz)));
      K2 += aa * rng.scalarNormal()*
        sin(scalar(2)*pi*(mesh.C().primitiveField()&vector(-fx[kx]/dx,-fy[ky]/dy,-fz[kz]/dz)));

      K1 += aa * rng.scalarNormal()*
        cos(scalar(2)*pi*(mesh.C().primitiveField()&vector(fx[kx]/dx,-fy[ky]/dy,-fz[kz]/dz)));
      K2 += aa * rng.scalarNormal()*
        sin(scalar(2)*pi*(mesh.C().primitiveField()&vector(fx[kx]/dx,-fy[ky]/dy,-fz[kz]/dz)));

    }
  }
}

K = K1 + K2;

// // rescale to match mean and variance (theoretical one)
// const scalar renorm
// ( Foam::sqrt
//   (
//     Foam::pow(scalar(2),dim)
//     /
//     (n*aasum)
//   )
// );
// K *= renorm;

// rescale to match mean and variance (brute-force)
{
  const scalar meank(gSum(K*mesh.V())/vol); // spatial average
  const scalar vark(gSum(pow(K,2)*mesh.V())/vol); // spatial variance
  K = (K-meank)/Foam::sqrt((-meank*meank+vark));
}

K = (sigma*K+mu);
