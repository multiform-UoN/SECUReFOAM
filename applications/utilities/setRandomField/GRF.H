// - keeps track of the total energy (no longer needed since we force the mean/variance)
// scalar aasum(0.0);

label rngI(0);

forAll(fx,kx)
{
  forAll(fy,ky)
  {
    forAll(fz,kz)
    {
      #include "corrLength.H"
      Info << "Frequency " << vector(fx[kx]/dx,fy[ky]/dy,fz[kz]/dz) << " weight " << aa << endl;

      if (aa<corrTol)
      {
        Info << "Skip frequency " << endl;
        break;
      }

      // - keeps track of the total energy (no longer needed since we force the mean/variance)
      // aasum += aa;
      
      // // - Sum using fields
      // {

      //   K1 += aa *
      //     (
      //       rngField[rngI+0]*Foam::cos((
      //       scalar(2)*pi*
      //       (
      //         C
      //         &
      //         vector(fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)
      //       )
      //     ))
      //       +
      //       rngField[rngI+1]*Foam::cos((
      //       scalar(2)*pi*
      //       (
      //         C
      //         &
      //         vector(-fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)
      //       )
      //     ))
      //       +
      //       rngField[rngI+2]*Foam::cos((
      //       scalar(2)*pi*
      //       (
      //         C
      //         &
      //         vector(fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)
      //       )
      //     ))
      //       +
      //       rngField[rngI+3]*Foam::cos((
      //       scalar(2)*pi*
      //       (
      //         C
      //         &
      //         vector(-fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)
      //       )
      //     ))
      //     );
      //   K2 += aa *
      //     (
      //       rngField[rngI+4]*Foam::sin((
      //       scalar(2)*pi*
      //       (
      //         C
      //         &
      //         vector(fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)
      //       )
      //     ))
      //       +
      //       rngField[rngI+5]*Foam::sin((
      //       scalar(2)*pi*
      //       (
      //         C
      //         &
      //         vector(-fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)
      //       )
      //     ))
      //       +
      //       rngField[rngI+6]*Foam::sin((
      //       scalar(2)*pi*
      //       (
      //         C
      //         &
      //         vector(fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)
      //       )
      //     ))
      //       +
      //       rngField[rngI+7]*Foam::sin((
      //       scalar(2)*pi*
      //       (
      //         C
      //         &
      //         vector(-fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)
      //       )
      //     ))
      //     );
      // }

      //- Loop over cells
      forAll(K1,celli)
      {
        const scalar f1
          (
            scalar(2)*pi*
            (
              C[celli]
              &
              vector(fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)
            )
          );
        const scalar f2
          (
            scalar(2)*pi*
            (
              C[celli]
              &
              vector(-fx[kx]/dx,fy[ky]/dy,fz[kz]/dz)
            )
          );
        const scalar f3
          (
            scalar(2)*pi*
            (
              C[celli]
              &
              vector(fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)
            )
          );
        const scalar f4
          (
            scalar(2)*pi*
            (
              C[celli]
              &
              vector(-fx[kx]/dx,-fy[ky]/dy,fz[kz]/dz)
            )
          );

        K1[celli] += aa *
          (
            rngField[rngI+0]*Foam::cos(f1)
            +
            rngField[rngI+1]*Foam::cos(f2)
            +
            rngField[rngI+2]*Foam::cos(f3)
            +
            rngField[rngI+3]*Foam::cos(f4)
          );
        K2[celli] += aa *
          (
            rngField[rngI+4]*Foam::sin(f1)
            +
            rngField[rngI+5]*Foam::sin(f2)
            +
            rngField[rngI+6]*Foam::sin(f3)
            +
            rngField[rngI+7]*Foam::sin(f4)
          );
      }

      rngI += 8;

      kz += skipfreqz;
    }
    ky += skipfreqy;
  }
  kx += skipfreq;
}

K = K1 + K2;

// rescale to match mean and variance (brute-force)
#include "setStats.H"

Info << endl << "Gaussian Random Field generated" << endl;

#include "printStats.H"
#include "printLogStats.H"

Info << endl << "Postprocessing" << endl;
