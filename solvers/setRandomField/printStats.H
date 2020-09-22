if (printStats)
{
  // ---- TESTS AND DEBUGGING
  const scalar meank(gSum((*ptrK)*mesh.V())/vol); // spatial average
  const scalar vark(gSum(pow((*ptrK),2)*mesh.V())/vol); // spatial variance
  Info << endl << "---- Statistics random field ---" << endl;
  Info << "CorrLengths: " << Lcorr << endl;   // TODO - compute these
  Info << "Mean: " << meank << endl;
  Info << "Var: " << (-meank*meank+vark) << endl;
  Info << "--------------------------------" << endl << endl;

  // ///////////
  // // -- TESTS POSTPROCESSING
  // ///////////
  // const scalar eps(1e-8);
  // const label nb(10);   // number of bins
  // const label nd(100);  // number of radii
  // const scalar delta(maxlx/nd + eps); // min radius
  // const scalar maxk(gMax(*ptrK) + eps);
  // const scalar mink(gMin(*ptrK) - eps);
  // const scalar deltak((maxk-mink)/nb + eps); // bin size
  // scalarField C(nd,scalar(0));
  // scalarField V(nd,scalar(0));
  // // scalarField H(nd,scalar(0));
  // scalarField HH(nd,scalar(0));
  // // scalarField HV(nd,scalar(0));
  // RectangularMatrix<scalar> Hpdf(nd,nb,scalar(0));
  // forAll(*ptrK,celli)
  // {
  //   forAll(*ptrK,cellj)
  //   {
  //     const scalar dist(mag( mesh.C()[celli] - mesh.C()[cellj] ));
  //     const label ind(floor(dist/delta));
  //     if (ind>=nd)
  //     {
  //       continue;
  //     }
  //     // two-point correlation
  //     C[ind] +=
  //               ((*ptrK)[celli]-meank)
  //               *
  //               ((*ptrK)[cellj]-meank)
  //               *
  //               mesh.V()[celli]*mesh.V()[cellj]
  //               /
  //               vark;
  //     V[ind] += mesh.V()[celli]*mesh.V()[cellj];
  //     // // entropy using the field itself normalised to 0-1
  //     // const scalar Hi =
  //     //           Foam::log(((*ptrK)[celli]-mink)/maxk)
  //     //           *
  //     //           ((*ptrK)[celli]-mink)/maxk
  //     //           *
  //     //           mesh.V()[celli]*mesh.V()[cellj];
  //     // for (int kk=0; kk<ind; kk++)
  //     // {
  //     //   H[kk] -= Hi;
  //     //   HV[kk] += mesh.V()[celli]*mesh.V()[cellj];
  //     // }
  //     // proper entropy
  //     const label indb(floor(((*ptrK)[cellj]-mink)/deltak));
  //     Hpdf(ind,indb) += mesh.V()[cellj]/vol;
  //   }
  //   for (int jj=0; jj<nd; jj++)
  //   {
  //   for (int kk=0; kk<nb; kk++)
  //   {
  //     if (Hpdf(jj,kk)<eps)
  //     {
  //       continue;
  //     }
  //     HH[jj] -= Hpdf(jj,kk)*Foam::log(Hpdf(jj,kk));
  //   }
  //     HH[jj] *= mesh.V()[celli]/vol;
  //   }
  // }
  // Info << "Two-point correlation" << endl << C/V  << endl;
  // Info << "Entropy" << endl << HH  << endl;
  // ///////////
  // // -- END
  // ///////////
}
