ff = pow(fx[kx],2) + pow(fy[ky],2) + pow(fz[kz],2) ;

// Using analytical Fourier Transform
// This would include prefactors that depends on dx,dy,dz
// but we remove then as everything will be scaled afterwards

if (corrtype=="exponential")
{
  aa =  Foam::pow((ff+scalar(1)),-scalar(dim+1)/scalar(2)) ;
}
else if (corrtype=="gaussian")
{
  aa = Foam::exp(-ff/pi) ;
}
else if (corrtype=="matern")
{
  aa = Foam::pow((scalar(2)*nu+ff),-(nu+dim/scalar(2))) ;
}
else
{
  FatalError
  << "I cannot recognise correlation " << corrtype << nl
  << "Avaliable values are: exponential, gaussian, matern " << nl
  << exit(FatalError);
}

aa = Foam::sqrt(aa);