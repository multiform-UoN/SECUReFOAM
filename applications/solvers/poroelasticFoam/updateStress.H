Sigma = mu*dev(twoSymm(gradD + T(gradD))) + lambda*tr(gradD)*I;

//- Update permeability to account for fractures
if(frakModel)
{
    forAll(Sigma,cellI)
    {
        if(mag(Sigma[cellI]) >frakSigmaThreshold)
        {
            K[cellI] = tensor::I*frakK;
        }
    }
}
