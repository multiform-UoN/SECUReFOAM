if (dispersionModel)
{
    //- Update cell-based dispersion
    D =
        beta
        *
        (
          (Dm + alphaT * mag(U)) * symmTensor::I
          +
          (alphaL-alphaT) * (sqr(U))/(mag(U)+USMALL)
        );

    //- Update face-based dispersion
    Df = inv
          (
              reverseLinear<SymmTensor<scalar>>
                  (
                      mesh
                  ).interpolate(inv(D))
          );
}

