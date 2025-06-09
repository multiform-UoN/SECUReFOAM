#include "../updateD.H"

if (dispersionModel)
{
    //- Update cell-based dispersion
    D_fr =
    (
      (Dm + alphaT_fr * mag(U_fr)) * symmTensor::I
      +
      (alphaL_fr-alphaT_fr) * (sqr(U_fr))/(mag(U_fr)+USMALL)
    );

    //- Update face-based dispersion
    Df_fr = inv
          (
              reverseLinear<SymmTensor<scalar>>
                  (
                      mesh
                  ).interpolate(inv(D_fr))
          );
}

