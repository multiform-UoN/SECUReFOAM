//- Ref:  A handy approximation for the error function and its inverse" by Sergei Winitzki

//- Percentiles
const scalarField PP01(1,0.5);
const scalarField PP02(1,0.5);
const scalarField PP
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("percentiles",PP01)
);
const scalarField PP2
(
  setRandomFieldDict.lookupOrDefault<scalarField> ("percentiles2",PP02)
);
Info << endl << "Percentiles " << PP << "    " << PP2 << endl;

const scalarField SGNTT(sign(2*PP-scalar(1))); 
const scalarField SGNTT2(sign(2*PP2-scalar(1)));

const scalarField LN( log( 4*PP*(scalar(1)-PP)) );
const scalarField LN2( log( 4*PP2*(scalar(1)-PP2)) );

const scalarField AA( scalar(2)/(pi*(0.147)) + 0.5 * LN );
const scalarField AA2( scalar(2)/(pi*(0.147)) + 0.5 * LN2 );

const scalarField BB( scalar(1)/(0.147) * LN );
const scalarField BB2( scalar(1)/(0.147) * LN2 );

const scalarField TT01( mu + sigma*Foam::sqrt(scalar(2))*SGNTT*sqrt(-AA + sqrt(AA*AA - BB)) );
const scalarField TT02( mu + sigma*Foam::sqrt(scalar(2))*SGNTT2*sqrt(-AA2 + sqrt(AA2*AA2 - BB2)) );

