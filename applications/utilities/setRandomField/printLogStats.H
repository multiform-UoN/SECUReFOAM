if (printLogStats)
{
  if (gMin(*ptrK)>scalar(0))
  {
    scalarField logField(log(*ptrK));
    // ---- TESTS AND DEBUGGING
    const scalar meank(gSum(logField*V)/vol); // spatial average
    const scalar vark(gSum(pow(logField,2)*V)/vol); // spatial variance
    Info << endl << "---- Statistics log of random field ---" << endl;
    Info << "CorrLengths: " << Lcorr << endl; // TODO - compute these
    Info << "Mean: " << meank << endl;
    Info << "Var: " << (-meank*meank+vark) << endl;
    Info <<  "------------------------------------" << endl << endl;
  }
  else
  {
    Info << endl
         << "Cannot print logStats as the field takes non-positive values"
         << endl;
  }
}
