if (printStats)
{
  // ---- TESTS AND DEBUGGING
  const scalar meank(gSum((*ptrK)*V)/vol); // spatial average
  const scalar vark(gSum(pow((*ptrK),2)*V)/vol); // spatial variance
  Info << endl << "---- Statistics random field ---" << endl;
  Info << "CorrLengths: " << Lcorr << endl;   // TODO - compute these
  Info << "Mean: " << meank << endl;
  Info << "Var: " << (-meank*meank+vark) << endl;
  Info << "--------------------------------" << endl << endl;
}