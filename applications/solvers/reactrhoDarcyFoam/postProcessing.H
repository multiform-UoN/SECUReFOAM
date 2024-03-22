//- THIS FILE IS FOR EITHER QUICK DEBUG OR POSTPROCESSING
//- COMMENT OUT FOR FINAL RELEASE

//- Comparison of Boussinesq vs non-Boussinesq
// Info << endl << "dcdt " << (Foam::sqrt(fvc::domainIntegrate(Foam::pow(beta*fvc::ddt(c),2)))/fvc::domainIntegrate(c)).value()
//         << " Ugrad(c) " << (Foam::sqrt(fvc::domainIntegrate(Foam::pow((U&fvc::grad(c)),2)))/fvc::domainIntegrate(c)).value()
//         << " sum " << (Foam::sqrt(fvc::domainIntegrate(Foam::pow((U&fvc::grad(c))+beta*fvc::ddt(c),2)))/fvc::domainIntegrate(c)).value() << endl;
