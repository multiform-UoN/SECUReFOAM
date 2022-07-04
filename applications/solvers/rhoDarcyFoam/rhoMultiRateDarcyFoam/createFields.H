#include "../createFields.H"

Info<< "Reading multiRateProperties\n" << endl;

IOdictionary multiRateProperties
(
    IOobject
    (
        "multiRateProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
