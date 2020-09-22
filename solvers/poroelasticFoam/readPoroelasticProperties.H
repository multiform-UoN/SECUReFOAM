Info << "Reading poroelastic properties\n" << endl;
IOdictionary poroelasticProperties
(
	IOobject
	(
		"poroelasticProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

//- Read Biot coupling coefficient
dimensionedScalar alpha
(
    "alpha",
    dimless,
    poroelasticProperties.lookup("alpha")
);


dictionary fluidProperties(poroelasticProperties.subDict("fluidProperties"));

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    fluidProperties.lookup("nu")
);

dimensionedScalar rhof
(
    "rhof",
    dimDensity,
    fluidProperties.lookup("rho")
);

dimensionedScalar stor
(
    "stor",
    dimless/dimPressure/dimDensity,
    fluidProperties.lookup("storativity")
);

dictionary solidProperties(poroelasticProperties.subDict("solidProperties"));

//- First Lame' coefficient
dimensionedScalar mu
(
    "mu",
    dimPressure,
    solidProperties.lookup("mu")
);

//- Second Lame' coefficient
dimensionedScalar lambda
(
    "lambda",
    dimPressure,
    solidProperties.lookup("lambda")
);

dimensionedScalar rhos
(
    "rhos",
    dimDensity,
    solidProperties.lookup("rho")
);
