Info << "Mean velocity  = "
<< gSum(U.primitiveField()*mesh.V())/gSum(mesh.V())
<< endl;

Info << "Mean velocity fracture = "
<< gSum(U_fr.primitiveField()*mesh.V())/gSum(mesh.V())
<< endl;

scalar MassTot
(
	gSum
	(
		(
			rho.primitiveField()*beta.primitiveField()
			+
			rho_fr.primitiveField()*beta_fr.primitiveField()
		)
		*
		mesh.V()
	)
);

Info << "Total mass  = "<< MassTot << endl;

Info << "Mass rate  = "<< (MassTot - oldMassTot)/runTime.deltaT().value()
<< endl;

scalar MassFracture
(
	gSum(rho_fr.primitiveField()*mesh.V())
);

Info << "Mass rate fracture  = "
<< (MassFracture - oldMassFracture)/runTime.deltaT().value()
<< nl << endl;

oldMassTot = MassTot;
oldMassFracture = MassFracture;
