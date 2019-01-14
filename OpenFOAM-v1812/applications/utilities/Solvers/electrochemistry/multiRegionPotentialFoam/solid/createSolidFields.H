    // Initialise solid field pointer lists
PtrList<coordinateSystem> coordinates(solidRegions.size());
PtrList<volScalarField> kss(solidRegions.size());
PtrList<volScalarField> fiss(solidRegions.size());

    // Populate solid field pointer lists
forAll(solidRegions, i)
{
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;


        Info<< "    Adding to kss\n" << endl;
        kss.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "ks",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        Info<< "    Adding to fiss\n" << endl;
        fiss.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fi",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
       );
}

