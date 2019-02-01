    // Initialise fluid field pointer lists
    PtrList<volScalarField> kfs(fluidRegions.size());
    PtrList<volScalarField> fifs(fluidRegions.size());


    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to kfs\n" << endl;
        kfs.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kf",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to fifs\n" << endl;
        fifs.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );


    }


