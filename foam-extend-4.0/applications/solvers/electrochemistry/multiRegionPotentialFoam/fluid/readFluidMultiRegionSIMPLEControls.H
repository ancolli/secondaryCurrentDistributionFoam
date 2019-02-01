    dictionary simple = fluidRegions[i].solutionDict().subDict("SIMPLE");

    int nNonOrthCorr = 0;
    if (simple.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(simple.lookup("nNonOrthogonalCorrectors"));
    }

