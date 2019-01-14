{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix fifEqn
        (
            (
               - fvm::laplacian(kf, fif) 
            )
          
        );

        fifEqn.relax();

        eqnResidual = fifEqn.solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);

     
    }
}


Info<< "Min/max fif:" << min(fif).value() << ' '
    << max(fif).value() << endl;
