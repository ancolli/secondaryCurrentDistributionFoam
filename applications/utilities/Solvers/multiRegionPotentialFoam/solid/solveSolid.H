{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix fisEqn
        (
            (
                 - fvm::laplacian(ks, fis)
            )
          
        );

        fisEqn.relax();
	
	eqnResidual = fisEqn.solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);

     
    }
}


Info<< "Min/max fis:" << min(fis).value() << ' '
    << max(fis).value() << endl;
