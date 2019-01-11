/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    wallCurrentFlux

Description
    Calculates and writes current distribution of fi and total current at the wall.

    Details: Electrochimica Acta 290 (2018) 676-685  https://doi.org/10.1016/j.electacta.2018.09.121
    See also: https://github.com/ancolli/secondaryCurrentDistributionFoam

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fiheader
        (
            "fi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject kfheader
        (
            "kf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check fi and kf exists
        if (fiheader.headerOk() && kfheader.headerOk())//
        {
            mesh.readUpdate();

            Info<< "    Reading fi" << endl;
            volScalarField fi(fiheader, mesh);

            Info<< "    Reading kf" << endl;
            volScalarField kf(kfheader, mesh);

            Info<< "    Calculating wallCurrentFlux" << endl;

            volScalarField wallCurrentFlux
            (
                IOobject
                (
                    "wallCurrentFlux",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar   
                (
                    "wallCurrentFlux",
                    fi.dimensions()*kf.dimensions()/dimLength,
                    0.0
                )
            );

        surfaceScalarField currentFlux =
            fvc::interpolate(kf)*fvc::snGrad(fi);

        const surfaceScalarField::GeometricBoundaryField& patchCurrentFlux =
            currentFlux.boundaryField();

        Info<< "\nWall current [A]" << endl;
        forAll(patchCurrentFlux, patchi)
        {
            if (mesh.boundary()[patchi].isWall())
            {
                Info<< mesh.boundary()[patchi].name()
                    << " "
                    << sum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchCurrentFlux[patchi]
                       )
                    << endl;
            }
        }
        Info<< endl;

            forAll(wallCurrentFlux.boundaryField(), patchi)
            {
                wallCurrentFlux.boundaryField()[patchi] = patchCurrentFlux[patchi];
            }


            wallCurrentFlux.write();
        }
        else
        {
            Info<< "    No fi or kf " << endl;
        }
    }

    Info<< "End" << endl;

    return 0;
}
// ************************************************************************* //
