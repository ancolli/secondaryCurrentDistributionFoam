/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "wallCurrentFlux.H"

#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallCurrentFlux, 0);
    addToRunTimeSelectionTable(functionObject, wallCurrentFlux, dictionary);
}
}

void Foam::functionObjects::wallCurrentFlux::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "Wall flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min current density (A/m2)");
    writeTabbed(os, "max current density (A/m2)");
    writeTabbed(os, "current (A)");
    writeTabbed(os, "Absolute current (A)");
    os << endl;
}

void Foam::functionObjects::wallCurrentFlux::calcFlux
(
    const volScalarField& kf_,
    const volScalarField& fi_, 
    volScalarField& wallCurrentFlux
)

{
    
    surfaceScalarField flux
    (
      fvc::interpolate(kf_)*fvc::snGrad(fi_)
    );

    volScalarField::Boundary& wallCurrentFluxBf =
        wallCurrentFlux.boundaryFieldRef();

    const surfaceScalarField::Boundary& fluxBf =
        flux.boundaryField();

    forAll(wallCurrentFluxBf, patchi)
    {
        wallCurrentFluxBf[patchi] = fluxBf[patchi];
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::wallCurrentFlux::wallCurrentFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    patchSet_(),


    kf_
    (
       IOobject
       (
           "kf",
           
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
    ),
 
    fi_
    (
       IOobject
       (
          "fi",
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
      
    )
    

    
{
    volScalarField* wallCurrentFluxPtr
    (
	new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", fi_.dimensions()*kf_.dimensions()/dimLength, 0.0)
        )
    );
 
    mesh_.objectRegistry::store(wallCurrentFluxPtr);  
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallCurrentFlux::~wallCurrentFlux()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallCurrentFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);
     
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::wallCurrentFlux::execute()
{
    volScalarField& wallCurrentFlux = lookupObjectRef<volScalarField>(type());

    volScalarField fi_
    (
       IOobject
       (
          "fi",
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
      
    );

    calcFlux(kf_, fi_, wallCurrentFlux);

    return true;
}


bool Foam::functionObjects::wallCurrentFlux::write()
{
    const volScalarField& wallCurrentFlux = lookupObject<volScalarField>(type()); 

       Log << type() << " " << name() << " write:" << nl
        << "    writing field " << wallCurrentFlux.name() << endl;

    wallCurrentFlux.write();

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    
    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallCurrentFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);
	const scalar currentAbs = gSum(mag(magSf[patchi]*hfp)); //absolute value

        if (Pstream::master())
        {

	    writeTime(file());

            file()
                << mesh_.time().value()
                << token::TAB << pp.name()
                << token::TAB << minHfp
                << token::TAB << maxHfp
                << token::TAB << integralHfp
		<< token::TAB << currentAbs
                << endl;
        }

        Log << "    min current density (A/m2)/max current density (A/m2)/current (A)/ Absolute current (A) (" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << ", " << currentAbs << endl;
    }

    return true;
}


// ************************************************************************* //
