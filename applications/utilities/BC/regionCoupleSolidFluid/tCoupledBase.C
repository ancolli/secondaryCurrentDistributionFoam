/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is NOT part of OpenFOAM.

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

#include "tCoupledBase.H"
#include "volFields.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::tCoupledBase::KMethodType,
        1
    >::names[] =
    {
        "lookup"
    };
}


const Foam::NamedEnum<Foam::tCoupledBase::KMethodType, 1> 
    Foam::tCoupledBase::KMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tCoupledBase::tCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const word& kappaName,
    const word& alphaAniName
)
:
    patch_(patch),
    method_(KMethodTypeNames_[calculationType]),
    kappaName_(kappaName),
    alphaAniName_(alphaAniName)
{}


Foam::tCoupledBase::tCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(KMethodTypeNames_.read(dict.lookup("kappaMethod"))),
    kappaName_(dict.lookupOrDefault<word>("kappa", "none")),
    alphaAniName_(dict.lookupOrDefault<word>("alphaAni","Anialpha"))
{}


Foam::tCoupledBase::tCoupledBase
(
    const fvPatch& patch,
    const tCoupledBase& base
)
:
    patch_(patch),
    method_(base.method_),
    kappaName_(base.kappaName_),
    alphaAniName_(base.alphaAniName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::tCoupledBase::kappa
(
    const scalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    
    switch (method_)
    {
        
        case mtLookup:
        {
            if (mesh.foundObject<volScalarField>(kappaName_))
            {
                return patch_.lookupPatchField<volScalarField, scalar>
                (
                    kappaName_
                );
            }
            
            else
            {
                FatalErrorInFunction
                    << "Did not find field " << kappaName_
                    << " on mesh " << mesh.name() << " patch " << patch_.name()
                    << nl
                    << "Please set 'kappa' to one of "
                    << KMethodTypeNames_.toc()
                    << " and 'kappaName' to the name of the volScalar"
                    << " or volSymmTensor field (if kappa=lookup)"
                    << exit(FatalError);
            }

            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << KMethodTypeNames_[method_] << nl
                << "Please set 'kappa' to one of " << KMethodTypeNames_.toc()
                << " and 'kappaName' to the name of the volScalar"
                << " or volSymmTensor field (if kappa=lookup)"
                << exit(FatalError);
        }
    }

    return scalarField(0);
}


void Foam::tCoupledBase::write(Ostream& os) const
{
    os.writeKeyword("kappaMethod") << KMethodTypeNames_[method_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphaAni") << alphaAniName_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
