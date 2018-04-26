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

Class
    Foam::temperatureCoupledBase

Description
    Common functions used in potential coupled boundaries.


SourceFiles
    tCoupledBase.C

\*---------------------------------------------------------------------------*/

#ifndef tCoupledBase_H
#define tCoupledBase_H

#include "scalarField.H"
#include "NamedEnum.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class temperatureCoupledBase Declaration
\*---------------------------------------------------------------------------*/

class tCoupledBase
{
public:

    // Public enumerations

        //- Type of supplied Kappa
        enum KMethodType
        {
            mtLookup
        };


protected:

    // Protected data

        static const NamedEnum<KMethodType, 1> KMethodTypeNames_;

        //- Underlying patch
        const fvPatch& patch_;

        //- How to get K
        const KMethodType method_;

        //- Name of thermal conductivity field (if looked up from database)
        const word kappaName_;

        //- Name of the non-Isotropic alpha (default: Anialpha)
        const word alphaAniName_;


public:

    // Constructors

        //- Construct from patch and K name
        tCoupledBase
        (
            const fvPatch& patch,
            const word& calculationMethod,
            const word& kappaName,
            const word& alphaAniName
        );

        //- Construct from patch and dictionary
        tCoupledBase
        (
            const fvPatch& patch,
            const dictionary& dict
        );

         //- Construct from patch and tCoupledBase
        tCoupledBase
        (
            const fvPatch& patch,
            const tCoupledBase& base
        );


    // Member functions

        //- Method to obtain K
        word KMethod() const
        {
            return KMethodTypeNames_[method_];
        }

        //- Name of thermal conductivity field
        const word& kappaName() const
        {
            return kappaName_;
        }

        //- Given patch temperature calculate corresponding K field
        tmp<scalarField> kappa(const scalarField& Tp) const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
