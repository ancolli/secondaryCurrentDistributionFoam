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
    Foam::regionCoupleSolidFvPatchScalarField

Description
    Mixed boundary condition for potential in a solid, to be used for secondary 
    distribution problems.

    See paper description: https://doi.org/10.1016/j.electacta.2018.09.121


Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | fi
        j0           | list of exchangeCurrentDensity	| yes | 
	b	     | list of tafel slope	| yes |     
        E0	     | list of equilibriumPotential	| yes |
        kappaMethod  | do not touch | inherited |
        kappa        | solid conductivity | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {

	    type	    regionCoupleSolid;
	    Tnbr            fi;
	    j0		    (-2.09 9.38e-2);
	    b		    (-0.0621347 0.0485);
            E0		    (0.0 1.23);
            kappaMethod     lookup;
            kappa           ks;
            value           $internalField;
        
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

See also
    Foam:: regionCoupleFluidFvPatchScalarField

SourceFiles
   regionCoupleSolidFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleSolidFvPatchScalarField_H
#define regionCoupleSolidFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "potentialCoupledBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupleSolidFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupleSolidFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public potentialCoupledBase
{
    // Private data

        //- Name of field on the neighbour region
        const word TnbrName_;
        //- Kinetic parameters
        scalarList exchangeCurrentDensity_;
        scalarList beta_;
        scalarList equilibriumPotential_;

public:

    //- Runtime type information
    TypeName("regionCoupleSolid");

    // Constructors

        //- Construct from patch and internal field
        regionCoupleSolidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupleSolidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  regionCoupled onto a
        //  new patch
        regionCoupleSolidFvPatchScalarField
        (
            const regionCoupleSolidFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupleSolidFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupleSolidFvPatchScalarField
        (
            const regionCoupleSolidFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupleSolidFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
