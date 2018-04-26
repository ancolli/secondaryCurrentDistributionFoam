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
    Foam::regionCoupleFluidFvPatchScalarField

Description
    Mixed boundary condition for potential in the fluid, to be used for secondary 
    distribution problems. 

    See paper description


Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | fi
        j0           | list of exchangeCurrentDensity	| yes | 
	b	     | list of tafel slope	| yes |     
        E0	     | list of equilibriumPotential	| yes |
	kappaMethod  | do not touch | inherited |
        kappa        | fluid conductivity | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {

	    type	    regionCoupleFluid;
	    Tnbr            fi;
	    j0		    (-2.09 9.38e-2);
	    b		    (-0.0621347 0.0485);
            E0		    (0.0 1.23);
            kappaMethod     lookup;
            kappa           kf;
            value           $internalField;
        
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

See also
    Foam::regionCoupleSolidFvPatchScalarField

SourceFiles
    regionCoupleFluidFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleFluidFvPatchScalarField_H
#define regionCoupleFluidFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "tCoupledBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupleFluidFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupleFluidFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public tCoupledBase
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
    TypeName("regionCoupleFluid");


    // Constructors

        //- Construct from patch and internal field
        regionCoupleFluidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupleFluidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        regionCoupleFluidFvPatchScalarField
        (
            const regionCoupleFluidFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupleFluidFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupleFluidFvPatchScalarField
        (
            const regionCoupleFluidFvPatchScalarField&,
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
                new regionCoupleFluidFvPatchScalarField
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
