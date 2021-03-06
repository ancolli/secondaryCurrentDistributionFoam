/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::regionCoupledSolidFluidFvPatchScalarField

Description
    Mixed boundary condition for potentials in a solids-fluids, to be used for secondary 
    distribution problems.

    Details: Electrochimica Acta 290 (2018) 676-685  https://doi.org/10.1016/j.electacta.2018.09.121
    See also: https://github.com/ancolli/secondaryCurrentDistributionFoam


Usage
    \table
        Property     | Description             | Required    | Default value
        side         | name of the side    | yes | solid or fluid
        nbrField     | name of the field    | no | fi
        j0           | list of exchangeCurrentDensity	| yes | 
	b	     | list of tafel slope	| yes |     
        E0	     | list of equilibriumPotential	| yes |
        kappaMethod  | do not touch | inherited |
        kappa        | electrical conductivity | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {

	    type	    regionCoupleSolidFluid;
	    side	    solid or fluid;
	    nbrField        fi;
	    j0		    (-2.09 9.38e-2);
	    b		    (-0.0621347 0.0485);
            E0		    (0.0 1.23);
            kappaMethod     lookup;
            kappa           ks or kf;
            value           $internalField;
        
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

See also
    Foam::compressible::turbulentTemperatureCoupledBaffleMixed

SourceFiles
   regionCoupledSolidFluidFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledSolidFluidFvPatchScalarField_H
#define regionCoupledSolidFluidFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupledSolidFluidFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledSolidFluidFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;

        //- Name of electrical conductivity field
        const word KappaName_;

        //- Name of the region side
        const word sideName_; 

        //- Kinetic parameters
        scalarList exchangeCurrentDensity_;
        scalarList beta_;
        scalarList equilibriumPotential_;


public:

    //- Runtime type information
    TypeName("regionCoupledSolidFluid");


    // Constructors

        //- Construct from patch and internal field
        regionCoupledSolidFluidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledSolidFluidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  regionCoupledSolidFluidFvPatchScalarField onto a
        //  new patch
        regionCoupledSolidFluidFvPatchScalarField
        (
            const regionCoupledSolidFluidFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupledSolidFluidFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledSolidFluidFvPatchScalarField
        (
            const regionCoupledSolidFluidFvPatchScalarField&,
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
                new regionCoupledSolidFluidFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding Kappa field
        tmp<scalarField> Kappa() const;

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
