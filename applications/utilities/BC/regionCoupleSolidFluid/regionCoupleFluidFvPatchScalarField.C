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

#include "regionCoupleFluidFvPatchScalarField.H"
#include "regionCoupleSolidFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regionCoupleFluidFvPatchScalarField::
regionCoupleFluidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    potentialCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    TnbrName_("undefined-Tnbr"),
    exchangeCurrentDensity_(0),
    beta_(0),
    equilibriumPotential_(0)
  
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


regionCoupleFluidFvPatchScalarField::
regionCoupleFluidFvPatchScalarField
(
    const regionCoupleFluidFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    potentialCoupledBase(patch(), ptf),
    TnbrName_(ptf.TnbrName_),
    exchangeCurrentDensity_(ptf.exchangeCurrentDensity_),
    beta_(ptf.beta_),
    equilibriumPotential_(ptf.equilibriumPotential_)

{}


regionCoupleFluidFvPatchScalarField::
regionCoupleFluidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    potentialCoupledBase(patch(), dict),
    TnbrName_(dict.lookup("Tnbr")),
    exchangeCurrentDensity_(0),
    beta_(0),
    equilibriumPotential_(0)
 
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    // search for kinetic parameters
    {
        dict.lookup("j0") >> exchangeCurrentDensity_;
        dict.lookup("b") >> beta_;
	dict.lookup("E0") >> equilibriumPotential_;
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


regionCoupleFluidFvPatchScalarField::
regionCoupleFluidFvPatchScalarField
(
    const regionCoupleFluidFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    potentialCoupledBase(patch(), wtcsf),
    TnbrName_(wtcsf.TnbrName_),
    exchangeCurrentDensity_(wtcsf.exchangeCurrentDensity_),
    beta_(wtcsf.beta_),
    equilibriumPotential_(wtcsf.equilibriumPotential_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void regionCoupleFluidFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    // Calculate potentials
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const regionCoupleSolidFvPatchScalarField& nbrField =
    refCast
    <
        const regionCoupleSolidFvPatchScalarField
    >
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            TnbrName_
        )
    );

// Swap to obtain full local values of neighbour internal field

    tmp<scalarField> nbrPotencial(new scalarField(nbrField.size(), 0.0));
  
    nbrPotencial.ref() = nbrField;

    tmp<scalarField> myKDelta = kappa(*this)*patch().deltaCoeffs();
    tmp<scalarField> potential = patchInternalField() + snGrad()/patch().deltaCoeffs();		
    
    mpp.distribute(nbrPotencial.ref());
  
    
    tmp<scalarField> Ai = 0*nbrPotencial();
    tmp<scalarField> Bi = 0*nbrPotencial();

// For more than one reaction at the same electrode
    if (exchangeCurrentDensity_.size() > 0)
    {
    	forAll(exchangeCurrentDensity_, iReaction)
        {
        	Ai = Ai() + exchangeCurrentDensity_[iReaction]*exp((nbrPotencial()-potential()-equilibriumPotential_[iReaction])/beta_[iReaction])/beta_[iReaction]/myKDelta();
		Bi = Bi() + exchangeCurrentDensity_[iReaction]*exp((nbrPotencial()-potential()-equilibriumPotential_[iReaction])/beta_[iReaction])/myKDelta();
        }
    }

// Explanation in the paper

    this->refValue() = potential()+Bi()/Ai();
    this->refGrad() = 0;
    this->valueFraction() = Ai()/(Ai()+1);

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :"
            << " current:" << Q
            << " wallPotential "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}

void regionCoupleFluidFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("Tnbr")<< TnbrName_
        << token::END_STATEMENT << nl;
    exchangeCurrentDensity_.writeEntry("j0", os);
    beta_.writeEntry("b", os);
    equilibriumPotential_.writeEntry("E0", os);
    potentialCoupledBase::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    regionCoupleFluidFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
