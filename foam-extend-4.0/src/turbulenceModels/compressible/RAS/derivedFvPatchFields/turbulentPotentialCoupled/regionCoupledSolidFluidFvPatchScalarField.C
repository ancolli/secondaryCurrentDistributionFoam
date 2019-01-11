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

    Details: Electrochimica Acta 290 (2018) 676-685  https://doi.org/10.1016/j.electacta.2018.09.121
    See also: https://github.com/ancolli/secondaryCurrentDistributionFoam

\*---------------------------------------------------------------------------*/

#include "regionCoupledSolidFluidFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "directMappedPatchBase.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regionCoupledSolidFluidFvPatchScalarField::
regionCoupledSolidFluidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_("undefined-nbrField"),
    KappaName_("undefined-Kappa"),
    sideName_("undefined-sideName"),
    exchangeCurrentDensity_(0),
    beta_(0),
    equilibriumPotential_(0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


regionCoupledSolidFluidFvPatchScalarField::
regionCoupledSolidFluidFvPatchScalarField
(
    const regionCoupledSolidFluidFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_),
    KappaName_(ptf.KappaName_),
    sideName_(ptf.sideName_),
    exchangeCurrentDensity_(ptf.exchangeCurrentDensity_),
    beta_(ptf.beta_),
    equilibriumPotential_(ptf.equilibriumPotential_)

{}


regionCoupledSolidFluidFvPatchScalarField::
regionCoupledSolidFluidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_(dict.lookup("nbrField")),
    KappaName_(dict.lookup("Kappa")),
    sideName_(dict.lookup("side")),
    exchangeCurrentDensity_(0),
    beta_(0),
    equilibriumPotential_(0)

{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "regionCoupledSolidFluidFvPatchScalarField::"
            "regionCoupledSolidFluidFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
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


regionCoupledSolidFluidFvPatchScalarField::
regionCoupledSolidFluidFvPatchScalarField
(
    const regionCoupledSolidFluidFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    KappaName_(wtcsf.KappaName_),
    sideName_(wtcsf.sideName_),
    exchangeCurrentDensity_(wtcsf.exchangeCurrentDensity_),
    beta_(wtcsf.beta_),
    equilibriumPotential_(wtcsf.equilibriumPotential_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
regionCoupledSolidFluidFvPatchScalarField::Kappa() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (KappaName_ == "none")
    {
       FatalErrorIn
        (
            ": KappaName can not be none"
            "\n"
        )   
            << exit(FatalError);

	return scalarField(0);
    }
    else if (mesh.objectRegistry::foundObject<volScalarField>(KappaName_))
    {
        return lookupPatchField<volScalarField, scalar>(KappaName_);
    }
    else if (mesh.objectRegistry::foundObject<volSymmTensorField>(KappaName_))
    {
        const symmTensorField& KappaWall =
            lookupPatchField<volSymmTensorField, scalar>(KappaName_);

        vectorField n = patch().nf();

        return n & KappaWall & n;
    }
    else
    {
        FatalErrorIn
        (
            "regionCoupledSolidFluidFvPatchScalarField::Kappa()"
            " const"
        )   << "Did not find field " << KappaName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'Kappa' to 'none', a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void regionCoupledSolidFluidFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the directMappedPatchBase
    const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    // Force recalculation of mapping and schedule
    const mapDistribute& distMap = mpp.map();

    tmp<scalarField> intFld = patchInternalField();


    // Calculate the temperature by harmonic averaging
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const regionCoupledSolidFluidFvPatchScalarField& nbrField =
    refCast
    <
        const regionCoupledSolidFluidFvPatchScalarField
    >
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            neighbourFieldName_
        )
    );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld = nbrField.patchInternalField();
    mapDistribute::distribute
    (
        Pstream::defaultComms(),
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nbrIntFld
    );

    tmp<scalarField> nbrPotencial(new scalarField(nbrField.size(), 0.0));
    nbrPotencial() = nbrField;

    tmp<scalarField> myKDelta = Kappa()*patch().deltaCoeffs();
    tmp<scalarField> potential = patchInternalField() + snGrad()/patch().deltaCoeffs();
   

    tmp<scalarField> Ai = 0*nbrPotencial();
    tmp<scalarField> Bi = 0*nbrPotencial();

    
if (sideName_ == "solid")
{
// For more than one reaction at the same electrode
    if (exchangeCurrentDensity_.size() > 0)
    {
    	forAll(exchangeCurrentDensity_, iReaction)
        {
        	Ai = Ai() + exchangeCurrentDensity_[iReaction]*exp((potential()-nbrPotencial()-equilibriumPotential_[iReaction])/beta_[iReaction])/beta_[iReaction]/myKDelta();
		Bi = Bi() + exchangeCurrentDensity_[iReaction]*exp((potential()-nbrPotencial()-equilibriumPotential_[iReaction])/beta_[iReaction])/myKDelta();
        }
    }

// Explanation in the paper

    this->refValue() = potential()-Bi()/Ai();

} 
else if(sideName_ == "fluid")
{
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
}

else 
{
	FatalErrorIn
        (
            ": In changeDictionaryDict, side should be either solid or fluid"
            "\n"
        )   
            << exit(FatalError);
}

    this->refGrad() = 0;
    this->valueFraction() = Ai()/(1+Ai());


    mixedFvPatchScalarField::updateCoeffs();


    if (debug)
    {
        scalar Q = gSum(Kappa()*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " -> "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heatFlux:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void regionCoupledSolidFluidFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nbrField")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("Kappa") << KappaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("side")<< sideName_ << token::END_STATEMENT << nl;
    exchangeCurrentDensity_.writeEntry("j0", os);
    beta_.writeEntry("b", os);
    equilibriumPotential_.writeEntry("E0", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    regionCoupledSolidFluidFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
