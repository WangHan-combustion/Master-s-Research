/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "WrayAgarwalDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDES<BasicTurbulenceModel>::chi() const
{
    return Rnu_/this->nu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDES<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Aplus_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDES<BasicTurbulenceModel>::blend
(
	const volScalarField& Switch,
	const dimensionedScalar& psi1,
	const dimensionedScalar& psi2
) const
{
	return Switch*(psi1 - psi2) + psi2;
}

template<class BasicTurbulenceModel>
void WrayAgarwalDES<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = Rnu_*fv1;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void WrayAgarwalDES<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}

template<class BasicTurbulenceModel>
void WrayAgarwalDES<BasicTurbulenceModel>::precalculations
(
	const volScalarField& S,
    const volTensorField& gradU
) {}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDES<BasicTurbulenceModel>::fdes
(
    const volScalarField& S,			
    const volTensorField& gradU
) const
{
    
    tmp<volScalarField> tdTilda(CDES_*this->delta());
    min(tdTilda().dimensionedInternalField(), tdTilda());
    
    return max
	(
		sqrt(Rnu_) / (sqrt(S) * tdTilda),
		scalar(1)
	);
}

template<class BasicTurbulenceModel>
void WrayAgarwalDES<BasicTurbulenceModel>::blendfactor()
{
	blendfactor_ = neg(scalar(1) - fdes_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwalDES<BasicTurbulenceModel>::WrayAgarwalDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    Aplus_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Aplus",
            this->coeffDict_,
            8.541426//13.0
        )
    ),

    C1ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1ke",
            this->coeffDict_,
			0.1127
        )
    ),

    C1kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1kw",
            this->coeffDict_,
			0.0829
        )
    ),

    sigmake_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmake",
            this->coeffDict_,
			1.0
        )
    ),

    sigmakw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmakw",
            this->coeffDict_,
			0.72
        )
    ),

    C2ke_(C1ke_/sqr(kappa_)+sigmake_),

    C2kw_(C1kw_/sqr(kappa_)+sigmakw_),
    
    CDES_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDES",
            this->coeffDict_,
            0.41
        )
    ),

    Rnu_
    (
        IOobject
        (
            "Rnu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    Switch1_
    (
        IOobject
        (
            "Switch1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    y_(wallDist::New(this->mesh_).y()),

	fdes_
	(
		IOobject
		(
			"fdes",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
	),

	outDelta_
	(
		IOobject
		(
		    "outDelta",
		    this->runTime_.timeName(),
		    this->mesh_,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("Length", dimensionSet(0, 1, 0, 0, 0), 0.0)
	),

	blendfactor_
	(
		IOobject
		(
			"blendfactor",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0)
	),

	UBlendingFactor_
	(
		IOobject
		(
		    "UBlendingFactor",
		    this->runTime_.timeName(),
		    this->mesh_,
		    IOobject::READ_IF_PRESENT,
		    IOobject::NO_WRITE
		),
		fvc::interpolate(blendfactor_)
	),

	RnuBlendingFactor_
	(
		IOobject
		(
		    "RnuBlendingFactor",
		    this->runTime_.timeName(),
		    this->mesh_,
		    IOobject::READ_IF_PRESENT,
		    IOobject::NO_WRITE
		),
		fvc::interpolate(blendfactor_)
	),

	pBlendingFactor_
	(
		IOobject
		(
		    "pBlendingFactor",
		    this->runTime_.timeName(),
		    this->mesh_,
		    IOobject::READ_IF_PRESENT,
		    IOobject::NO_WRITE
		),
		fvc::interpolate(blendfactor_)
	),

	KBlendingFactor_
	(
		IOobject
		(
		    "KBlendingFactor",
		    this->runTime_.timeName(),
		    this->mesh_,
		    IOobject::READ_IF_PRESENT,
		    IOobject::NO_WRITE
		),
		fvc::interpolate(blendfactor_)
	),

	eBlendingFactor_
	(
		IOobject
		(
		    "eBlendingFactor",
		    this->runTime_.timeName(),
		    this->mesh_,
		    IOobject::READ_IF_PRESENT,
		    IOobject::NO_WRITE
		),
		fvc::interpolate(blendfactor_)
	),

	hBlendingFactor_
	(
		IOobject
		(
		    "hBlendingFactor",
		    this->runTime_.timeName(),
		    this->mesh_,
		    IOobject::READ_IF_PRESENT,
		    IOobject::NO_WRITE
		),
		fvc::interpolate(blendfactor_)
	)
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WrayAgarwalDES<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        sigma_.readIfPresent(this->coeffDict());
        CDES_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDES<BasicTurbulenceModel>::DRnuEff(volScalarField Switch) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRnuEff", Rnu_*sigma(Switch) + this->nu())
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDES<BasicTurbulenceModel>::k() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDES<BasicTurbulenceModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
        )
    );
}


template<class BasicTurbulenceModel>
void WrayAgarwalDES<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;

    LESeddyViscosity<BasicTurbulenceModel>::correct();
    
	// Calculate strain rate magnitude
    volScalarField S2(2.0*magSqr(symm(fvc::grad(U))));
	volScalarField S = sqrt(S2);
	bound(S, dimensionedScalar("0", S.dimensions(), SMALL));
	bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));

	// Calculate chi and fv1
    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

	// Calculate the switch function
	volScalarField eta = y_*sqrt(Rnu_*S)/(20.0*this->nu());
	volScalarField arg1 = (1.0+18.0*eta)/(1.0+sqr(eta));
	Switch1_ = tanh(pow((1.0+20.0*eta)/(1.0+sqr(y_*max(sqrt(Rnu_*S),dimensionedScalar("Min", dimensionSet(0, 1, -1, 0, 0), 1.5))/(20.0*this->nu()))),4.0));
	Switch1_ = min(Switch1_,0.9);
	bound(Switch1_,SMALL);

////////////////////////////////////////////////////////////////

	// Calculate velocity gradient
    tmp<volTensorField> tgradU = fvc::grad(U);

	// Do some pre-calculations
   	precalculations(S, tgradU());

	// Calculate fdes
	const volScalarField fdes(this->fdes(S, tgradU()));
    const volScalarField fdes2 = sqr(fdes);
	fdes_ = fdes;
    
    // Blend Scheme
	this->blendfactor();

	// Assigned the blendfactor
    UBlendingFactor_ = fvc::interpolate(blendfactor_);
    RnuBlendingFactor_ = UBlendingFactor_;
    pBlendingFactor_ = UBlendingFactor_;
    KBlendingFactor_ = UBlendingFactor_;
    eBlendingFactor_ = UBlendingFactor_;
    hBlendingFactor_ = UBlendingFactor_;
    
////////////////////////////////////////////////////////////////

    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*DRnuEff(Switch1_), Rnu_)
     ==
        C1(Switch1_)*alpha*rho*S*Rnu_
	  + C2kw_*Switch1_*alpha*rho*fvm::Sp((fvc::grad(Rnu_) & fvc::grad(S))/S/fdes2,Rnu_)
	  - (1.0-Switch1_)*alpha*rho*fvm::Sp((C2ke_)*Rnu_*magSqr(fvc::grad(S))/S2/fdes2,Rnu_)
    );

    RnuEqn().relax();
    solve(RnuEqn);
    bound(Rnu_, dimensionedScalar("0", Rnu_.dimensions(), 0.0));
    Rnu_.correctBoundaryConditions();

    correctNut(fv1);
	
	// Calculate and bound delta
	outDelta_ = this->delta();
}


// LESRegion = 1 in LES region and =0 in RANS region
template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDES<BasicTurbulenceModel>::LESRegion() const
{
    tmp<volScalarField> tLESRegion
    (
        new volScalarField
        (
            IOobject
            (
                "DES::LESRegion",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //neg(dTilda(chi, fv1, fvc::grad(this->U_)) - y_)
            neg(scalar(1) - fdes_)
        )
    );

    return tLESRegion;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
