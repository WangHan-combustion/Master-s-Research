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

#include "WADESDIT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WADESDIT<BasicTurbulenceModel>::chi() const
{
    return Rnu_/this->nu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WADESDIT<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Aplus_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WADESDIT<BasicTurbulenceModel>::blend
(
	const volScalarField& Switch,
	const dimensionedScalar& psi1,
	const dimensionedScalar& psi2
) const
{
	return Switch*(psi1 - psi2) + psi2;
}

template<class BasicTurbulenceModel>
void WADESDIT<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = Rnu_*fv1;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void WADESDIT<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WADESDIT<BasicTurbulenceModel>::fdes
(
    const volScalarField& S
) const
{
    return sqrt(Rnu_) / (sqrt(S) * CDES_*this->delta());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WADESDIT<BasicTurbulenceModel>::WADESDIT
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
            13.0
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
            0.0833
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

    Switch_
    (
        IOobject
        (
            "Switch",
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
			IOobject::READ_IF_PRESENT,
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
		    IOobject::NO_READ,//MUST_READ,
		    IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("Length", dimensionSet(0, 1, 0, 0, 0), 0.0)
	),
    
    old_
	(
		IOobject
		(
		    "old",
		    this->runTime_.timeName(),
		    this->mesh_,
		    IOobject::NO_READ,//MUST_READ,
		    IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("Zero", dimensionSet(0, -2, 0, 0, 0), 0.0)
	),
    
    new_
	(
		IOobject
		(
		    "new",
		    this->runTime_.timeName(),
		    this->mesh_,
		    IOobject::NO_READ,//MUST_READ,
		    IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("Zero", dimensionSet(0, -2, 0, 0, 0), 0.0)
	)
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WADESDIT<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WADESDIT<BasicTurbulenceModel>::DRnuEff(volScalarField Switch) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRnuEff", Rnu_*sigma(Switch) + this->nu())
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WADESDIT<BasicTurbulenceModel>::k() const
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
tmp<volScalarField> WADESDIT<BasicTurbulenceModel>::epsilon() const
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
void WADESDIT<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    //const volVectorField& U = this->U_;

    LESeddyViscosity<BasicTurbulenceModel>::correct();
    
    volScalarField S2(2.0*magSqr(symm(fvc::grad(this->U_))));
	volScalarField S = sqrt(S2);
	bound(S, dimensionedScalar("0", S.dimensions(), SMALL));
	bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

//////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////

    volScalarField Lvk
    (
            sqrt(S2)
            /(
                    mag(fvc::laplacian(this->U_))
                    + dimensionedScalar
                    (
                            "ROOTVSMALL",
                            dimensionSet(0, -1, -1, 0, 0),
                            ROOTVSMALL
                    )
            )
    );
    
    volScalarField boundEke = sqr(Rnu_) * sqr(1.0/Lvk);
    old_ = magSqr(fvc::grad(S))/S2;
    new_ = sqr(1.0/Lvk);

////////////////////////////////////////////////////////////////
    
	// Calculate and bound delta
	outDelta_ = this->delta();

	// Calculate fdes
    volScalarField fdes(this->fdes(S));
	fdes_ = fdes;
    volScalarField fdes2 = sqr(fdes);
    bound(fdes2,SMALL);
    
////////////////////////////////////////////////////////////////

    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*DRnuEff(Switch_), Rnu_)
     ==
        C1(Switch_)*alpha*rho*S*Rnu_
	  //+ C2kw_*Switch_*alpha*rho*fvm::Sp((fvc::grad(Rnu_) & fvc::grad(S))/S/fdes2,Rnu_)
	  //- (1.0-Switch_)*alpha*rho*fvm::Sp((C2ke_)*Rnu_*magSqr(fvc::grad(S))/S2/fdes2,Rnu_)
      - C2ke_*alpha*rho*boundEke/fdes2
    );

    RnuEqn().relax();
    solve(RnuEqn);
    bound(Rnu_, dimensionedScalar("0", Rnu_.dimensions(), 0.0));
    Rnu_.correctBoundaryConditions();

    correctNut(fv1);
}


// LESRegion = 1 in LES region and =0 in RANS region
template<class BasicTurbulenceModel>
tmp<volScalarField> WADESDIT<BasicTurbulenceModel>::LESRegion() const
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
