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

#include "SSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SSTDES<BasicTurbulenceModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SSTDES<BasicTurbulenceModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SSTDES<BasicTurbulenceModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SSTDES<BasicTurbulenceModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23() *= F3();
    }

    return f23;
}


template<class BasicTurbulenceModel>
void SSTDES<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SSTDES<BasicTurbulenceModel>::rd
(
    const volScalarField& magGradU
) const
{
    tmp<volScalarField> tr
    (
        min
        (
            this->nuEff()
           /(
               max
               (
                   magGradU,
                   dimensionedScalar("SMALL", magGradU.dimensions(), SMALL)
               )
              *sqr(this->kappa_*this->y_)
            ),
            scalar(10)
        )
    );
    tr().boundaryField() == 0.0;

    return tr;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SSTDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    //return 1 - tanh(pow3(8*rd(magGradU)));
	return 1 - tanh(pow3(20*rd(magGradU)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SSTDES<BasicTurbulenceModel>::fdes
(
    const volTensorField& gradU
)
{
    tmp<volScalarField> tdTilda(CDES_*this->delta());
    min(tdTilda().dimensionedInternalField(), tdTilda());
   
	fd_ = fd(mag(gradU));

	if(DDES_)
	{
		return scalar(1) - fd_ * min
                                (
                                    scalar(1) - sqrt(k_) / (betaStar_ * omega_ * tdTilda),
                                    scalar(0)
                                );
	}
	else
	{
		return max
		(
			sqrt(k_) / (betaStar_ * omega_ * tdTilda),
			scalar(1)
		);
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SSTDES<BasicTurbulenceModel>::SSTDES
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

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),
	RASDict_(this->subOrEmptyDict("RAS")),
	omegaMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "omegaMin",
            RASDict_,
            dimless/dimTime,
            SMALL
        )
    ),
    CkeDES_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CkeDES",
            this->coeffDict_,
            0.61
        )
    ),
    CkwDES_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CkwDES_",
            this->coeffDict_,
            0.78
        )
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
	DDES_
    (
        Switch::lookupOrAddToDict
        (
            "DDES",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
	CDES_
	(
		IOobject
		(
			"CDES",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
	),
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
	F1_
	(
		IOobject
		(
			"F1",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
	),
	LESregion_
	(
		IOobject
		(
			"LESregion",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
	),
	fd_
	(
		IOobject
		(
			"fd",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
	)

{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

	if (DDES_)
	{
		Info<< "Delayed Detached-Eddy Simulation (DDES) mode on" << endl << endl;
	}
	else
	{
		Info<< "Delayed Detached-Eddy Simulation (DDES) mode off" << endl << endl;
	}
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SSTDES<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void SSTDES<BasicTurbulenceModel>::correct()
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
    volScalarField& nut = this->nut_;

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField GbyNu((tgradU() && dev(twoSymm(tgradU()))));
    volScalarField G(this->GName(), nut*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
	F1_ = F1;

    
    volScalarField gamma(this->gamma(F1));
    volScalarField beta(this->beta(F1));

	// Calculate fdes_
	CDES_ = blend(F1, CkwDES_, CkeDES_); // CDES_ = F1*CkwDES_ + (1-F1)*CkeDES_
	const volScalarField fdes(this->fdes(fvc::grad(U)));
	fdes_ = fdes;

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
     ==
        alpha*rho*gamma
       *min
        (
            GbyNu,
            (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, b1_*F23()*sqrt(S2))
        )
      - fvm::SuSp((2.0/3.0)*alpha*rho*gamma*divU, omega_)
      - fvm::Sp(alpha*rho*beta*omega_, omega_)
      - fvm::SuSp
        (
            alpha*rho*(F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, this->omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        min(alpha*rho*G, (c1_*betaStar_)*alpha*rho*k_*omega_)
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*betaStar_*omega_*fdes_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, this->kMin_);

    correctNut(S2);

	LESregion_ = neg(scalar(1) - fdes_);
}


// LESRegion = 1 in LES region and =0 in RANS region
template<class BasicTurbulenceModel>
tmp<volScalarField> SSTDES<BasicTurbulenceModel>::LESRegion() const
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
