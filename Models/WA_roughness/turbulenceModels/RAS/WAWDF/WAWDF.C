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

#include "WAWDF.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::chi() const
{
    return Rnu_/this->nu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cw_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::WDF_R
(
    const volScalarField& S,
    const volScalarField& W   
) const
{
    return mag(W/S);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::WDF_omega
(
    const volScalarField& S 
) const
{
    return S/sqrt(Cmu_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::WDF_k
(
    const volScalarField& omega
) const
{
    return this->nut_*omega;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::arg1
(
    const volScalarField& S,
    const volScalarField& W
) const
{
    const volScalarField R = WDF_R(S, W);
    const volScalarField omega = WDF_omega(S);
    const volScalarField k = WDF_k(omega);

    const volScalarField eta = S*max(1.0, R);

    return (this->nu()+Rnu_)/2 * sqr(eta)/max(Cmu_*k*omega,
                                              dimensionedScalar("SMALL", 
                                                                dimensionSet(0, 2, -3, 0, 0), 
                                                                SMALL)
                                             );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::calcSwitch
(
    const volScalarField& S,
    const volScalarField& W    
) const
{
    return tanh(pow(Cs1_*arg1(S, W), Cs2_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::blend
(
	const volScalarField& Switch,
	const dimensionedScalar& psi1,
	const dimensionedScalar& psi2
) const
{
	return Switch*(psi1 - psi2) + psi2;
}

template<class BasicTurbulenceModel>
void WAWDF<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = Rnu_*fv1;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void WAWDF<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WAWDF<BasicTurbulenceModel>::WAWDF
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
    eddyViscosity<RASModel<BasicTurbulenceModel> >
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
    Cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw",
            this->coeffDict_,
            12.5
        )
    ),

    C1ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1ke",
            this->coeffDict_,
			0.338
        )
    ),

    C1kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1kw",
            this->coeffDict_,
			0.083
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

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
			0.09
        )
    ),

    Cs1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs1",
            this->coeffDict_,
			2.0
        )
    ),

    Cs2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs2",
            this->coeffDict_,
			4.0
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

    f1_
    (
        IOobject
        (
            "f1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    S_
    (
        IOobject
        (
            "S",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WAWDF<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {   
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::DRnuEff(volScalarField Switch) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRnuEff", Rnu_*sigma(Switch) + this->nu())
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::k() const
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
tmp<volScalarField> WAWDF<BasicTurbulenceModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << nl;

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
void WAWDF<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    // Calculate strain rate magnitude S
	volScalarField S2(2.0*magSqr(symm(fvc::grad(this->U_))));
	volScalarField S = sqrt(S2);
	bound(S, dimensionedScalar("0", S.dimensions(), SMALL)); // SMALL = 1e-15
	bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));
	S_ = S;

    // Calculate vorticity magnitude W
    volScalarField W2(2.0*magSqr(skew(fvc::grad(this->U_))));
    volScalarField W = sqrt(W2);
	bound(W, dimensionedScalar("0", W.dimensions(), SMALL));
	bound(W2, dimensionedScalar("0", W2.dimensions(), SMALL));

    // Calculate switch function (f1)
    f1_ = calcSwitch(S, W);

    // Define and solve R-Equation
    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*DRnuEff(f1_), Rnu_)
     ==
        C1(f1_)        *alpha*rho*S*Rnu_
      + f1_      *C2kw_*alpha*rho*fvm::Sp((fvc::grad(Rnu_)&fvc::grad(S))/S, Rnu_)
	  - (1.0-f1_)*C2ke_*alpha*rho*fvm::Sp(Rnu_*magSqr(fvc::grad(S))/S2, Rnu_)
    );

    RnuEqn().relax();
    solve(RnuEqn);
    bound(Rnu_, dimensionedScalar("0", Rnu_.dimensions(), 0.0));
    Rnu_.correctBoundaryConditions();

    correctNut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
