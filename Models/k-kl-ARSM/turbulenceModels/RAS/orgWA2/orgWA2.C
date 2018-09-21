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

#include "orgWA2.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> orgWA2<BasicTurbulenceModel>::chi() const
{
    return Rnu_/this->nu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> orgWA2<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cw_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> orgWA2<BasicTurbulenceModel>::sigmaR
(
    const volScalarField& Switch
) const
{
    return Switch*(sigmakw_-sigmake_) + sigmake_;
}

template<class BasicTurbulenceModel>
void orgWA2<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = Rnu_*fv1;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void orgWA2<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
orgWA2<BasicTurbulenceModel>::orgWA2
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
            13.0
        )
    ),

    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            0.144
        )
    ),

    Cb_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb",
            this->coeffDict_,
            1.66
        )
    ),

    C2ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2ke",
            this->coeffDict_,
			1.86
        )
    ),

    C2kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2kw",
            this->coeffDict_,
			1.36 // 1.0
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
			 2.0 //0.5
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
    ),

    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool orgWA2<BasicTurbulenceModel>::read()
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
tmp<volScalarField> orgWA2<BasicTurbulenceModel>::DRnuEff(const volScalarField& Switch) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRnuEff", Rnu_*sigmaR(Switch) + this->nu())
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> orgWA2<BasicTurbulenceModel>::k() const
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
tmp<volScalarField> orgWA2<BasicTurbulenceModel>::epsilon() const
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
void orgWA2<BasicTurbulenceModel>::correct()
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

	volScalarField S2(2.0*magSqr(symm(fvc::grad(this->U_))));
	volScalarField S = sqrt(S2);
	bound(S, dimensionedScalar("0", S.dimensions(), SMALL));
	bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));

	// output S_
	S_ = S;

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

//////////////////////////////////////////////////////////////
    const volScalarField arg1 = min(
                                    Cb_*Rnu_/(S*sqr(kappa_*y_)),
                                    sqr((Rnu_+this->nu())/this->nu())
                                    //sqr(Rnu_/y_*dimensionedScalar("Min", dimensionSet(0, -1, 1, 0, 0), 1.0))
                                   );
    f1_ = tanh(pow(arg1, 4.0));

////////////////////////////////////////////////////////////////

    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*DRnuEff(f1_), Rnu_)
     ==
        C1_*alpha*rho*S*Rnu_
	  + f1_*C2kw_*rho*fvm::Sp((fvc::grad(Rnu_) & fvc::grad(S))/S, Rnu_)
	  - (1.0-f1_)*C2ke_*rho*fvm::Sp(Rnu_*magSqr(fvc::grad(S))/S2, Rnu_)
    );

    RnuEqn().relax();
    solve(RnuEqn);
    bound(Rnu_, dimensionedScalar("0", Rnu_.dimensions(), 0.0));
    Rnu_.correctBoundaryConditions();

    correctNut(fv1);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
