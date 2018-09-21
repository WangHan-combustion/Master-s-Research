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

#include "kkl.H"
#include "bound.H"
#include "wallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
void kkl<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = pow(Cmu_,0.25)*kl_/sqrt(k_);
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kkl<BasicTurbulenceModel>::kkl
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

    eta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "eta1",
            this->coeffDict_,
            1.2         
        )
    ),
    
    eta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "eta2",
            this->coeffDict_,
            0.97
        )
    ),
    
    eta3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "eta3",
            this->coeffDict_,
            0.13
        )
    ),
    
    Cphi2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cphi2",
            this->coeffDict_,
            eta3_.value()
        )
    ),
    
    C11_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C11",
            this->coeffDict_,
            10.0
        )
    ),
    
    C12_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C12",
            this->coeffDict_,
            1.3
        )
    ),

    Cd1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd1",
            this->coeffDict_,
            4.7
        )
    ),
    
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    
    sigmaphi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaphi",
            this->coeffDict_,
            1.0
        )
    ),
           
    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    
    kl_
    (
        IOobject
        (
            "kl",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    
    
    
/*   nuTilda_
    (
        IOobject
        (
            "nuTilda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),*/
    
    y_(wallDist::New(this->mesh_).y())
    
{
    
    if (type == typeName)
    {
        this->printCoeffs(type);
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kkl<BasicTurbulenceModel>::read()
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
tmp<volScalarField> kkl<BasicTurbulenceModel>::DklEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DklEff", this->nut_*sigmaphi_ + this->nu())
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkl<BasicTurbulenceModel>::DkEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DkEff", this->nut_*sigmak_ + this->nu())
    );
}

template<class BasicTurbulenceModel>
void kkl<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }  
    volScalarField S2(2.0*magSqr(symm(fvc::grad(this->U_))));
    volScalarField S = sqrt(S2);
    bound(S, dimensionedScalar("0", S.dimensions(), SMALL));
    bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));
    
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
       
    volScalarField& nut_ = this->nut_;

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();
    
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));
        
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField GP(this->GName(), nut_*(dev(twoSymm(tgradU())) && tgradU()));
    tgradU.clear();
    
    //limiter on G
    volScalarField G = min(GP-(2.0/3.0)*divU*k_,scalar(20.0)*nut_*S2);    
                                                                                                                        
    volScalarField fp = max(min(G/(pow(Cmu_,0.75)*pow(k_,2.5)/kl_),0.5),1.0);
    volScalarField lvkmax = C12_*kappa_*y_*fp;
    volScalarField lvkmin = kl_/(k_*C11_);
       
    volScalarField Ux=U.component(0); 
    volScalarField Uy=U.component(1);
    volScalarField Uz=U.component(2);
    volScalarField U2=sqrt(
                            sqr(fvc::laplacian(Ux))
                           +sqr(fvc::laplacian(Uy))
                           +sqr(fvc::laplacian(Uz))
                          );
                 
    volScalarField lvkl = kappa_*mag(S/U2);
    //limiter on lvk
    volScalarField lvk = max(lvkmin,min(lvkmax,lvkl));
    
    volScalarField Cphi1_= eta1_-eta2_*sqr(kl_/(k_*lvk));
    volScalarField ep=y_*sqrt(scalar(0.3)*k_)/(scalar(20.0)*(this->nu()));
    volScalarField fphi=(1+Cd1_*ep)/(1+pow(ep,4.0));
        
   
    
    //kl equation
    tmp<fvScalarMatrix> klEqn
    (
        fvm::ddt(alpha, rho, kl_)
      + fvm::div(alphaRhoPhi, kl_)
      - fvm::laplacian(alpha*rho*DklEff(), kl_)
     ==
        alpha*rho*Cphi1_*G*kl_/k_
      - alpha*rho*fvm::SuSp((2.0/3.0)*divU, kl_)
      - alpha*rho*Cphi2_*pow(k_,1.5)
      - alpha*rho*6*(this->nu())*kl_*fphi/sqr(y_)
    );
    
    klEqn().relax();
    solve(klEqn);
    bound(kl_, dimensionedScalar("0", kl_.dimensions(), 1e-15));
    kl_.correctBoundaryConditions();
    
    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - alpha*rho*fvm::SuSp((2.0/3.0)*divU, k_)
      - alpha*rho*pow(Cmu_,0.75)*pow(k_,2.5)/kl_
      - alpha*rho*2*(this->nu())*k_/sqr(y_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, dimensionedScalar("0", k_.dimensions(), 1e-15));
    k_.correctBoundaryConditions();

    correctNut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
