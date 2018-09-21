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

#include "kklone.H"
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
tmp<volScalarField> kklone<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/(this->nu());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kklone<BasicTurbulenceModel>::ka
(
    const volScalarField& S
) const
{
    return nuTilda_*S/a1_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kklone<BasicTurbulenceModel>::len
(
    const volScalarField& S
) const
{
    return pow(nuTilda_,0.5)*pow(S,-0.5);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kklone<BasicTurbulenceModel>::lvkmin
(
    const volScalarField& len
) const
{
    return len/C11_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kklone<BasicTurbulenceModel>::ep
(
    const volScalarField& ka
) const
{
    return y_*sqrt(scalar(0.3)*ka)/(scalar(20.0)*(this->nu()));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kklone<BasicTurbulenceModel>::fphi
(
    const volScalarField& ep
) const
{
    return (1+Cd1_*ep)/(1+pow(ep,4.0));
}

template<class BasicTurbulenceModel>
void kklone<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = nuTilda_;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kklone<BasicTurbulenceModel>::kklone
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
    
    f1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "f1",
            this->coeffDict_,
            6.0
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
    
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            7.0
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
    
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.3
        )
    ),
    
    pl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "pl",
            this->coeffDict_,
            1.0
        )
    ),
    
    pk_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "pk",
            this->coeffDict_,
            1.0
        )
    ),
    
    ph_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ph",
            this->coeffDict_,
            1.0
        )
    ),
    
    pf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "pf",
            this->coeffDict_,
            1.0
        )
    ),

    nuTilda_
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
bool kklone<BasicTurbulenceModel>::read()
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
tmp<volScalarField> kklone<BasicTurbulenceModel>::DnuTildaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DnuTildaEff", nuTilda_*pk_ + this->nu())
    );
}

template<class BasicTurbulenceModel>
void kklone<BasicTurbulenceModel>::correct()
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
    const volScalarField ka(this->ka(S));
    const volScalarField len(this->len(S));
    const volScalarField lvkmin(this->lvkmin(len));
    const volScalarField ep(this->ep(ka));
    const volScalarField fphi(this->fphi(ep));
       
    volScalarField& nut = this->nut_;

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();
    
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));
        
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField GP(this->GName(), nut*(dev(twoSymm(tgradU())) && tgradU()));
    tgradU.clear();
    
    //limiter on G
    volScalarField G = min(GP-(2.0/3.0)*divU*ka,scalar(20.0)*nuTilda_*S2);
    
    volScalarField fp = min(max(G/(max(
                                        nuTilda_*S2,
                                                    dimensionedScalar("SMALL", dimensionSet(0, 2, -3, 0, 0), 1e-20))),
                                                                                                                        0.5),1.0);
    volScalarField lvkmax = C12_*kappa_*y_*fp;
       
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
    volScalarField Cphi1_= eta1_-eta2_*sqr(len/lvk);
    
    //E1e
    volScalarField Eke = sqr(nuTilda_) * magSqr(fvc::grad(S)) / S2;
    volScalarField Ebb = max(magSqr(fvc::grad(nuTilda_)),
                         dimensionedScalar("EbbMin", dimensionSet(0, 2, -2, 0, 0), 1e-15));
    volScalarField E1e = C3_ * Ebb * tanh(Eke/(C3_*Ebb));
        

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, nuTilda_)
      + fvm::div(alphaRhoPhi, nuTilda_)
      - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
     ==
        alpha*rho*a1_*(Cphi1_-0.5)*G/S
      + alpha*rho*(0.5*a1_-pow(a1_,-0.5)*Cphi2_)*fvm::Sp(S, nuTilda_)
      + alpha*rho*(this->nu())*nuTilda_*(1-f1_*fphi)/sqr(y_)
      + alpha*rho*pl_*0.5*fvm::Sp((fvc::grad(nuTilda_) & fvc::grad(S))/S, nuTilda_)
      + alpha*rho*ph_*0.75*magSqr(fvc::grad(nuTilda_))
      - alpha*rho*pf_*0.25*E1e
    );
    
    nuTildaEqn().relax();
    solve(nuTildaEqn);
    bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

    correctNut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
