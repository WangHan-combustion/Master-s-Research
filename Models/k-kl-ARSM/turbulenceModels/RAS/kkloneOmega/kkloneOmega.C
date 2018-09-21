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

#include "kkloneOmega.H"
#include "bound.H"
#include "wallDist.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/(this->nu());
}
//kkloneOmega member functions*****************************************
template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::fv(const volScalarField& Ret) const
{
    return(1.0 - exp(-sqrt(Ret)/Av_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::fINT() const
{
    return
    (
        min
        (
            kt_/(Cint_*(kl_ + kt_ + this->kMin_)),
            dimensionedScalar("1.0", dimless, 1.0)
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::fSS(const volScalarField& Omega) const
{
    return(exp(-sqr(Css_*this->nu()*Omega/(kt_ + this->kMin_))));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::Cmu(const volScalarField& S) const
{
    return(1.0/(A0_ + As_*(S/(omega_ + this->omegaMin_))));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::BetaTS(const volScalarField& ReOmega) const
{
    return(scalar(1) - exp(-sqr(max(ReOmega - CtsCrit_, scalar(0)))/Ats_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& ktL,
    const volScalarField& Omega
) const
{
    return
    (
        scalar(1)
      - exp
        (
            -CtauL_*ktL
          /
            (
                sqr
                (
                    lambdaEff*Omega
                  + dimensionedScalar
                    (
                        "ROOTVSMALL",
                        dimLength*inv(dimTime),
                        ROOTVSMALL
                    )
                )
            )
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fv,
    const volScalarField& ktS
) const
{
    return(fv*CmuStd_*sqrt(ktS)*lambdaEff);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::fOmega
(
    const volScalarField& lambdaEff,
    const volScalarField& lambdaT
) const
{
    return
    (
        scalar(1)
      - exp
        (
           -0.41
           *pow4
            (
                lambdaEff
              / (
                    lambdaT
                  + dimensionedScalar
                    (
                        "ROTVSMALL",
                        lambdaT.dimensions(),
                        ROOTVSMALL
                    )
                )
            )
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::phiBP(const volScalarField& Omega) const
{
    return
    (
        min
        (
            max
            (
                kt_/this->nu()
             / (
                    Omega
                  + dimensionedScalar
                    (
                        "ROTVSMALL",
                        Omega.dimensions(),
                        ROOTVSMALL
                    )
                )
              - CbpCrit_,
                scalar(0)
            ),
            scalar(50.0)
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::phiNAT
(
    const volScalarField& ReOmega,
    const volScalarField& fNatCrit
) const
{
    return
    (
        max
        (
            ReOmega
          - CnatCrit_
          / (
                fNatCrit + dimensionedScalar("ROTVSMALL", dimless, ROOTVSMALL)
            ),
            scalar(0)
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::D(const volScalarField& k) const
{
    return this->nu()*magSqr(fvc::grad(sqrt(k)));
}


//kkloneOmega member functions*****************************************
template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::ka
(
    const volScalarField& S
) const
{
    return nuTilda_*S/a1_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::len
(
    const volScalarField& S
) const
{
    return pow(nuTilda_,0.5)*pow(S,-0.5);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::lvkmin
(
    const volScalarField& len
) const
{
    return len/C11_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::ep
(
    const volScalarField& ka
) const
{
    return y_*sqrt(scalar(0.3)*ka)/(scalar(20.0)*(this->nu()));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::fphi
(
    const volScalarField& ep
) const
{
    return (1.0+Cd1_*ep)/(1.0+pow(ep,4.0));
}

template<class BasicTurbulenceModel>
void kkloneOmega<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = nuTilda_;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kkloneOmega<BasicTurbulenceModel>::kkloneOmega
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
            0.97IOobject::groupName("omega", U.group())
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
    
    //for the omega model
     A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.04
        )
    ),
    As_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "As",
            this->coeffDict_,
            2.12
        )
    ),
    
    Av_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Av",
            this->coeffDict_,
            6.75
        )
    ),
    
    Abp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Abp",
            this->coeffDict_,
            0.6
        )
    ),
    Anat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anat",
            this->coeffDict_,
            200
        )
    ),
    Ats_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ats",
            this->coeffDict_,
            200
        )
    ),
    CbpCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CbpCrit",
            this->coeffDict_,
            1.2
        )
    ),
    Cnc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cnc",
            this->coeffDict_,
            0.1
        )
    ),
    CnatCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CnatCrit",
            this->coeffDict_,
            1250
        )
    ),
    Cint_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            this->coeffDict_,
            0.75
        )
    ),
    CtsCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtsCrit",
            this->coeffDict_,
            1000
        )
    ),
    CrNat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CrNat",
            this->coeffDict_,
            0.02
        )
    ),
    /*
    C111_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C111",
            this->coeffDict_,
            3.4e-6
        )
    ),
    C112_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C112",
            this->coeffDict_,
            1.0e-10IOobject::groupName("omega", U.group())
        )
    ),
    */
    CR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CR",
            this->coeffDict_,
            0.12
        )
    ),
    CalphaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CalphaTheta",
            this->coeffDict_,
            0.035
        )
    ),
    Css_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Css",
            this->coeffDict_,
            1.5
        )
    ),
    CtauL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtauL",
            this->coeffDict_,
            4360
        )
    ),
    Cw1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw1",
            this->coeffDict_,
            0.44
        )
    ),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.92
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            0.3
        )
    ),
    CwR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CwR",
            this->coeffDict_,
            1.5
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clambda",
           this->coeffDict_,
            2.495
        )
    ),
    CmuStd_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuStd",
            this->coeffDict_,
            0.09
        )
    ),
    Prtheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prtheta",
            this->coeffDict_,
            0.85
        )
    ),
    Sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmak",
            this->coeffDict_,
            1
        )
    ),
    Sigmaw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmaw",
            this->coeffDict_,
            1.17
        )
    ),
    kt_
    (
        IOobject
        (
           // IOobject::groupName("kt", U.group()),
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
            //IOobject::groupName("kl", U.group()),
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
            //IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    
    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_
        ),
        kt_*omega_ + D(kl_) + D(kt_)
    ),
    
    y_(wallDist::New(this->mesh_).y())
{
    bound(kt_, this->kMin_);
    bound(kl_, this->kMin_);
    bound(omega_, this->omegaMin_);
    bound(epsilon_, this->epsilonMin_);
    if (type == typeName)
    {
        this->nut_.correctBoundaryConditions();
        this->printCoeffs(type);
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kkloneOmega<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
    /*
        A0_.readIfPresent(coeffDict());
        As_.readIfPresent(coeffDict());
        Av_.readIfPresent(coeffDict());
        Abp_.readIfPresent(coeffDict());
        Anat_.readIfPresent(coeffDict());
        Abp_.readIfPresent(coeffDict());
        Ats_.readIfPresent(coeffDict());
        CbpCrit_.readIfPresent(coeffDict());
        Cnc_.readIfPresent(coeffDict());
        CnatCrit_.readIfPresent(coeffDict());
        Cint_.readIfPresent(coeffDict());
        CtsCrit_.readIfPresent(coeffDict());
        CrNat_.readIfPresent(coeffDict());
        C11_.readIfPresent(coeffDict());
        C12_.readIfPresent(coeffDict());
        CR_.readIfPresent(coeffDict());
        CalphaTheta_.readIfPresent(coeffDict());
        Css_.readIfPresent(coeffDict());
        CtauL_.readIfPresent(coeffDict());
        Cw1_.readIfPresent(coeffDict());
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        CwR_.readIfPresent(coeffDict());
        Clambda_.readIfPresent(coeffDict());
        CmuStd_.readIfPresent(coeffDict());
        Prtheta_.readIfPresent(coeffDict());
        Sigmak_.readIfPresent(coeffDict());
        Sigmaw_.readIfPresent(coeffDict());
        */
        return true;
    }
    else
    {
        return false;
    }
}



template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneOmega<BasicTurbulenceModel>::DnuTildaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DnuTildaEff", nuTilda_*pk_ + this->nu())
    );
}

template<class BasicTurbulenceModel>
void kkloneOmega<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    // New
    const volScalarField lambdaT(sqrt(kt_)/(omega_ + this->omegaMin_));

    const volScalarField lambdaEff(min(Clambda_*y_, lambdaT));

    const volScalarField fw
    (
        pow
        (
            lambdaEff
           /(lambdaT + dimensionedScalar("SMALL", dimLength, ROOTVSMALL)),
            2.0/3.0
        )
    );
    
    tmp<volTensorField> tU(fvc::grad(this->U_));
    const volTensorField& uGrad = tU();
    
    const volScalarField Omega(sqrt(2.0)*mag(skew(uGrad)));

    const volScalarField S22(2.0*magSqr(dev(symm(uGrad))));

    const volScalarField ktS(fSS(Omega)*fw*kt_);

    const volScalarField nuts
    (
        fv(sqr(fw)*kt_/this->nu()/(omega_ + this->omegaMin_))
       *fINT()
       *Cmu(sqrt(S22))*sqrt(ktS)*lambdaEff
    );
    const volScalarField Pkt(nuts*S22);

    const volScalarField ktL(kt_ - ktS);
    const volScalarField ReOmega(sqr(y_)*Omega/this->nu());
    const volScalarField nutl
    (
        min
        (
            C11_*fTaul(lambdaEff, ktL, Omega)*Omega*sqr(lambdaEff)
           *sqrt(ktL)*lambdaEff/this->nu()
          + C12_*BetaTS(ReOmega)*ReOmega*sqr(y_)*Omega
        ,
            0.5*(kl_ + ktL)/(sqrt(S22) + this->omegaMin_)
        )
    );

    const volScalarField Pkl(nutl*S22);

    const volScalarField alphaTEff
    (
        alphaT(lambdaEff, fv(sqr(fw)*kt_/this->nu()/(omega_ + this->omegaMin_)), ktS)
    );

    // By pass source term divided by kl_

    const dimensionedScalar fwMin("SMALL", dimless, ROOTVSMALL);

    const volScalarField Rbp
    (
        CR_*(1.0 - exp(-phiBP(Omega)()/Abp_))*omega_
       /(fw + fwMin)
    );

    const volScalarField fNatCrit(1.0 - exp(-Cnc_*sqrt(kl_)*y_/this->nu()));

    // Natural source term divided by kl_
    const volScalarField Rnat
    (
        CrNat_*(1.0 - exp(-phiNAT(ReOmega, fNatCrit)/Anat_))*Omega
    );
  
    // Calculate Strain rate magnitude S
    volScalarField S2(2.0*magSqr(symm(fvc::grad(this->U_))));
    volScalarField S = sqrt(S2);
    bound(S, dimensionedScalar("0", S.dimensions(), SMALL));
    bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));
    S_ = S;
    
	// Calculate vorticity magnitude W
    volScalarField W2(2.0*magSqr(skew(fvc::grad(this->U_))));
    volScalarField W = sqrt(W2);
	bound(W, dimensionedScalar("1e-15", W.dimensions(), SMALL));
	bound(W2, dimensionedScalar("1e-15", W2.dimensions(), SMALL));
	
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
    const volScalarField Nt(this->Nt());
    const volScalarField Nv(this->Nv(S));
    const volScalarField Tu(this->Tu(S));
    const volScalarField lambda(this->lambda());
    const volScalarField FPG(this->FPG(lambda));
    const volScalarField nuTildaThetac(this->nuTildaThetac(FPG, Tu));
    const volScalarField Fturb(this->Fturb(Nt));
    const volScalarField Fonsetorg(this->Fonsetorg(Nv, nuTildaThetac, Nt));
    

       
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
    
//omega Equation*****************************************
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(this->phi_, omega_)
      - fvm::laplacian(DomegaEff(alphaTEff), omega_)
     ==
        Cw1_*Pkt*omega_/(kt_ + this->kMin_)
      - fvm::SuSp
        (
            (1.0 - CwR_/(fw + fwMin))*kl_*(Rbp + Rnat)/(kt_ + this->kMin_)
          , omega_
        )
      - fvm::Sp(Cw2_*sqr(fw)*omega_, omega_)
      + (
            Cw3_*fOmega(lambdaEff, lambdaT)*alphaTEff*sqr(fw)*sqrt(kt_)
        )/pow3(y_)
    );

    omegaEqn().relax();
    solve(omegaEqn);
    bound(omega_, this->omegaMin_);
    omega_ = max (omega_, this->omegaMin_);
//omega Equation*****************************************
    
    
    
    
//nuTildaEqn Equation********************************************
    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, nuTilda_)
      + fvm::div(alphaRhoPhi, nuTilda_)
      - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
     ==
        alpha*rho*a1_*(Cphi1_-0.5)*G/S
      + alpha*rho*(0.5*a1_-pow(a1_,-0.5)*Cphi2_)*0.8*omega_*fvm::Sp(S, nuTilda_)
      + alpha*rho*(this->nu())*nuTilda_*(1-f1_*fphi)/sqr(y_)
      + alpha*rho*pl_*0.5*fvm::Sp((fvc::grad(nuTilda_) & fvc::grad(S))/S, nuTilda_)
      + alpha*rho*ph_*0.75*magSqr(fvc::grad(nuTilda_))
      - alpha*rho*pf_*0.25*E1e
    );
//nuTildaEqn Equation********************************************
    
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
