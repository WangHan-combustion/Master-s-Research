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

#include "kkloneTranMenter.H"
#include "bound.H"
#include "wallDist.H"
//#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/(this->nu());
}
//intermittency member functions*****************************************

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::Nt
(
   
) const
{
    return  this->rho_*this->nut_/this->mu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::Nv
(
  const volScalarField& S
) const
{
    return (this->rho_*sqr(y_)*S)/(this->mu());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::Tu
(
  const volScalarField& S
) const
{
  volScalarField omeg = max(S_/sqrt(Cmu_), dimensionedScalar("1e-15", dimensionSet(0, 0, -1, 0, 0), 1e-15));

    return max(min(scalar(100.0)*(sqrt(nuTilda_/scalar(1.5))/(sqrt(omeg)*y_)), scalar(100.0)), TuLim_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::lambda() const
{
   // volTensorField tG = fvc::grad(this->U_);
   // volScalarField tG =fvc::grad(this->U_);
   volScalarField tG = (fvc::grad((n_ & this->U_)/mag(n_)) & n_)/mag(n_);
    
    return min(max((scalar(-7.57e-3)*(tG)*sqr(y_)/this->nu())+scalar(0.0128), scalar(-1.0)), scalar(1.0));
}
 
template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::FPG
(
  const volScalarField& lambda
) const
{
   volScalarField FPG 
        (
            IOobject
            (
                "FPG",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
            gamma_ 
    );
    
   forAll(FPG, cellI)
    {
      if(lambda[cellI]  < scalar(0.0))
        FPG[cellI] = min(scalar(1.0)+(scalar(-7.34)*lambda[cellI])+(0.0*min(lambda[cellI]+scalar(0.0681), scalar(0.0))), scalar(3.0));
      else
        FPG[cellI] = min(scalar(1.0)+(14.68*lambda[cellI]),scalar(1.5));
    }  
    
    return max(FPG, scalar(0.0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::nuTildaThetac
(
  const volScalarField& FPG,
  const volScalarField& Tu
) const
{
    return (CTU1_ + scalar(1000.0)*exp(CTU3_*Tu*FPG));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::Fonsetorg
(
  const volScalarField& Nv,
  const volScalarField& nuTildaThetac,
  const volScalarField& Nt
) const
{
    return max(min((Nv/(2.2*nuTildaThetac)), scalar(2.0)) - max(1.0-pow3(Nt/scalar(3.5)), scalar(0.0)), scalar(0.0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::Fturb 
(
  const volScalarField& Nt
) const
{
    return exp(-(pow4(Nt/2.0)));
}

//intermittency member functions*****************************************
template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::ka
(
    const volScalarField& S
) const
{
    return nuTilda_*S/a1_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::len
(
    const volScalarField& S
) const
{
    return pow(nuTilda_,0.5)*pow(S,-0.5);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::lvkmin
(
    const volScalarField& len
) const
{
    return len/C11_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::ep
(
    const volScalarField& ka
) const
{
    return y_*sqrt(scalar(0.3)*ka)/(scalar(20.0)*(this->nu()));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::fphi
(
    const volScalarField& ep
) const
{
    return (1.0+Cd1_*ep)/(1.0+pow(ep,4.0));
}

template<class BasicTurbulenceModel>
void kkloneTranMenter<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = nuTilda_;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kkloneTranMenter<BasicTurbulenceModel>::kkloneTranMenter
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
    
    //for the transition model
    Flength_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Flength",
            this->coeffDict_,
            100.0
        )
    ),
    
    ca2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ca2",
            this->coeffDict_,
            0.06
        )
    ),
    
    ce2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce2",
            this->coeffDict_,
            50.0
        )
    ),

    sigmaf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaf",
            this->coeffDict_,
            1.0
        )
    ),
    
    CTU1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU1",
            this->coeffDict_,
            100.0
        )
    ),
    
    CTU3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU3",
            this->coeffDict_,
            -1.0
        )
    ),
    /*
    dimensionedScalar TuLim_;
    dimensionedScalar gLim1_;
    dimensionedScalar gLim2_;
    */
    NtLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gLim2",
            this->coeffDict_,
            3.5
        )
    ),            

    CP3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CP3",
            this->coeffDict_,
            1.0
        )
    ), 
    
    Ct1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct1",
            this->coeffDict_,
            1.0
        )
    ), 
    
    Ct2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct2",
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
       /*
    gamma_
    (
        IOobject
        (
            "gamma",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
*/
    gamma_
    (
        IOobject
        (
            "gamma",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    this->mesh_
    ),

    
    y_(wallDist::New(this->mesh_).y()),
    n_(wallDist::New(this->mesh_).n())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kkloneTranMenter<BasicTurbulenceModel>::read()
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
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::DgammaEff() const
{
    return tmp<volScalarField>
    (
        //new volScalarField("DgammaEff", (this->nut_/sigmaf_) + this->nu())
        new volScalarField("DgammaEff", (nuTilda_/sigmaf_) + this->nu())
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkloneTranMenter<BasicTurbulenceModel>::DnuTildaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DnuTildaEff", nuTilda_*pk_ + this->nu())
    );
}

template<class BasicTurbulenceModel>
void kkloneTranMenter<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    
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
    
//intermittency Equation*****************************************
    tmp<fvScalarMatrix> gammaEqn
    (
        fvm::ddt(alpha, rho, gamma_)
      + fvm::div(alphaRhoPhi, gamma_)
      - fvm::laplacian(alpha*rho*DgammaEff(), gamma_)
     ==
        alpha*rho*Flength_*Fonsetorg*S*gamma_
        - fvm::Sp(alpha*rho*Flength_*S*gamma_*Fonsetorg, gamma_)
        + alpha*rho*ca2_*Fturb*W*gamma_
        - fvm::Sp(alpha*rho*ca2_*W*Fturb*(ce2_*gamma_), gamma_)
    ); 

    gammaEqn().relax();
    solve(gammaEqn);

    bound(gamma_,scalar(0));
    gamma_ = min(gamma_,1.0);    //1.0
//intermittency Equation*****************************************
    
    
    
    
//nuTildaEqn Equation********************************************
    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, nuTilda_)
      + fvm::div(alphaRhoPhi, nuTilda_)
      - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
     ==
        alpha*rho*a1_*Ct1_*gamma_*(Cphi1_-0.5)*G/S
      + alpha*rho*(0.5*a1_-pow(a1_,-0.5)*Cphi2_)*gamma_*fvm::Sp(S, nuTilda_)
      + alpha*rho*(this->nu())*nuTilda_*(1-f1_*fphi)*gamma_/sqr(y_)
      + alpha*rho*pl_*0.5*fvm::Sp((fvc::grad(nuTilda_) & fvc::grad(S))/S, nuTilda_)
      + alpha*rho*ph_*0.75*magSqr(fvc::grad(nuTilda_))
      +CP3_*5*max(gamma_-0.2, 0.0)*(1.0-gamma_)*min(max((Nv/2420)-1.0, 0.0), 3.0)*max(3.0*this->nu()-nuTilda_, dimensionedScalar("0", dimensionSet(0, 2, -1, 0, 0), 0))*W
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
