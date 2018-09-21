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

#include "GReWAmm.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::chi() const
{
    //return Rnu_/(this->nu())+Cr0_*RH_*dimensionedScalar("1",dimensionSet(0,1,0,0,0),1)/(y_);
    return Rnu_/(this->nu());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::Rt() const
{
    return (this->rho_*Rnu_)/(this->mu());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::Rv
(
  const volScalarField& S
) const
{
    return (this->rho_*sqr(y_)*S)/(this->mu());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::Tu
(
  const volScalarField& S
) const
{
  //volTensorField tgradU = fvc::grad(this->U_);
  //volTensorField Symm = (tgradU + tgradU.T());
  //volScalarField Pb = Rnu_*mag(dev(twoSymm(tgradU)));
  //volScalarField UU = max(mag(this->U_), dimensionedScalar("1e-15", dimensionSet(0, 1, -1, 0, 0), 1e-15));
  //volScalarField dw = max(mag(y_), dimensionedScalar("1e-6", dimensionSet(0, 1, 0, 0, 0), 1e-15));
  volScalarField omeg = max(S_/sqrt(Cmu_), dimensionedScalar("1e-15", dimensionSet(0, 0, -1, 0, 0), 1e-15));

    //return max(min(scalar(100)*(sqrt(Rnu_/scalar(1.5))/(sqrt(omeg)*y_)+(Cr1_*sqrt(Rnu_/(this->nu()))*RH_*dimensionedScalar("1",dimensionSet(0,1,0,0,0),1)/y_)), scalar(100)), TuLim_);
    //return max(min(scalar(100)*sqrt((Rnu_/(scalar(1.5)*omeg))+(Rnu_*sqr(RH_*Cr1_)*dimensionedScalar("1",dimensionSet(0,2,0,0,0),1)/((this->nu()))))/y_, scalar(100)), TuLim_);
    return max(
                min(   
                     scalar(100)*(
                                  /*sqrt(
                                       (sqr(Cr1_)*sqr(1+Rnu_/(this->nu()))*sqr(RH_)*dimensionedScalar("1",dimensionSet(0,2,0,0,0),1)*(omeg))
                                            )+*/
                                             sqrt(Rnu_/scalar(1.5))
                                         )/
                                               (sqrt(omeg)*y_), 
                        scalar(100)), 
                                      TuLim_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::lambda() const
{
    //volTensorField tgradU = fvc::grad(this->U_);
    volScalarField tgradU = fvc::grad(n_ & this->U_) & n_;
    
    
    return min(max(((-7.57e-3)*tgradU*sqr(y_)/this->nu())+scalar(0.0128), scalar(-1.0)), scalar(1.0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::FPG
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
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::ReThetac
(
  const volScalarField& FPG,
  const volScalarField& Tu
) const
{
    return (scalar(100.0) + scalar(1000.0)*exp(-1.0*Tu*FPG));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::Fonset1
(
    const volScalarField& Rv,
    const volScalarField& ReThetac   
) const
{
    return (Rv/(2.2*ReThetac));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::Fonset2
(
    const volScalarField& Fonset1  
) const
{
    return min(Fonset1, scalar(2.0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::Fonset3
(
    const volScalarField& Rt  
) const
{
    return max(1.0-pow3(Rt/RtLim_), scalar(0.0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::Fonset
(
  const volScalarField& Fonset2,
  const volScalarField& Fonset3
) const
{
    return max(Fonset2-Fonset3, scalar(0.0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::Fturb 
(
  const volScalarField& Rt
) const
{
    return exp(-(pow4(Rt/2.0)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cw_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::WDF_R
(
    const volScalarField& S,
    const volScalarField& W   
) const
{
    return mag(W/S);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::WDF_omega
(
    const volScalarField& S 
) const
{
    return S/sqrt(Cmu_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::WDF_k
(
    const volScalarField& omegaa
) const
{
    return this->nut_*omegaa;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::arg1
(
    const volScalarField& S,
    const volScalarField& W
) const
{
    const volScalarField R = WDF_R(S, W);
    const volScalarField omegaa = WDF_omega(S);
    const volScalarField k = WDF_k(omegaa);

    const volScalarField eta = S*max(1.0, R);

    return (this->nu()+Rnu_)/2 * sqr(eta)/max(Cmu_*k*omegaa,
                                              dimensionedScalar("SMALL", 
                                                                dimensionSet(0, 2, -3, 0, 0), 
                                                                SMALL)
                                             );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::RnuWall
(
    const volScalarField& S
) const
{

    return (this->nut_*S/sqrt(Cmu_))*sqr(RH_)*dimensionedScalar("1",dimensionSet(0,2,0,0,0),1)
                /
            (Cr0_*(this->nu()));
}

/*template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::arg1
(
    const volScalarField& S,
    const volScalarField& W 
) const
{
    const volScalarField eta = y_*sqrt(Rnu_*S)/(20.0*this->nu());

    return (1.0+20.0*eta)/(1.0+sqr(max(y_*sqrt(Rnu_*S),1.5*Rnu_)/(20.0*this->nu())));
}*/

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::calcSwitch
(
    const volScalarField& S,
    const volScalarField& W    
) const
{
    return tanh(pow(arg1(S, W), 4.0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::blend
(
	const volScalarField& f1,
	const dimensionedScalar& psi1,
	const dimensionedScalar& psi2
) const
{
	return f1*(psi1 - psi2) + psi2;
}

template<class BasicTurbulenceModel>
void GReWAmm<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = Rnu_*fv1;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void GReWAmm<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
GReWAmm<BasicTurbulenceModel>::GReWAmm
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
    
    TuLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "TuLim",
            this->coeffDict_,
            0.0
        )
    ),
    
    gLim1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gLim1",
            this->coeffDict_,
            1.0
        )
    ),
    
    gLim2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gLim2",
            this->coeffDict_,
            1.0
        )
    ),
    
    RtLim_
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
    
    Clim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim",
            this->coeffDict_,
            8.0
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
    
    Cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw",
            this->coeffDict_,
            8.54
        )
    ),

    C1ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1ke",
            this->coeffDict_,
			0.1284//0.16 //0.1127
        )
    ),

    C1kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1kw",
            this->coeffDict_,
			0.0829 // 0.0829 doesnt change too much
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
    
    Cr0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr0",
            this->coeffDict_,
			1
        )
    ),
    
    /*Cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr1",
            this->coeffDict_,
			0.03
        )
    ),
    */
    RH_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "RH",
            this->coeffDict_,
			1e-3
        )
    ),

    C2ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2ke",
            this->coeffDict_,
			C1ke_.value()/sqr(kappa_.value()) + sigmake_.value()
        )
    ),

    C2kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2kw",
            this->coeffDict_,
			C1kw_.value()/sqr(kappa_.value()) + sigmakw_.value()
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
 
    dummy_
    (
        IOobject
        (
            "dummy",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("1.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    dummy2_
    (
        IOobject
        (
            "dummy2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("1.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    dummy_fonset1_
    (
        IOobject
        (
            "dummy_fonset1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("1.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    dummy_fonset2_
    (
        IOobject
        (
            "dummy_fonset2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("1.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    dummy_fonset3_
    (
        IOobject
        (
            "dummy_fonset3",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("1.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    dummy_fonset_
    (
        IOobject
        (
            "dummy_fonset",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("1.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    dummyTu_
    (
        IOobject
        (
            "dummyTu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("1.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
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

    /*SWdiff_
    (
        IOobject
        (
            "SWdiff",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),

    Cras_
    (
        IOobject
        (
            "Cras",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),*/
    
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
bool GReWAmm<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
        //sigma_.readIfPresent(this->coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::DgammaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DgammaEff", (this->nut_/sigmaf_) + this->nu())
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::DRnuEff(volScalarField f1) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRnuEff", Rnu_*sigma(f1) + this->nu())
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::k() const
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
tmp<volScalarField> GReWAmm<BasicTurbulenceModel>::epsilon() const
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
void GReWAmm<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    
	  // Calculate Strain rate magnitude S
	  volScalarField S2(2.0*magSqr(symm(fvc::grad(this->U_))));
	  volScalarField S = sqrt(S2);
	  volScalarField Sk(2.0*magSqr(skew(fvc::grad(this->U_))));
	  bound(S, dimensionedScalar("1e-15", S.dimensions(), SMALL)); // SMALL = 1e-15
	  bound(S2, dimensionedScalar("1e-15", S2.dimensions(), SMALL));
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
    const volScalarField Rt(this->Rt());
    const volScalarField Rv(this->Rv(S));
    const volScalarField Tu(this->Tu(S));
    const volScalarField lambda(this->lambda());
    const volScalarField FPG(this->FPG(lambda));
    const volScalarField ReThetac(this->ReThetac(FPG, Tu));
    const volScalarField Fonset1(this->Fonset1(Rv, ReThetac));
    const volScalarField Fonset2(this->Fonset2(Fonset1));
    const volScalarField Fonset3(this->Fonset3(Rt));
    const volScalarField Fturb(this->Fturb(Rt));
    const volScalarField Fonset(this->Fonset(Fonset2, Fonset3));
    const volScalarField RnuWall(this->RnuWall(S));
    
    dummy_ = lambda;
    dummy2_ = ReThetac;
    dummy_fonset1_ = Fonset1;
    dummy_fonset2_ = Fonset2;   
    dummy_fonset3_ = Fonset3;
    dummy_fonset_ = Fonset;
    dummyTu_ = Tu;
    
    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    // Eras
    /*Cras_ = 0;
    if (EnableEras_)
    {
        SWdiff_ = mag(S - W);
        forAll(Cras_, cellI)
        {
            if (SWdiff_[cellI] > 1e-14)
            {
                Cras_[cellI] = 9 * Cmu_.value(); 
            }
        }
    }*/
       
    // Intermittency equation
    
    tmp<fvScalarMatrix> gammaEqn
    (
        fvm::ddt(alpha, rho, gamma_)
      + fvm::div(alphaRhoPhi, gamma_)
      - fvm::laplacian(alpha*rho*DgammaEff(), gamma_)
     ==
        alpha*rho*Flength_*S*Fonset*gamma_
        - fvm::Sp(alpha*rho*Flength_*S*gamma_*Fonset, gamma_)
        + alpha*rho*ca2_*W*Fturb*gamma_
        - fvm::Sp(alpha*rho*ca2_*W*Fturb*(ce2_*gamma_), gamma_)
    ); 

    gammaEqn().relax();
    solve(gammaEqn);

    bound(gamma_,scalar(0));
    //gamma_ = min(gamma_,1.0);

    // Switch function (f1)
    f1_ = calcSwitch(S, W);
    //Switch1_ = max(min(Switch1_, 0.9), exp(-(pow((y_*sqrt(Rnu_*S/0.3)/this->nut_)/120.0, scalar(8)))));
    //f1_ = min(f1_, 0.9);
    //bound(f1_,SMALL);

    // R-Equation
    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*DRnuEff(f1_), Rnu_)
     ==
        C1(f1_)*alpha*rho*S*Rnu_*gamma_
	  + max(gamma_, gLim1_)*C2kw_*f1_*alpha*rho*fvm::Sp((fvc::grad(Rnu_) & fvc::grad(S))/S, Rnu_)
	  + CP3_*1.5*max(gamma_-0.2, 0.0)*(1.0-gamma_)*min(max((Rv/2420)-1.0, 0.0), 3.0)*max(3.0*this->nu()-this->nut_, dimensionedScalar("0", dimensionSet(0, 2, -1, 0, 0), 0))*W
	  - max(gamma_, gLim2_)*(1.0-f1_)*alpha*rho*min(C2ke_*sqr(Rnu_)*magSqr(fvc::grad(S))/S2, 
	                                          Clim_*magSqr(fvc::grad(Rnu_)))
    );

    RnuEqn().relax();
    solve(RnuEqn);
    bound(Rnu_, dimensionedScalar("1e-15", Rnu_.dimensions(), 1e-15));
    Rnu_.correctBoundaryConditions();

    correctNut();
    
    // Apply wall roughness boundary condition
    volScalarField RnuWalltemp = RnuWall;
    
    const fvPatchList& patches = this->mesh_.boundary();
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];
        if (isType<wallFvPatch>(currPatch))
        {       
            forAll(currPatch, facei)
            {
            
                label faceCelli = currPatch.faceCells()[facei];
          
                Rnu_[faceCelli] = RnuWalltemp[faceCelli];
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
