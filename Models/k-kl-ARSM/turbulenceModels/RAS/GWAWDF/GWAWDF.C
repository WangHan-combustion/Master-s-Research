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

#include "GWAWDF.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::chi() const
{
    return Rnu_/this->nu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::Rt() const
{
    return (this->nut_)/(this->nu());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::Rv
(
  const volScalarField& W
) const
{
    return (this->rho_*sqr(y_)*W)/(scalar(2.188)*this->mu());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::Tw
(
  const volScalarField& S,
  const volScalarField& W,
  const volScalarField& Rt
) const
{
  volScalarField omeg = max(S/sqrt(Cmu_), dimensionedScalar("1e-15", dimensionSet(0, 0, -1, 0, 0), 1e-15));
  
    return Rt*W/omeg;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::Rc
(
  const volScalarField& Tw
) const
{
    return scalar(400.0)-scalar(360.0)*min(Tw*0.5, 1.0);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::Rs
(
  const volScalarField& S
) const
{
  volScalarField omeg = max(S/sqrt(Cmu_), dimensionedScalar("1e-15", dimensionSet(0, 0, -1, 0, 0), 1e-15));
  volVectorField SS = fvc::grad(S);
  
    return y_*(n_&SS)*omeg/(sqrt(2.0)*sqr(S));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::gammaEff
(
    const volScalarField& Rv,
    const volScalarField& Rt,
    const volScalarField& Rs  
) const
{
    return max(gamma_, min(2.0, exp(-(pow3(Rt/10.0))))*max(Rv-200.0,0.0)*(min(1.0, max(10+(5.0*Rs), 0.0))*min(1.0, max(10-(5.0*Rs), 0.0))));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::Fgamma
(
    const volScalarField& Rv,
    const volScalarField& Rc  
) const
{
    return 2.0*max(0.0, min(100.0-(0.7*Rv), 1.0))*min(max(Rv-Rc,0.0),4.0);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::Ggamma
(
    const volScalarField& Rv
) const
{
    return 7.5*max(0.0, min(100.0-Rv, 1.0))*min(max(Rv-18,0.0),1.0);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::Fturb 
(
  const volScalarField& Rv,
  const volScalarField& Rt
) const
{
    return exp(-(pow(Rv*Rt, 1.2)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cw_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::WDF_R
(
    const volScalarField& S,
    const volScalarField& W   
) const
{
    return mag(W/S);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::WDF_omega
(
    const volScalarField& S 
) const
{
    return S/sqrt(Cmu_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::WDF_k
(
    const volScalarField& omegaa
) const
{
    //return max(this->nut_*omega, dimensionedScalar("SMALL", dimensionSet(0, 2, -2, 0, 0), SMALL));
    return this->nut_*omegaa;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::arg1
(
    const volScalarField& S,
    const volScalarField& W
) const
{
    const volScalarField R = WDF_R(S, W);
    const volScalarField omegaa = WDF_omega(S);
    const volScalarField k = WDF_k(omegaa);

    const volScalarField eta = S*max(1.0, R);

    // bound A and B may not bound A*B
    //return (this->nu()+Rnu_)/2 * sqr(eta)/(Cmu_*k*omega);
    return (this->nu()+Rnu_)/2 * sqr(eta)/max(Cmu_*k*omegaa,
                                              dimensionedScalar("SMALL", 
                                                                dimensionSet(0, 2, -3, 0, 0), 
                                                                SMALL)
                                             );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::calcSwitch
(
    const volScalarField& S,
    const volScalarField& W    
) const
{
    return tanh(pow(2.0*arg1(S, W), 4.0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::blend
(
	const volScalarField& f1,
	const dimensionedScalar& psi1,
	const dimensionedScalar& psi2
) const
{
	return f1*(psi1 - psi2) + psi2;
}

template<class BasicTurbulenceModel>
void GWAWDF<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = Rnu_*fv1;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void GWAWDF<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
GWAWDF<BasicTurbulenceModel>::GWAWDF
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
    
    
    sigmal_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaf",
            this->coeffDict_,
            5.0
        )
    ),
    
    sigmag_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaf",
            this->coeffDict_,
            0.2
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
            12.5
        )
    ),

    C1ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1ke",
            this->coeffDict_,
			0.338//0.16 //0.1127
        )
    ),

    C1kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1kw",
            this->coeffDict_,
			0.083 // 0.0829 doesnt change too much
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

    EnableEras_
    (
        Switch::lookupOrAddToDict
        (
            "EnableEras",
            this->coeffDict_,
            false
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
            "gammaEff",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("1.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    dummyTw_
    (
        IOobject
        (
            "dummyTw",
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
    
    W_
    (
        IOobject
        (
            "W",
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

    if (EnableEras_)
    {
        Info << "Eras term is enabled" << nl << nl;
    }
    else
    {
        Info << "Eras term is disabled" << nl << nl;
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool GWAWDF<BasicTurbulenceModel>::read()
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
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::DgammaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DgammaEff", (this->nut_/sigmag_) + (this->nu()/sigmal_))
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::DRnuEff(volScalarField f1) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRnuEff", Rnu_*sigma(f1) + this->nu())
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::k() const
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
tmp<volScalarField> GWAWDF<BasicTurbulenceModel>::epsilon() const
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
void GWAWDF<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    
	  // Calculate Strain rate magnitude S
	  volScalarField S2(2.0*magSqr(symm(fvc::grad(this->U_))));
	  volScalarField S = sqrt(S2);
	  volScalarField Sk(2.0*magSqr(skew(fvc::grad(this->U_))));
	  bound(S, dimensionedScalar("0", S.dimensions(), SMALL)); // SMALL = 1e-15
	  bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));
	  S_ = S;
	  
	  // Calculate vorticity magnitude W
    volScalarField W2(2.0*magSqr(skew(fvc::grad(this->U_))));
    volScalarField W = sqrt(W2);
	  bound(W, dimensionedScalar("0", W.dimensions(), SMALL));
	  bound(W2, dimensionedScalar("0", W2.dimensions(), SMALL));
	  W_=W;
	
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volScalarField Rt(this->Rt());
    const volScalarField Rv(this->Rv(W));
    const volScalarField Tw(this->Tw(S,W,Rt));
    const volScalarField Rc(this->Rc(Tw));
    const volScalarField Rs(this->Rs(S));
    const volScalarField gammaEff(this->gammaEff(Rv,Rt,Rs));
    const volScalarField Fgamma(this->Fgamma(Rv,Rc));
    const volScalarField Ggamma(this->Ggamma(Rv));
    const volScalarField Fturb(this->Fturb(Rv,Rt));
    
    dummy_ = Rv;
    dummy2_ = gammaEff;
    dummyTw_ = Tw;
    
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
        alpha*rho*Fgamma*W*sqrt(gamma_)*scalar(1.1)
        -fvm::Sp(alpha*rho*Fgamma*W*sqrt(gamma_), gamma_)
        -fvm::Sp(alpha*rho*Ggamma*Fturb*W*sqrt(gamma_), gamma_)
    ); 

    gammaEqn().relax();
    solve(gammaEqn);

    bound(gamma_,scalar(0));
    gamma_ = min(gamma_,1.0);

    // Switch function (f1)
    f1_ = calcSwitch(S, W);
    //Switch1_ = max(min(Switch1_, 0.9), exp(-(pow((y_*sqrt(Rnu_*S/0.3)/this->nut_)/120.0, scalar(8)))));
    //f1_ = min(f1_, 0.9);
    bound(f1_,SMALL);

    // R-Equation
    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*DRnuEff(f1_), Rnu_)
     ==
        C1(f1_)*alpha*rho*S*Rnu_*gammaEff
	  + C2kw_*f1_*rho*fvm::Sp((fvc::grad(Rnu_) & fvc::grad(S))/S, Rnu_)
	  - (1.0-f1_)*rho*(fvm::Sp((C2ke_)*Rnu_*magSqr(fvc::grad(S))/S2, Rnu_)
                            //+
                            //Cras_*magSqr(fvc::grad(Rnu_)) 
                           )
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
