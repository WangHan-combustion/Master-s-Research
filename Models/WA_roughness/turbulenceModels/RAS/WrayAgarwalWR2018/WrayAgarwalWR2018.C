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

#include "WrayAgarwalWR2018.H"
#include "bound.H"
#include "wallFvPatch.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::chi() const
{
    //return Rnu_/this->nu();
   // return Rnu_/this->nu()+Cr1_*0.03*mag(this->U_)*ks_*dimensionedScalar("ks",dimensionSet(0,1,0,0,0),1.0)/this->nu();//+Cr1_*ks_/(y_+0.03*ks_);
    
    
        return (Rnu_/(this->nu()))/(0.41*0.03*ks_);//+Cr1_*0.03*mag(this->U_)*ks_*dimensionedScalar("ks",dimensionSet(0,1,0,0,0),1.0)/this->nu();

   
}




template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::fmu
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return (chi3/(chi3 + pow3(Cw_)))*(1.0+ks_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::Rt() const
{
    return Rnu_/this->nu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::WDF_R
(
    const volScalarField& S,
    const volScalarField& W   
) const
{
    return mag(W/S);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::WDF_omega
(
    const volScalarField& S 
) const
{
    return S/sqrt(Cmu_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::WDF_k
(
    const volScalarField& omega
) const
{
    return (this->nut_+this->nu())  * omega;//*(1.0/(1.0+mag(this->U_)*ks_*dimensionedScalar("ks",dimensionSet(0,1,0,0,0),1.0)/this->nu()));
	//return Rnu_*omega;

}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::arg1
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
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::calcSwitch
(
    const volScalarField& S,
    const volScalarField& W    
) const
{
    return tanh(pow(Cs1_*arg1(S, W), Cs2_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::blend
(
    const volScalarField& Switch,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
) const
{
    return Switch*(psi1 - psi2) + psi2;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::sigmaR(const volScalarField& Switch) const
{
    return blend(Switch, sigmakw_, sigmake_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::C1(const volScalarField& Switch) const
{
    return blend(Switch, C1kw_, C1ke_);
}         


template<class BasicTurbulenceModel>
void WrayAgarwalWR2018<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fmu
)
{
    this->nut_ = Rnu_*fmu;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void WrayAgarwalWR2018<BasicTurbulenceModel>::correctNut()
{
    correctNut(fmu(this->chi()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwalWR2018<BasicTurbulenceModel>::WrayAgarwalWR2018
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
    
    ks_
   /*
    (
        dimensioned<scalar>
        (
            "ks",
            dimensionSet(0,1,0,0,0),
            0.0015//0.00086655//0.00025
        )
    ),
   */
 
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ks",
            this->coeffDict_,
            0.00125
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
    
    Cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr1",
            this->coeffDict_,
            1.0
        )
    ),
    
    Cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr2",
            this->coeffDict_,
            0.5
        )
    ),
    
    Cr3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr3",
            this->coeffDict_,
            1.0
        )
    ),
    
    C1ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1ke",
            this->coeffDict_,
			0.1284
        )
    ),

    C1kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1kw",
            this->coeffDict_,
			0.0829
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

    Cs1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs1",
            this->coeffDict_,
			1.0
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

    Cm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cm",
            this->coeffDict_,
			8.0
        )
    ),

    Rnu_
    (
        IOobject
        (
            "Rnu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    
    dnew_
    (
        IOobject
        (
            "dnew",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 1, 0, 0, 0), 0.0)
    ),

    utau_
    (
        IOobject
        (
            "utau",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 1, -1, 0, 0), 0.0)
    ),
    
    wGUx_
    (
        IOobject
        (
            "wGUx",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),
    
    wallGradUTest_
    (
        IOobject
        (
            "wallGradUTest",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedVector
        (
            "wallGradU",
            U.dimensions()/dimLength,
            vector::zero
        )
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
    
    y_(wallDist::New(this->mesh_).y())
    
    
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WrayAgarwalWR2018<BasicTurbulenceModel>::read()
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
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::DRnuEff(volScalarField Switch) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRnuEff", (Rnu_*sigmaR(Switch) + (this->nu()))*(1.0+Cr1_*0.03*mag(this->U_)*ks_*dimensionedScalar("ks",dimensionSet(0,1,0,0,0),1.0)/this->nu()))
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::k() const
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
tmp<volScalarField> WrayAgarwalWR2018<BasicTurbulenceModel>::epsilon() const
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
void WrayAgarwalWR2018<BasicTurbulenceModel>::correct()
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
	W_ = W;

    // Calculate switch function (f1)
    f1_ = calcSwitch(S, W);
    /*f1_ = min(f1_, 0.9);
    bound(f1_,SMALL);*/

    const volScalarField Rt(this->Rt());
    const volScalarField chi(this->chi());
    const volScalarField fmu(this->fmu(chi));
    const volVectorField& U = this->U_;
    const volScalarField nut = this->nut_;
    const volScalarField nu = this->nu();
           
    // Apply wall roughness boundary condition


    volScalarField Ux=U.component(0);
    surfaceScalarField gradnU = fvc::snGrad(Ux);
    surfaceScalarField gradR = fvc::snGrad(Rnu_);
    dnew_ = y_ + 0.03 * ks_*dimensionedScalar("ks",dimensionSet(0,1,0,0,0),1.0);
    const fvPatchList& patches = this->mesh_.boundary();
    volScalarField Rnutemp = Rnu_;
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];
        if (isType<wallFvPatch>(currPatch))
        {       
            forAll(currPatch, facei)
            {
                wallGradUTest_.boundaryField()[patchi] = -this->U_.boundaryField()[patchi].snGrad();
             
                wGUx_.boundaryField()[patchi] = wallGradUTest_.boundaryField()[patchi].component(vector::X);
            
                utau_.boundaryField()[patchi] =sqrt((this->nu()+this->nut_)*mag(wGUx_.boundaryField()[patchi]));
                
                //label faceCelli = currPatch.faceCells()[facei];
                
                /*Rnu_[faceCelli] = (nut[facei] + nu[facei]) * Rnutemp[faceCelli] * 0.41 * (y_[faceCelli]+ks_.value());///fmu[faceCelli];*/
               /* Rnu_[faceCelli] = utau_[faceCelli] * 0.41 * (y_[faceCelli]+0.03*ks_.value());*/
               //gradR[faceCelli] = Rnutemp[faceCelli]/(0.03 * ks_.value());
          /*
                Rnu_[faceCelli] = (Cr3_.value()*
                                                pow(
                                                    mag(gradnU[faceCelli])*(1.0+Rt[faceCelli]*fmu[faceCelli]),
                                                                                              Cr2_.value())*
                                                                                                            (pow(ks_.value(),(2*Cr2_.value())))
          
                                                                                                                                             )/fmu[faceCelli];*/  
            }
        }
    }
    
    
        
        
    // Define and solve R-Equation
    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*DRnuEff(f1_), Rnu_)
     ==
        alpha*rho*C1(f1_)*fvm::Sp(S, Rnu_)
      + alpha*rho*f1_*C2kw_*Cr2_*fvm::Sp((fvc::grad(Rnu_)&fvc::grad(S))/S, Rnu_)
      - alpha*rho*(1.0-f1_)*min(C2ke_*Rnu_*Rnu_*magSqr(fvc::grad(S))/S2,
                                Cm_*magSqr(fvc::grad(Rnu_)))
    );
    RnuEqn().relax();
    solve(RnuEqn);
    bound(Rnu_, dimensionedScalar("0", Rnu_.dimensions(), 0.0));
    Rnu_.correctBoundaryConditions();
    correctNut();

    




    
//calc Utau
  	//volScalarField Rnutemp = Rnu_;
  
   //volScalarField ksPlus;
  /*
  	const fvPatchList& patches = this->mesh_.boundary();
    forAll(patches, patchi)
    {
       const fvPatch& currPatch = patches[patchi];
        if (isType<wallFvPatch>(currPatch))
        {
        
            wallGradUTest_.boundaryField()[patchi] =
             -this->U_.boundaryField()[patchi].snGrad();
            wGUx_.boundaryField()[patchi] = 
            wallGradUTest_.boundaryField()[patchi].component(vector::X);
        utau_.boundaryField()[patchi] =sqrt((this->nu()+this->nut_)*mag(wGUx_.boundaryField()[patchi]));
	// Rnu_.boundaryField()[patchi] = utau_.boundaryField()[patchi]*wGUx_.boundaryField()[patchi]/(0.41*0.03*ks_.value());
		}
    }
 */
 
    // Apply wall roughness boundary condition
 /* 
    surfaceScalarField gradnR = fvc::snGrad(Rnu_);

    //const fvPatchList& patches = this->mesh_.boundary();
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];
        if (isType<wallFvPatch>(currPatch))
        {
            
            forAll(currPatch, facei)
            {
            
                label faceCelli = currPatch.faceCells()[facei];
          
				Rnu_[faceCelli] = utau_[faceCelli]*min(ksPlus[faceCelli]/90.0,1.0)
									/
									(0.41*0.03*ks_.value());
			
				
           }
        }
    }
*/

                
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
