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

Class
    Foam::RASModels::kkloneOmega

Group
    grpRASTurbulence

Description
    Wray-Agarwal one equation Turbulence Model 
    WA-2017 version on NASA Turbulence Modeling Resource (TMR) website
    
    Reference: 
    https://turbmodels.larc.nasa.gov/wray_agarwal.html

    The default model coefficients are
       
        kkloneOmegaCoeffs
        {
            A0             4.04
            As             2.12
            Av             6.75
            Abp            0.6
            Anat           200
            Ats            200
            CbpCrit        1.2
            Cnc            0.1
            CnatCrit       1250
            Cint           0.75
            CtsCrit        1000
            CrNat          0.02
            C11            3.4e-6
            C12            1.0e-10
            CR             0.12
            CalphaTheta    0.035
            Css            1.5
            CtauL          4360
            Cw1            0.44
            Cw2            0.92
            Cw3            0.3
            CwR            1.5
            Clambda        2.495
            CmuStd         0.09
            Prtheta        0.85
            Sigmak         1
            Sigmaw         1.17
        
        }   
    
SourceFiles
    kkloneOmega.C

\*---------------------------------------------------------------------------*/

#ifndef kkloneOmega_H
#define kkloneOmega_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class kkloneOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kkloneOmega
:
    public eddyViscosity<RASModel<BasicTurbulenceModel> >
{
    // Private Member Functions
        
        tmp<volScalarField> fv(const volScalarField& Ret) const;

        tmp<volScalarField> fINT() const;

        tmp<volScalarField> fSS(const volScalarField& omega) const;

        tmp<volScalarField> Cmu(const volScalarField& S) const;

        tmp<volScalarField> BetaTS(const volScalarField& Rew) const;

        tmp<volScalarField> fTaul
        (
            const volScalarField& lambdaEff,
            const volScalarField& ktL,
            const volScalarField& omega
        ) const;

        tmp<volScalarField> alphaT
        (
            const volScalarField& lambdaEff,
            const volScalarField& fv,
            const volScalarField& ktS
        ) const;

        tmp<volScalarField> fOmega
        (
            const volScalarField& lambdaEff,
            const volScalarField& lambdaT
        ) const;

        tmp<volScalarField> phiBP(const volScalarField& omega) const;

        tmp<volScalarField> phiNAT
        (
            const volScalarField& ReOmega,
            const volScalarField& fNatCrit
        ) const;

        tmp<volScalarField> D(const volScalarField& k) const;
        // Disallow default bitwise copy construct and assignment
        kkloneOmega(const kkloneOmega&);
        kkloneOmega& operator=(const kkloneOmega&);


protected:

    // Protected data

        // Model coefficients
            dimensionedScalar kappa_;
            dimensionedScalar Aplus_;
            dimensionedScalar eta1_;
            dimensionedScalar eta2_;
            dimensionedScalar eta3_;
            dimensionedScalar Cphi2_;
            dimensionedScalar f1_;
            dimensionedScalar C11_;
            dimensionedScalar C12_;
            dimensionedScalar Cd1_;
            dimensionedScalar C3_;
            dimensionedScalar Cmu_;
            dimensionedScalar a1_;
            dimensionedScalar pl_;
            dimensionedScalar pk_;
            dimensionedScalar ph_;
            dimensionedScalar pf_;
            
        // gamma-ReThetat-model coefficients
    
            dimensionedScalar A0_;
            dimensionedScalar As_;
            dimensionedScalar Av_;
            dimensionedScalar Abp_;
            dimensionedScalar Anat_;
            dimensionedScalar Ats_;
            dimensionedScalar CbpCrit_;
            dimensionedScalar Cnc_;
            dimensionedScalar CnatCrit_;
            dimensionedScalar Cint_;
            dimensionedScalar CtsCrit_;
            dimensionedScalar CrNat_;
            //dimensionedScalar C111_;
            //dimensionedScalar C112_;
            dimensionedScalar CR_;
            dimensionedScalar CalphaTheta_;
            dimensionedScalar Css_;
            dimensionedScalar CtauL_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar CwR_;
            dimensionedScalar Clambda_;
            dimensionedScalar CmuStd_;
            dimensionedScalar Prtheta_;
            dimensionedScalar Sigmak_;
            dimensionedScalar Sigmaw_;
           

        // Fields
            volScalarField nuTilda_;
            volScalarField S_;
            volScalarField kt_;
            volScalarField kl_;
            volScalarField omega_;
            volScalarField epsilon_;


            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            
            const volScalarField& y_;


// Model functions of gamma-ReThetat-model with Crossflow
        tmp<volScalarField> Nt() const;
        tmp<volScalarField> Nv(const volScalarField& S) const;
        tmp<volScalarField> Tu(const volScalarField& S) const;
        tmp<volScalarField> lambda() const;       
        tmp<volScalarField> FPG(const volScalarField& lambda) const;
        tmp<volScalarField> nuTildaThetac(const volScalarField& FPG, const volScalarField& Tu) const;
        tmp<volScalarField> Fonsetorg(const volScalarField& Nv, const volScalarField& nuTildaThetac, const volScalarField& Nt) const;
        tmp<volScalarField> Fturb(const volScalarField& Nt) const;

    // Protected Member Functions

        tmp<volScalarField> chi() const;
        
        tmp<volScalarField> ka(const volScalarField& S) const;

        tmp<volScalarField> len(const volScalarField& S) const;
        
        tmp<volScalarField> lvkmin(const volScalarField& len) const;
        
        tmp<volScalarField> ep(const volScalarField& ka) const;
        
        tmp<volScalarField> fphi(const volScalarField& ep) const;
        

        

        virtual void correctNut();

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kkloneOmega");


    // Constructors

        //- Construct from components
        kkloneOmega
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kkloneOmega()
    {}


    // Member Functions


        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaT/Sigmaw_ + this->nu())
            );
        }

        //- Return the laminar kinetic energy
        virtual tmp<volScalarField> kl() const
        {
            return kl_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> kt() const
        {
            return kt_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }
        
        //- Return the total fluctuation kinetic energy
        virtual tmp<volScalarField> kf() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "kf",
                       this-> mesh_.time().timeName(),
                        this->mesh_
                    ),
                    kt_ + kl_,
                    omega_.boundaryField().types()
                )
            );
        }
    
        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
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

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
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
        
        //- Return the effective diffusivity (gamma eqn)
        tmp<volScalarField> DgammaEff() const;
        
        //- Return the effective diffusivity for nuTilda
        virtual tmp<volScalarField> DnuTildaEff() const;
        
        //- Re-read model coefficients if they have changed
        virtual bool read();
   
        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kkloneOmega.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
