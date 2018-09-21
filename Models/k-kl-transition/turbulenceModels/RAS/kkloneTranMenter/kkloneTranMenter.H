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
    Foam::RASModels::kkloneTranMenter

Group
    grpRASTurbulence

Description
    Wray-Agarwal one equation Turbulence Model 
    WA-2017 version on NASA Turbulence Modeling Resource (TMR) website
    
    Reference: 
    https://turbmodels.larc.nasa.gov/wray_agarwal.html

    The default model coefficients are
        kkloneCoeffs
        {
            kappa       0.41
            Cw          8.54
            C1ke        0.1127
            C1kw        0.0829
            sigmake     1.0
            sigmakw     0.72
        }   
    
SourceFiles
    kkloneTranMenter.C

\*---------------------------------------------------------------------------*/

#ifndef kkloneTranMenter_H
#define kkloneTranMenter_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class kkloneTranMenter Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kkloneTranMenter
:
    public eddyViscosity<RASModel<BasicTurbulenceModel> >
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kkloneTranMenter(const kkloneTranMenter&);
        kkloneTranMenter& operator=(const kkloneTranMenter&);


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
    
            dimensionedScalar Flength_;
            dimensionedScalar ca2_;
            dimensionedScalar ce2_;
            dimensionedScalar sigmaf_;
            dimensionedScalar CTU1_;
            dimensionedScalar CTU3_;
            dimensionedScalar TuLim_;
            dimensionedScalar gLim1_;
            dimensionedScalar gLim2_;
            dimensionedScalar NtLim_;
            dimensionedScalar CP3_;  
            dimensionedScalar Ct1_;   
            dimensionedScalar Ct2_;
                      

        // Fields
            volScalarField nuTilda_;
            volScalarField S_;
            volScalarField gamma_;


            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            
            const volScalarField& y_;
            const volVectorField& n_;  

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
    TypeName("kkloneTranMenter");


    // Constructors

        //- Construct from components
        kkloneTranMenter
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
    virtual ~kkloneTranMenter()
    {}


    // Member Functions
    
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
#   include "kkloneTranMenter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
