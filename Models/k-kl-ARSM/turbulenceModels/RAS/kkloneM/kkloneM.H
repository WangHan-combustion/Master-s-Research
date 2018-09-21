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
    Foam::RASModels::kkloneM

Group
    grpRASTurbulence

Description
    Wray-Agarwal one equation Turbulence Model 
    WA-2017 version on NASA Turbulence Modeling Resource (TMR) website
    
    Reference: 
    https://turbmodels.larc.nasa.gov/wray_agarwal.html

    The default model coefficients are
        kkloneMCoeffs
        {
            kappa       0.41
            Cw          8.54
            C1ke        0.1127
            C1kw        0.0829
            sigmake     1.0
            sigmakw     0.72
        }   
    
SourceFiles
    kkloneM.C

\*---------------------------------------------------------------------------*/

#ifndef kkloneM_H
#define kkloneM_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class kkloneM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kkloneM
:
    public eddyViscosity<RASModel<BasicTurbulenceModel> >
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kkloneM(const kkloneM&);
        kkloneM& operator=(const kkloneM&);


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
            

        // Fields

            volScalarField nuTilda_;
            volScalarField k_;
            volScalarField kl_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


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
    TypeName("kkloneM");


    // Constructors

        //- Construct from components
        kkloneM
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
    virtual ~kkloneM()
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
        
        //- Return the effective diffusivity for nuTilda
        virtual tmp<volScalarField> DnuTildaEff() const;
        virtual tmp<volScalarField> DklEff() const;
        virtual tmp<volScalarField> DkEff() const;
        
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
#   include "kkloneM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
