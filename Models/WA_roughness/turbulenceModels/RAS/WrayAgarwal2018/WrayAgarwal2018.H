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
    Foam::RASModels::WrayAgarwal2018

Group
    grpRASTurbulence

Description
	Wray-Agarwal Turbulence Model with Rahman's wall distance free modification
with type-2 switch function, and Han's coefficient adjustment. 

SourceFiles
    WrayAgarwal2018.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwal2018_H
#define WrayAgarwal2018_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class WrayAgarwal2018 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WrayAgarwal2018
:
    public eddyViscosity<RASModel<BasicTurbulenceModel> >
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        WrayAgarwal2018(const WrayAgarwal2018&);
        WrayAgarwal2018& operator=(const WrayAgarwal2018&);

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar kappa_;
            dimensionedScalar ks_;
            dimensionedScalar Cw_;
            dimensionedScalar Cr1_;
            dimensionedScalar Cr2_;
            dimensionedScalar Cr3_;
            dimensionedScalar C1ke_;
			dimensionedScalar C1kw_;
            dimensionedScalar sigmake_;
            dimensionedScalar sigmakw_;
			dimensionedScalar C2ke_;
            dimensionedScalar C2kw_;
            dimensionedScalar Cmu_;
            dimensionedScalar Cs1_;
            dimensionedScalar Cs2_;
            dimensionedScalar Cm_;
        
        // Fields

            volScalarField Rnu_;
            volScalarField utau_;
            volScalarField wGUx_;
            volVectorField wallGradUTest_;
			volScalarField f1_;
            volScalarField S_;
            volScalarField W_;
           

    // Protected Member Functions

        tmp<volScalarField> chi() const;

        tmp<volScalarField> fmu(const volScalarField& chi) const;

        tmp<volScalarField> WDF_R
        (
            const volScalarField& S,
            const volScalarField& W    
        ) const;

        tmp<volScalarField> WDF_omega
        (
            const volScalarField& S 
        ) const;

        tmp<volScalarField> WDF_k
        (
            const volScalarField& omega 
        ) const;

        tmp<volScalarField> arg1
        (
            const volScalarField& S,
            const volScalarField& W    
        ) const;

        tmp<volScalarField> calcSwitch
        (
            const volScalarField& S,
            const volScalarField& W    
        ) const;

        tmp<volScalarField> blend
        (
            const volScalarField& Switch,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const;

        tmp<volScalarField> sigmaR(const volScalarField& Switch) const;

        tmp<volScalarField> C1(const volScalarField& Switch) const;

        void correctNut(const volScalarField& fmu);
        virtual void correctNut();

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("WrayAgarwal2018");

    // Constructors

        //- Construct from components
        WrayAgarwal2018
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
    virtual ~WrayAgarwal2018()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for Rnu
        tmp<volScalarField> DRnuEff(volScalarField Switch) const;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WrayAgarwal2018.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
