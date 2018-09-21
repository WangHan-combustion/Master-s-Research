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
    Foam::incompressible::RASModels::SpalartAllmaras

Group
    grpIcoRASTurbulence

SourceFiles
    SpalartAllmarasNoft2.C

\*---------------------------------------------------------------------------*/
#ifndef SpalartAllmarasNoft2_H
#define SpalartAllmarasNoft2_H

#include "RASModel.H"
#include "wallDist.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                    Class SpalartAllmarasNoft2 Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmarasNoft2
:
    public RASModel
{

protected:
    // Protected data
    
        // Model coefficients
            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;
            
            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;

        // Fields
        
            volScalarField nuTilda_;
            volScalarField nut_;    
            
            wallDist d_;
            
    // Protected Member Functions
    
        tmp<volScalarField> chi() const;

        tmp<volScalarField> fv1(const volScalarField& chi) const;

        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        
        tmp<volScalarField> fw(const volScalarField& Stilda) const;

public:

    //- Runtime type information
    TypeName("SpalartAllmarasNoft2");
    
    // Constructors
    
        //- Construct from components
        SpalartAllmarasNoft2
        (
            const volVectorField& U,
            const surfaceScalarField& phi,transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );
        
    //- Destructor
    virtual ~SpalartAllmarasNoft2()
    {}
    
    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }
        
    //- Return the effective diffusivity for nuTilda
    tmp<volScalarField> DnuTildaEff() const;
    
    //- Return the turbulence kinetic energy
    virtual tmp<volScalarField> k() const;

    //- Return the turbulence kinetic energy dissipation rate
    virtual tmp<volScalarField> epsilon() const;

    //- Return the Reynolds stress tensor
    virtual tmp<volSymmTensorField> R() const;

    //- Return the effective stress tensor including the laminar stress
    virtual tmp<volSymmTensorField> devReff() const;

    //- Return the source term for the momentum equation
    virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

    //- Return the source term for the momentum equation
    virtual tmp<fvVectorMatrix> divDevRhoReff
    (
        const volScalarField& rho,
        volVectorField& U
    )const;
    
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();

    //- Read RASProperties dictionary
    virtual bool read();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
