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
    Foam::LESModels::WrayAgarwalDESfdes

Group
    grpDESTurbulence

Description
    WrayAgarwalDESfdes DES turbulence model for incompressible and
    compressible flows

    Reference:
    \verbatim
        Spalart, P. R., Jou, W. H., Strelets, M., & Allmaras, S. R. (1997).
        Comments on the feasibility of LES for wings, and on a hybrid
        RANS/LES approach.
        Advances in DNS/LES, 1, 4-8.
    \endverbatim

SourceFiles
    WrayAgarwalDESfdes.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwalDESfdes_H
#define WrayAgarwalDESfdes_H

#include "LESeddyViscosity.H"

#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class WrayAgarwalDESfdes Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WrayAgarwalDESfdes
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        WrayAgarwalDESfdes(const WrayAgarwalDESfdes&);
        WrayAgarwalDESfdes& operator=(const WrayAgarwalDESfdes&);


protected:

    // Protected data

        // Model constants
            // RANS constants
            dimensionedScalar sigma_;
            dimensionedScalar kappa_;
            dimensionedScalar Aplus_;
            dimensionedScalar C1ke_;
			dimensionedScalar C1kw_;
            dimensionedScalar sigmake_;
            dimensionedScalar sigmakw_;
			dimensionedScalar C2ke_;
            dimensionedScalar C2kw_;
            // LES constants
            dimensionedScalar CDES_;

        // Fields
            // RANS fields
            volScalarField Rnu_;
            volScalarField fdes_;
			volScalarField Switch1_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;
            
            // LES fields
            volScalarField outDelta_;
            volScalarField blendfactor_;
            surfaceScalarField UBlendingFactor_;
            surfaceScalarField RnuBlendingFactor_;
            surfaceScalarField pBlendingFactor_;
            surfaceScalarField KBlendingFactor_;
            surfaceScalarField eBlendingFactor_;
            surfaceScalarField hBlendingFactor_;
            

    // Protected Member Functions
        // RANS functions
        tmp<volScalarField> chi() const;

        tmp<volScalarField> fv1(const volScalarField& chi) const;

        tmp<volScalarField> blend
        (
            const volScalarField& Switch,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const;
        
        tmp<volScalarField> sigma(const volScalarField& Switch) const
        {
            return blend(Switch, sigmakw_, sigmake_);
        }

        tmp<volScalarField> C1(const volScalarField& Switch) const
        {
            return blend(Switch, C1kw_, C1ke_);
        }

        void correctNut(const volScalarField& fv1);
        virtual void correctNut();
        
        // LES functions
        tmp<volScalarField> fdes
        (
            const volScalarField& S
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("WrayAgarwalDESfdes");


    // Constructors

        //- Construct from components
        WrayAgarwalDESfdes
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
    virtual ~WrayAgarwalDESfdes()
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

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WrayAgarwalDESfdes.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
