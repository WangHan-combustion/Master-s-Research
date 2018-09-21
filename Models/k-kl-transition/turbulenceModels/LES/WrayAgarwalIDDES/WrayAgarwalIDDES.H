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
    Foam::LESModels::WrayAgarwalIDDES

Group
    grpDESTurbulence

Description
    WrayAgarwalIDDES DES turbulence model for incompressible and
    compressible flows

    Reference:
    \verbatim
        Spalart, P. R., Jou, W. H., Strelets, M., & Allmaras, S. R. (1997).
        Comments on the feasibility of LES for wings, and on a hybrid
        RANS/LES approach.
        Advances in DNS/LES, 1, 4-8.
    \endverbatim

SourceFiles
    WrayAgarwalIDDES.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwalIDDES_H
#define WrayAgarwalIDDES_H

#include "WrayAgarwalDES.H"
#include "IDDESDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class WrayAgarwalIDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WrayAgarwalIDDES
:
    public WrayAgarwalDES<BasicTurbulenceModel>
{
    // Private data

        // Model constants

            dimensionedScalar fwStar_;
            dimensionedScalar cl_;
            dimensionedScalar ct_;
            dimensionedScalar Cd1_;
			Switch IDES_;

        // Fields
			
			volScalarField fdtilda_;
			volScalarField fe_;
			volScalarField fd_;
            const IDDESDelta& IDDESDelta_;


    // Private Member Functions

        tmp<volScalarField> alpha() const;
        tmp<volScalarField> ft(const volScalarField& magGradU) const;
        tmp<volScalarField> fl(const volScalarField& magGradU) const;

        tmp<volScalarField> rd
        (
            const volScalarField& nur,
            const volScalarField& magGradU
        ) const;

        //- Delay function
        tmp<volScalarField> fd(const volScalarField& magGradU) const;

        // Disallow default bitwise copy construct and assignment
        WrayAgarwalIDDES(const WrayAgarwalIDDES&);
        WrayAgarwalIDDES& operator=(const WrayAgarwalIDDES&);


protected:

    // Protected Member Functions

        //- Length scale

		virtual void precalculations
		(
			const volScalarField& S,
			const volTensorField& gradU
		);

        virtual tmp<volScalarField> fdes
        (
			const volScalarField& S,
            const volTensorField& gradU
        ) const;

		//- Blend factor in div scheme
		virtual void blendfactor();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("WrayAgarwalIDDES");


    // Constructors

        //- Construct from components
        WrayAgarwalIDDES
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
    virtual ~WrayAgarwalIDDES()
    {}


    // Member Functions

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WrayAgarwalIDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
