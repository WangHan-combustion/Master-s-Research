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

#include "WrayAgarwalDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDDES<BasicTurbulenceModel>::rd
(
    const volScalarField& nur,
    const volScalarField& magGradU
) const
{
    return min
    (
        nur
       /(
           max
           (
               magGradU,
               dimensionedScalar("SMALL", magGradU.dimensions(), SMALL)
           )*sqr(this->kappa_*this->y_)
       ),
       scalar(10)
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    return 1 - tanh(pow3(Cd1_*rd(this->nuEff(), magGradU)));
}

template<class BasicTurbulenceModel>
void WrayAgarwalDDES<BasicTurbulenceModel>::precalculations
(
	const volScalarField& S,
    const volTensorField& gradU
)
{
    const volScalarField magGradU(mag(gradU));
	
	fd_ = fd(magGradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalDDES<BasicTurbulenceModel>::fdes
(
	const volScalarField& S,
    const volTensorField& gradU
) const
{
	// lrans with SMALL bound
	const volScalarField lrans = max
								(
									dimensionedScalar("SMALL", dimLength, SMALL),
									sqrt(this->Rnu_/S)
								);
	// lDDES with SMALL bound
    const volScalarField lDDES = max
								(
									dimensionedScalar("SMALL", dimLength, SMALL),
									lrans - fd_*max(dimensionedScalar("ZERO", dimLength, 0), lrans - this->CDES_*this->delta())
								);

	return lrans / lDDES;
}


template<class BasicTurbulenceModel>
void WrayAgarwalDDES<BasicTurbulenceModel>::blendfactor()
{
	this->blendfactor_ = fd_;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwalDDES<BasicTurbulenceModel>::WrayAgarwalDDES
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
    WrayAgarwalDES<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cd1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd1",
            this->coeffDict_,
            8.0 // SA=8, SST=20
        )
    ),

	fd_
	(
		IOobject
		(
			"fd",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0)
	)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WrayAgarwalDDES<BasicTurbulenceModel>::read()
{
    if (WrayAgarwalDDES<BasicTurbulenceModel>::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
