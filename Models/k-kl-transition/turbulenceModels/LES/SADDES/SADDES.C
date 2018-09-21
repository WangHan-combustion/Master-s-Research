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

#include "SADDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SADDES<BasicTurbulenceModel>::rd
(
    const volScalarField& magGradU
) const
{
    tmp<volScalarField> tr
    (
        min
        (
            this->nuEff()
           /(
               max
               (
                   magGradU,
                   dimensionedScalar("SMALL", magGradU.dimensions(), SMALL)
               )
              *sqr(this->kappa_*this->y_)
            ),
            scalar(10)
        )
    );
    tr().boundaryField() == 0.0;

    return tr;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SADDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    return 1 - tanh(pow3(8*rd(magGradU)));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SADDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
)
{
	fd_ = fd(mag(gradU));

    return max
    (
        this->y_
      - fd_
       *max
        (
            this->y_ - this->CDES_*this->delta(),
            dimensionedScalar("zero", dimLength, 0)
        ),
        dimensionedScalar("small", dimLength, SMALL)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SADDES<BasicTurbulenceModel>::SADDES
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
    SADES<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
	fd_
	(
		IOobject
		(
			"fd",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
	)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
