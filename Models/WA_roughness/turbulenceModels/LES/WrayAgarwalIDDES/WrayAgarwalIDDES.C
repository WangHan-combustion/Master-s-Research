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

#include "WrayAgarwalIDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalIDDES<BasicTurbulenceModel>::alpha() const
{
    return max
    (
        0.25 - this->y_/static_cast<const volScalarField&>(IDDESDelta_.hmax()),
        scalar(-5)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalIDDES<BasicTurbulenceModel>::ft
(
    const volScalarField& magGradU
) const
{
    return tanh(pow3(sqr(ct_)*rd(this->nut_, magGradU)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalIDDES<BasicTurbulenceModel>::fl
(
    const volScalarField& magGradU
) const
{
    return tanh(pow(sqr(cl_)*rd(this->nu(), magGradU), 10));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalIDDES<BasicTurbulenceModel>::rd
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
tmp<volScalarField> WrayAgarwalIDDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
	if(IDES_)
	{
		return max(scalar(1), 1 - tanh(pow3(Cd1_*rd(this->nuEff(), magGradU))));
	}
	else
	{
    	return 1 - tanh(pow3(Cd1_*rd(this->nuEff(), magGradU)));
	}
}

template<class BasicTurbulenceModel>
void WrayAgarwalIDDES<BasicTurbulenceModel>::precalculations
(
	const volScalarField& S,
    const volTensorField& gradU
)
{
	// alpha
    const volScalarField alpha(this->alpha());
	// e^(alpha^2)
    const volScalarField expTerm(exp(sqr(alpha)));

    const volScalarField magGradU(mag(gradU));

	// fe1
    tmp<volScalarField> fHill =
        2*(pos(alpha)*pow(expTerm, -11.09) + neg(alpha)*pow(expTerm, -9.0));
	
	// fb
    tmp<volScalarField> fStep = min(2*pow(expTerm, -9.0), scalar(1));
	
	// Store fd_
	fd_ = fd(magGradU);

	// fd_tilda
    const volScalarField fHyb(max(1 - fd_, fStep));

	// fe2
    tmp<volScalarField> fAmp = 1 - max(ft(magGradU), fl(magGradU));
  
	// fe
	tmp<volScalarField> fRestore = max(fHill - 1, scalar(0))*fAmp;

	// Store fdtilda_ and fe_
	fdtilda_ = fHyb;
	fe_ = fRestore;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalIDDES<BasicTurbulenceModel>::fdes
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
	// liddes with SMALL bound
    const volScalarField liddes = max
								(
									dimensionedScalar("SMALL", dimLength, SMALL),
									fdtilda_*(1 + fe_)*lrans
								  + (1 - fdtilda_)*this->CDES_*this->delta()
								);

	return lrans / liddes;
}


template<class BasicTurbulenceModel>
void WrayAgarwalIDDES<BasicTurbulenceModel>::blendfactor()
{
	this->blendfactor_ = 1 - fdtilda_;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwalIDDES<BasicTurbulenceModel>::WrayAgarwalIDDES
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

    fwStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fwStar",
            this->coeffDict_,
            0.424
        )
    ),

    cl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cl",
            this->coeffDict_,
            3.55
        )
    ),

    ct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ct",
            this->coeffDict_,
            1.63
        )
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

	IDES_
    (
        Switch::lookupOrAddToDict
        (
            "IDES",
            this->coeffDict_,
            false
        )
    ),

	fdtilda_
	(
		IOobject
		(
			"fdtilda",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0)
	),

	fe_
	(
		IOobject
		(
			"fe",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		this->mesh_,
		dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0)
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
	),

    IDDESDelta_(refCast<IDDESDelta>(this->delta_()))
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }

	if(IDES_)
	{
		Info<< "IDES mode on, running IDES" << endl << endl;
	}
	else
	{
		Info<< "IDES mode off, running IDDES" << endl << endl;
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WrayAgarwalIDDES<BasicTurbulenceModel>::read()
{
    if (WrayAgarwalIDDES<BasicTurbulenceModel>::read())
    {
        //fwStar_.readIfPresent(this->coeffDict());
        //cl_.readIfPresent(this->coeffDict());
        //ct_.readIfPresent(this->coeffDict());
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
