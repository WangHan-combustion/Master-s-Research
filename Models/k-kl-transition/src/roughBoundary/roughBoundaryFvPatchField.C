/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "roughBoundaryFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::roughBoundaryFvPatchField<Type>::roughBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    ks_(1e-6)
{}


template<class Type>
Foam::roughBoundaryFvPatchField<Type>::roughBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    ks_(readScalar(dict.lookup("ks")))
{
    evaluate();
}


template<class Type>
Foam::roughBoundaryFvPatchField<Type>::roughBoundaryFvPatchField
(
    const roughBoundaryFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    ks_(ptf.ks_)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::roughBoundaryFvPatchField<Type>::roughBoundaryFvPatchField
(
    const roughBoundaryFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    ks_(ptf.ks_)
{}


template<class Type>
Foam::roughBoundaryFvPatchField<Type>::roughBoundaryFvPatchField
(
    const roughBoundaryFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    ks_(ptf.ks_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::roughBoundaryFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type> media(Field<Type>(this->patch().size(),pTraits<Type>::one) - 
            Field<Type>(this->patch().size(),pTraits<Type>::one)/
                        (0.03 * ks_ * this->patch().deltaCoeffs())
                        );

    label maxIter = 100;

    label it = 0;

    for (it = 0; it < maxIter; it++)
    {
        Field<Type> grad = (*this) / (0.03 * ks_); 

        Type sum0 = pTraits<Type>::zero;

        forAll(this->patch(),fI)
        {
            sum0 += (*this)[fI];
        }

        Field<Type> patchInternal(this->patchInternalField());

        Field<Type>::operator=
            (
             patchInternal + grad / this->patch().deltaCoeffs()
            ); 

        Type sum = pTraits<Type>::zero;

        forAll(this->patch(),fI)
        {
            sum += (*this)[fI];
        }

        scalar residual = mag(sum0) - mag(sum);

        if (mag(residual) / (mag(sum) + 1e-6) < 0.001)
        {
            break;
        }

    }

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::roughBoundaryFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    //return tmp<Field<Type>>(new Field<Type>(this->size(), pTraits<Type>::one));

    return *this;

}

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::roughBoundaryFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::roughBoundaryFvPatchField<Type>::gradientInternalCoeffs() const
{
//  return -pTraits<Type>::one*this->patch().deltaCoeffs();
//  return tmp<Field<Type>>
//  (
//      new Field<Type>(this->size(), Zero)
//  );
    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::roughBoundaryFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return *this;
//  return this->patch().deltaCoeffs()*(*this);
}


template<class Type>
void Foam::roughBoundaryFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("ks", os);
}


// ************************************************************************* //
