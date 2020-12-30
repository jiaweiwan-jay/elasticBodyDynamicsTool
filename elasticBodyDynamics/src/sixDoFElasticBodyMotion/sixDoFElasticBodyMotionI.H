/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::tensor Foam::sixDoFElasticBodyMotion::rotationTensorX
(
    scalar phi
) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


inline Foam::tensor Foam::sixDoFElasticBodyMotion::rotationTensorY
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


inline Foam::tensor Foam::sixDoFElasticBodyMotion::rotationTensorZ
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}

inline Foam::tensor Foam::sixDoFElasticBodyMotion::rotationTensor
(
    vector omega
) const
{
    tensor R = rotationTensorX(0.5*omega.x());
    tensor Q = tensor::I;
    
    Q = Q & R;

    R = rotationTensorY(0.5*omega.y());
    Q = Q & R;

    R = rotationTensorZ(omega.z());
    Q = Q & R;

    R = rotationTensorY(0.5*omega.y());
    Q = Q & R;

    R = rotationTensorX(0.5*omega.x());
    Q = Q & R;
    
    return Q;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::generalizedMotion& Foam::sixDoFElasticBodyMotion::gMotion()
{
    return gMotion_;
}

inline const Foam::point& Foam::sixDoFElasticBodyMotion::origin() const
{
    return origin_;
}

inline const Foam::vector& Foam::sixDoFElasticBodyMotion::direction() const
{
    return direction_;
}

inline const Foam::label& Foam::sixDoFElasticBodyMotion::nNode() const
{
    return nNode_;
}

inline const Foam::label& Foam::sixDoFElasticBodyMotion::nMode() const
{
    return nMode_;
}

inline const Foam::List<Foam::scalar>& Foam::sixDoFElasticBodyMotion::frequency() const
{
    return frequency_;
}

inline const Foam::List<Foam::scalar>& Foam::sixDoFElasticBodyMotion::dampingRatio() const
{
    return dampingRatio_;
}

inline const Foam::List<Foam::vector6D>& Foam::sixDoFElasticBodyMotion::mass() const
{
    return mass_;
}

inline const Foam::List<Foam::scalar>& Foam::sixDoFElasticBodyMotion::length() const
{
    return length_;
}

inline const Foam::List<Foam::List<Foam::vector6D> >& Foam::sixDoFElasticBodyMotion::mode() const
{
    return mode_;
}

inline const Foam::generalizedMotion& Foam::sixDoFElasticBodyMotion::gMotion() const
{
    return gMotion_;
}

inline const Foam::List<Foam::vector6D>& Foam::sixDoFElasticBodyMotion::u() const
{
    return u_;
}

inline const Foam::List<Foam::vector6D>& Foam::sixDoFElasticBodyMotion::v() const
{
    return v_;
}

inline const Foam::List<Foam::vector6D>& Foam::sixDoFElasticBodyMotion::a() const
{
    return a_;
}

inline Foam::point Foam::sixDoFElasticBodyMotion::transform
(
    const point& initialPoint
) const
{
    return
    (
        initialPoint
    );
}
// ************************************************************************* //
