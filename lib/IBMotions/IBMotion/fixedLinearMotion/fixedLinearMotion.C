/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - fixedLinearMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "fixedLinearMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fixedLinearMotion, 0);
    addToRunTimeSelectionTable(IBMotion, fixedLinearMotion, dictionary);
}
//---------------------------------Constructors------------------------------//
Foam::fixedLinearMotion::fixedLinearMotion
(
	const word& typeName,
	IBObject& obj,
	const dictionary& dict
)
:
	IBMotion(typeName, obj, dict),
	uTrans_(dict.lookup("U")),
	uRot_(vector::zero)
{}

//-------------------------------Member functions----------------------------//
Foam::vector& Foam::fixedLinearMotion::translationalVelocity()
{
	return uTrans_;
}

Foam::vector& Foam::fixedLinearMotion::rotationalVelocity()
{
	return uRot_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
