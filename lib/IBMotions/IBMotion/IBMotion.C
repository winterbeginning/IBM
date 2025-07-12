/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "IBMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBMotion, 0);
    defineRunTimeSelectionTable(IBMotion, dictionary);
}

//---------------------------------Constructors------------------------------//
Foam::IBMotion::IBMotion
(
	const word& typeName,
	IBObject& obj,
	const dictionary& dict
)
:
	motionType_(typeName)
{}

//------------------------------------Selectors------------------------------//

Foam::autoPtr<Foam::IBMotion> Foam::IBMotion::New
(
	const word& typeName,
	IBObject& obj,
	const dictionary& dict
)
{
	dictionaryConstructorTable::iterator dictCstrIter = 
		dictionaryConstructorTablePtr_->find(typeName);

	if(dictCstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"IBMotion::New(const fvMesh& mesh, "
			"const word& name, const dictionary& dict)"
		)	<< "Unknown object type "<< typeName 
			<< endl
			<< "Valid object type are: "<<endl
			<<dictionaryConstructorTablePtr_->toc()
			<<exit(FatalError);
	}

	return autoPtr<IBMotion>(dictCstrIter()(typeName, obj, dict));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
