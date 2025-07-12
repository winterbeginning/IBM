/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBObject class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/
#include "error.H"
// #include "IBMotion.H"
#include "IBObject.H"
#include "GeometricField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBObject, 0);
    defineRunTimeSelectionTable(IBObject, dictionary);
}
//---------------------------------Constructors------------------------------//
Foam::IBObject::IBObject
(
	const word& typeName,
	eulerMesh& emesh,
	const dictionary& dict
)
{}

//------------------------------------Selectors------------------------------//

Foam::autoPtr<Foam::IBObject> Foam::IBObject::New
(
	const word& typeName,
	eulerMesh& emesh,
	const dictionary& dict
)
{
	dictionaryConstructorTable::iterator dictCstrIter = 
		dictionaryConstructorTablePtr_->find(typeName);

	if(dictCstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"IBObject::New(const eulerMesh& emesh, "
			"const word& name, const dictionary& dict)"
		)	<< "Unknown object type "<< typeName 
			<< endl
			<< "Valid object type are: "<<endl
			<<dictionaryConstructorTablePtr_->toc()
			<<exit(FatalError);
	}

	return autoPtr<IBObject>(dictCstrIter()(typeName, emesh, dict));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
