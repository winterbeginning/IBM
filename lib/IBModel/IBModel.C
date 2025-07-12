/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBModel class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#include "IBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(IBModel, 0);
	defineRunTimeSelectionTable(IBModel, dictionary);
}

//---------------------------------Constructors------------------------------//

Foam::IBModel::IBModel
(
    eulerMesh& emesh,
	const dictionary& dict,
    PtrList<IBObject>& ibo
)
:
    emesh_(emesh)
{}

// -------------------------------Member Functions----------------------------//
Foam::scalar Foam::IBModel::deltaFunc(point p_Eul, point p_Lag)
{
    scalar deltaX = deltaFunc1D(p_Eul.x(), p_Lag.x());
    scalar deltaY = deltaFunc1D(p_Eul.y(), p_Lag.y());
    scalar deltaZ = deltaFunc1D(p_Eul.z(), p_Lag.z());

    if (emesh_.mesh().nGeometricD() == 3)
        return (deltaX*deltaY*deltaZ);
    else
        return (deltaX*deltaY);
}

Foam::scalar Foam::IBModel::deltaFunc1D(scalar x_Eul, scalar x_Lag)
{
    const scalar r = mag((x_Eul - x_Lag) / emesh_.h());
    scalar phiR = 0.0;

    if( r < 2.0 )
    {
        phiR = (1.0/4.0)*(1.0 + Foam::cos(constant::mathematical::pi*0.5*r));
    }
    // if (r < 0.5)
    // {
    //     phiR = 1.0/3.0 * (1.0 + Foam::sqrt(-3.0 * (r*r) + 1.0));
    // }
    // else if (r <= 1.5)
    // {
    //     phiR = 1.0/6.0 * (5.0 - 3.0*r - 
    //         Foam::sqrt(-3.0* (1.0-r)*(1.0-r) + 1.0));    
    // }
    else
    {
        phiR = 0.0;
    }

    scalar delta1D = phiR / emesh_.h();

    return delta1D;
}

Foam::autoPtr<Foam::IBModel> Foam::IBModel::New
(
    eulerMesh& emesh,
	const dictionary& dict,
    PtrList<IBObject>& ibo
)
{
	word typeName = dict.lookup("model");
	
	Info<< "Selecting Immersed Boundary approaching method: "<<typeName<< endl;
	
	dictionaryConstructorTable::iterator cstrIter = 
		dictionaryConstructorTablePtr_->find(typeName);
		
	if (cstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"IBModel::New(const word& name, const dynamicFvMesh& mesh,"
			"const dictionary& dict,PtrList<IBObject>& ibo,)"
		)	<< "Unknown approaching method "<< typeName <<endl <<endl
			<< "Valid approaching method are : " <<endl
			<< dictionaryConstructorTablePtr_->toc()
			<< exit(FatalError);
	}
	
	return autoPtr<IBModel>(cstrIter()(emesh, dict, ibo));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
