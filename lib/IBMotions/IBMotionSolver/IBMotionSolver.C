/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMotionSolver class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#include "IBMotionSolver.H"
#include "wallPolyPatch.H"
#include "IBParticle.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(IBMotionSolver, 0);
    defineRunTimeSelectionTable(IBMotionSolver, dictionary);
}

//---------------------------------Constructors------------------------------//
Foam::IBMotionSolver::IBMotionSolver
(
	word typeName,
    eulerMesh& emesh,
    PtrList<IBObject>& ibo,
    autoPtr<IBModel>& ibmodelPtr
)
:
    emesh_(emesh),
    ibo_(ibo)
{}

// -------------------------------Member Functions----------------------------//

Foam::vector Foam::IBMotionSolver::objMutualRepulsive
(
    label obj1ID,   //- the object is considered
    label obj2ID    
)
{
    vector force;

    IBParticle& particle1 = refCast<IBParticle>(ibo_[obj1ID]);
    IBParticle& particle2 = refCast<IBParticle>(ibo_[obj2ID]);
	scalar R1 = particle1.R();
	scalar R2 = particle2.R();
	point& C1 = particle1.CG();
	point& C2 = particle2.CG();

    scalar d12 = mag(C1 - C2);
    scalar xi = 2.0*emesh_.h();
    scalar epsilonP = particle1.epsilonP();
    scalar epsilonPP = particle1.epsilonP();

    if (d12 <= R1+R2)
    {
        force = ((C1-C2)*(R1+R2-d12)/epsilonP);
    }
    if ( (d12 > R1+R2) && (d12 <= R1+R2+xi) )
    {
        force = ((C1-C2)*Foam::pow(R1+R2+xi-d12, 2.0)/epsilonPP);
    }
    if ( d12 > R1+R2+xi ) 
    {
        force = vector::zero;
    }

    Info<<"    Calculating repulsive force of object "<<obj1ID<<" and "<<obj2ID
                            <<" : "<<force<<endl;
    return force;
}

Foam::autoPtr<Foam::IBMotionSolver> Foam::IBMotionSolver::New
(
	word typeName,
    eulerMesh& emesh,
    PtrList<IBObject>& ibo,
    autoPtr<IBModel>& ibmodelPtr
)
{
	Info<< "Selecting motion solver: "<<typeName<< endl;
	
	dictionaryConstructorTable::iterator cstrIter = 
		dictionaryConstructorTablePtr_->find(typeName);
		
	if (cstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"IBMotionSolver::New(word typeName, dynamicFvMesh& mesh,"
			"eulerMesh& emesh,,PtrList<IBObject>& ibo,autoPtr<IBModel>& ibmodelPtr)"
		)	<< "Unknown motion solver "<< typeName <<endl <<endl
			<< "Valid solvers are : " <<endl
			<< dictionaryConstructorTablePtr_->toc()
			<< exit(FatalError);
	}
	
	return autoPtr<IBMotionSolver>(cstrIter()(typeName, emesh, ibo, ibmodelPtr));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
