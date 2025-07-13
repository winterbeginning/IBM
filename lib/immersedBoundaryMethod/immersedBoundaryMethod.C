/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       immersedBoundaryMethod class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(immersedBoundaryMethod, 0);
}
//----------------------------Private Member Functions-----------------------//
void Foam::immersedBoundaryMethod::createObjects(const dictionary& dict)
{
    label i=0;
    scalar nObjects = dict.size();
    objects_.setSize(nObjects);
    Info<<"Found "<<nObjects<<" IB objects!!!"<<endl;
    forAllConstIter(IDLList<entry>, dict, iter)
    {
        if (iter().isDict())
        {
            objects_.set
            (
                i++,
                IBObject::New 
                (
                    //iter().keyword(),
					iter().dict().lookup("type"),
                    emesh_,
                    iter().dict()
                )
            );
        }
    }
    objects_.setSize(i);
}

//---------------------------------Constructors------------------------------//
Foam::immersedBoundaryMethod::immersedBoundaryMethod
(
	dynamicFvMesh& mesh
)
:	
	IOdictionary
	(
		IOobject
		(
			"IBMDict",
			mesh.time().constant(),
			mesh,
        	IOobject::MUST_READ,
        	IOobject::NO_WRITE
		)
	),
	emesh_(mesh)
{

	// Create immersed objects
	createObjects(subDict("IBObjects"));

	// Create immersed boundary model
	modelPtr_  = IBModel::New
				 (
					emesh_, 
					subDict("IBModel"), 
					objects_
				 );
	
	// Create motion solver for ibobjects motions
	motionSolverPtr_ = IBMotionSolver::New
				 (
					subDict("IBModel").lookup("motionSolver"), 
					emesh_, 
					objects_, 
					modelPtr_
				 );
}

// -------------------------------Member Functions----------------------------//

Foam::volVectorField Foam::immersedBoundaryMethod::ibForce
(
	const volVectorField& U
)
{
   return modelPtr_->ibForce(U);
}

Foam::volVectorField Foam::immersedBoundaryMethod::ibForceInt()
{
   return modelPtr_->ibForceInt();
}

void Foam::immersedBoundaryMethod::multiDirectForcing
(
	volVectorField& u,
	volVectorField& ibForce
)
{
	return modelPtr_->multiDirectForcing(u, ibForce);
}

void Foam::immersedBoundaryMethod::update()
{
	return motionSolverPtr_->moveObjects();
}

Foam::volScalarField Foam::immersedBoundaryMethod::calculateAlphaSolid()
{
	return motionSolverPtr_->calculateAlphaSolid();
}


void Foam::immersedBoundaryMethod::write()
{
	return modelPtr_->write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
