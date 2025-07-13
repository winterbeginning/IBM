/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       immersedBoundaryMethod class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryMethod_H
#define immersedBoundaryMethod_H

#include "dynamicFvMesh.H"
#include "autoPtr.H"
#include "eulerMesh.H"
#include "IBObject.H"
#include "IBModel.H"
#include "IBMotionSolver.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class immersedBoundaryMethod declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryMethod
:
	IOdictionary
{
	// Private data

		//- Euler mesh including cell size, cell volume, rho fluid, ...
		eulerMesh emesh_;

		//- IB Objects
		PtrList<IBObject> objects_;

		//- Pointer to selected IB model 
		autoPtr<IBModel> modelPtr_;

		//- Pointer to selected motion solver		
		autoPtr<IBMotionSolver> motionSolverPtr_;

	// Private member functions

		// Create ib objects
		void createObjects(const dictionary& dict);
	
public:
	
	TypeName("immersedBoundaryMethod");
	
	immersedBoundaryMethod(dynamicFvMesh& mesh);

	virtual ~immersedBoundaryMethod()
	{}

	virtual volVectorField ibForce(const volVectorField& U);

	virtual volVectorField ibForceInt();

	virtual void multiDirectForcing
	(
		volVectorField& u,
		volVectorField& ibForce
	);

	virtual void update();

	virtual volScalarField calculateAlphaSolid();

	virtual void write();

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
