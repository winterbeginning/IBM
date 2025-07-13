/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMotionSolver class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#ifndef IBMotionSolver_H
#define IBMotionSolver_H

// #include "IOmanip.H"
#include "autoPtr.H"
#include "dynamicFvMesh.H"
#include "eulerMesh.H"
#include "IBObject.H"
#include "IBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class IBMotionSolver declaration
\*---------------------------------------------------------------------------*/
// class IBModel;

class IBMotionSolver
{
	eulerMesh& emesh_;
	
	PtrList<IBObject>& ibo_;
	
public:
	//- Runtime type information
	TypeName("IBMotionSolver");

		declareRunTimeSelectionTable
		(
			autoPtr,
			IBMotionSolver,
			dictionary,
			(
				word typeName,
				eulerMesh& emesh,
				PtrList<IBObject>& ibo,
				autoPtr<IBModel>& ibmodelPtr
			),
			(typeName, emesh, ibo, ibmodelPtr)
		);

	//- Constructor
		IBMotionSolver
		(
			word typeName,
			eulerMesh& emesh,
			PtrList<IBObject>& ibo,
			autoPtr<IBModel>& ibmodelPtr
		);

	//- Selector 
		static autoPtr<IBMotionSolver> New
		(
			word typeName,
			eulerMesh& emesh,
			PtrList<IBObject>& ibo,
			autoPtr<IBModel>& ibmodelPtr
		);

	//- Destructor
		virtual ~IBMotionSolver()
		{}
	
	//- Member functions
		//- Calculate particle mutual interactions, Wan, Turek 2007
		vector objMutualRepulsive(label obj1ID, label obj2ID);

		//- Move objects
		virtual void moveObjects() = 0;

		//- Calculate alphaSolid field
		virtual volScalarField calculateAlphaSolid() = 0;

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //