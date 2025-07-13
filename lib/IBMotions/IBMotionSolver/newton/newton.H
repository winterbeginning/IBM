/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       newton class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
newton class:
    - Solve newtons equations for 6DoF motion.
    - Cannot handle objects that have density closed to fluids density
    - M.Uhlmann, J. Comput. Phys. 209, 448–476 (2005).
\*---------------------------------------------------------------------------*/

#ifndef newton_H
#define newton_H

// #include "IOmanip.H"
#include "autoPtr.H"
#include "dynamicFvMesh.H"
#include "eulerMesh.H"
#include "IBObject.H"
#include "IBModel.H"
#include "IBMotionSolver.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class newton declaration
\*---------------------------------------------------------------------------*/
class IBModel;

class newton
:
	public IBMotionSolver
{
	eulerMesh& emesh_;

	PtrList<IBObject>& ibo_;
	
	autoPtr<IBModel>& ibmodelPtr_;

public:
	
	TypeName("newton");

	//- Constructor
		newton
		(
			word typeName,
			eulerMesh& emesh,
			PtrList<IBObject>& ibo,
			autoPtr<IBModel>& ibmodelPtr
		);
	//- Destructor
		virtual ~newton()
		{}
	
	//- Member functions
		//- Move objects
		virtual void moveObjects();

		virtual volScalarField calculateAlphaSolid();

		scalar volFraction
		(
			IBObject& ibobj, 
			label cellID
		);
		
		scalar HeavisideFunc
		(
			scalar phi
		);
		
		scalar LevelSetFunc
		(
			scalar r, 
			point p, 
			point c
		);



}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //