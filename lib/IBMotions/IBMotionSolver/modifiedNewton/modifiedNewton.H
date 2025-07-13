/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       modifiedNewton class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
modifiedNewton class:
    - Solve newtons equations for 6DoF motion.
    - Integrate the motion of fluids inside objects.
	- Can handle objects with density close to fluids density
    - Kempe, T. J. Comput. Phys. 231, 3663–3684 (2012).
\*---------------------------------------------------------------------------*/
#ifndef modifiedNewton_H
#define modifiedNewton_H

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
							Class modifiedNewton declaration
\*---------------------------------------------------------------------------*/
class IBModel;

class modifiedNewton
:
	public IBMotionSolver
{
	eulerMesh& emesh_;

	PtrList<IBObject>& ibo_;
	
	autoPtr<IBModel>& ibmodelPtr_;

public:
	
	TypeName("modifiedNewton");
	
	modifiedNewton
	(
		word typeName,
		eulerMesh& emesh,
		PtrList<IBObject>& ibo,
		autoPtr<IBModel>& ibmodelPtr
	);

	virtual ~modifiedNewton()
	{}
	
	//- Move objects
	virtual void moveObjects();

	//- Calculate alphaSolid field
	virtual volScalarField calculateAlphaSolid();
	
	vector volIntegralU
	(
	    IBObject& ibobj,
	    const volVectorField& U
	);
	
	vector volIntegralRxU
	(
	    IBObject& ibobj,
	    const volVectorField& U
	);

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