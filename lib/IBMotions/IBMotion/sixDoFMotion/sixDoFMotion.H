/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - sixDoFMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef sixDoFMotion_H
#define sixDoFMotion_H

#include "IBMotion.H"
#include "IBObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
						Class sixDoFMotion declaration
\*---------------------------------------------------------------------------*/
class sixDoFMotion
:
	public IBMotion
{
	//- Private data

		//- Reference to Immersed object
		IBObject& obj_;

		//- Translational velocity
		vector uTrans_;

		//- Rotational velocity
		vector uRot_;

		//- Constraint
		vector constraintDir_;

		void readConstraint(const dictionary& dict);

public:

	//- Runtime type information
	TypeName("sixDoFMotion");

	//- Constructor
		sixDoFMotion
		(
			const word& typeName,
			IBObject& obj,
			const dictionary& dict
		);

	//- Destructor
		virtual ~sixDoFMotion()
		{}

	//- Member functions

		//- Solve Newton equation of motion, M. Uhlmann 2005
		void updateMotion
		(
			dynamicFvMesh& mesh,
			const vectorField& Fk, 
			vector repulsiveForce,
			scalar rhoFluid,
			vector g
		);

		//- Solve Newton equation of motion, Tobias Kempe & J.Frohlich 2012
		void updateMotion
		(
			dynamicFvMesh& mesh,
			const vectorField& Fk, 
			vector repulsiveForce,
			scalar rhoF,
			vector g,
			vector volIntegralU,
			vector volIntegralRxU
		);

		//- Return translational velocity
		virtual vector& translationalVelocity();

		//- Return rotational velocity
		virtual vector& rotationalVelocity();

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //