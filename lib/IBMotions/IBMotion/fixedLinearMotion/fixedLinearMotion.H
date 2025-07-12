/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - fixedLinearMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef fixedLinearMotion_H
#define fixedLinearMotion_H

#include "IBMotion.H"
#include "dictionary.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
						Class fixedLinearMotion declaration
\*---------------------------------------------------------------------------*/
class fixedLinearMotion
:
	public IBMotion
{
	//- Private data

		//- Translational velocity
		vector uTrans_;

		//- Rotational velocity, uRot = 0
		vector uRot_;

public:

	//- Runtime type information
	TypeName("fixedLinearMotion");

	//- Constructor
		fixedLinearMotion
		(
			const word& typeName,
			IBObject& obj,
			const dictionary& dict
		);

	//- Destructor
		virtual ~fixedLinearMotion()
		{}

	//- Member functions

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