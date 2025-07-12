/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBMotion_H
#define IBMotion_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
// #include "IBObject.H"
#include "vector.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	class IBObject;
/*---------------------------------------------------------------------------*\
                        	Class  IBMotion Declaration
\*---------------------------------------------------------------------------*/

class IBMotion
{

	word motionType_;

public:
	//- Runtime type information
	TypeName("IBMotion");

	//- Declare runtime constructor selection table
		declareRunTimeSelectionTable
		(
			autoPtr,
			IBMotion,
			dictionary,
			(
				const word& name,
				IBObject& obj,
				const dictionary& dict
			),
			(name, obj, dict)
		);

	//- Constructor
		IBMotion
		(
			const word& name,
			IBObject& obj,
			const dictionary& dict
		);

	//- Selectors
		static autoPtr<IBMotion> New 
		(
			const word& name,
			IBObject& obj,
			const dictionary& dict
		);

	//- Destructor
		virtual ~IBMotion()
		{}

	//- Member functions

		//- Return object's motion type
		const word& motionType(){return motionType_;}

		//- Return translational velocity
		virtual vector& translationalVelocity() = 0; 
		
		//- Return rotational velocity
		virtual vector& rotationalVelocity() = 0;

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
