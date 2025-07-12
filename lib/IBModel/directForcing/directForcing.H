/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       directForcing class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#ifndef directForcing_H
#define directForcing_H

#include "IBModel.H"
#include "IBObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class IBModel declaration
\*---------------------------------------------------------------------------*/

class directForcing
: 
	public IBModel
{
	//- Private data
		
		IBOStream ibOStream_;

		//- Reference to euler mesh
		eulerMesh& emesh_;

		dynamicFvMesh& mesh_;

		//- Reference to IB objects
		PtrList<IBObject>& objects_;

		//- Number of directForcing loop
		label nMDF_;
		
		List<vectorField> uLagrang_;

		List<vectorField> fLagrang_;
		
	//- Disallow default bitwise copy construct
	directForcing(const directForcing&);
	
	//- Disallow default bitwis assignment
	void operator=(const directForcing&);

	
public:

	//- Runtime type information
	TypeName("directForcing");

	//- Constructor
		directForcing
		(
			eulerMesh& emesh,
			const dictionary& dict,
			PtrList<IBObject>& ibo
		);

	//- Destructor
		virtual ~directForcing()
		{}

	//- Member functions

		virtual const List<vectorField>& ULagrange(){return uLagrang_;}
		virtual const List<vectorField>& FLagrange(){return fLagrang_;}
		virtual volVectorField ibForce(const volVectorField& U);

		virtual void write();
		virtual volVectorField ibForceInt();
		
		// //- Perform multi-Firect forcing
		virtual void multiDirectForcing
		(
			volVectorField& u,
			volVectorField& ibForce
		);

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
