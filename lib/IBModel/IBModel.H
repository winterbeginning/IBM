/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBModel class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#ifndef IBModel_H
#define IBModel_H

#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "fvPatchField.H"
#include "eulerMesh.H"
#include "IBOStream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class IBModel declaration
\*---------------------------------------------------------------------------*/

class IBModel
{
	//- Private data
	eulerMesh& emesh_;

	//- Disallow default bitwise copy construct
	IBModel(const IBModel&);
	
	//- Disallow default bitwise assignment
	void operator=(const IBModel&);
		
public:		
	//- Runtime type information
	TypeName("IBModel");

		declareRunTimeSelectionTable
		(
			autoPtr,
			IBModel,
			dictionary,
			(
				eulerMesh& emesh,
				const dictionary& dict,
				PtrList<IBObject>& ibo
			),
			(emesh, dict, ibo)
		);
		
	//- Constructor
		IBModel
		(
			eulerMesh& emesh,
			const dictionary& dict,
			PtrList<IBObject>& ibo
		);
		
	//- Selector
		static autoPtr<IBModel> New
		(
			eulerMesh& emesh,
			const dictionary& dict,
			PtrList<IBObject>& ibo
		);
		
	//- Destructor
		virtual ~IBModel()
		{}

	//- Member functions

		scalar deltaFunc(point pE, point pL);

		scalar deltaFunc1D(scalar xE, scalar xL);

		virtual void write() = 0;

		virtual volVectorField ibForce(const volVectorField& U) = 0;

		virtual volVectorField ibForceInt() = 0;

		virtual void multiDirectForcing
		(
			volVectorField& u,
			volVectorField& ibForce
		) = 0;

		virtual const List<vectorField>& ULagrange() = 0;

		virtual const List<vectorField>& FLagrange() = 0;

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
