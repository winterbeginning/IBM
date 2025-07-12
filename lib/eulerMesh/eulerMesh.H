/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       eulerMesh class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#ifndef eulerMesh_H
#define eulerMesh_H

#include "dynamicFvMesh.H"
#include "volMesh.H"
#include "GeometricField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class eulerMesh declaration
\*---------------------------------------------------------------------------*/

class eulerMesh
{
	
	//- Private data
		dynamicFvMesh& mesh_;
		const IOdictionary transportProperties_;
		scalar rhoFluid_;
		scalar rhoAir_;		
		vector g_;
		scalar h_;
		scalar dV_;

	//- Private funcs
		//- Read variables from transportProperties directory
		void readEnvironmentVariables();

		//- Calculate the lagrange cell size
		scalar cellSize(label cellID);

public:
	
	TypeName("eulerMesh");
	
	//-Constructor
		eulerMesh(dynamicFvMesh& mesh);

	//- Destructor
		virtual ~eulerMesh()
		{}

	//- Member funcs

		//- Manipulate
			void updateEulerMeshInfo();

		//- Access
			inline dynamicFvMesh& mesh() {return mesh_;}
			inline const scalar& rhoFluid() {return rhoFluid_;}
			inline const vector& g() {return g_;}
			inline scalar& h() {return h_;}			
			inline scalar& dV(){return dV_;}

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
