/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBOStream class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBOStream_H
#define IBOStream_H

#include "dynamicFvMesh.H"
#include "IBObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
							Class IBOStream declaration
\*---------------------------------------------------------------------------*/
class IBOStream
{
	//- Private data

		dynamicFvMesh& mesh_;

		label writeInterval_;

		bool writeLagrPoints_;

		bool writeShadPoints_;

		bool writeLagrForces_;

		bool writeIBForce_;
		
		bool writeObjVTU_;

		bool writeObjData_;

	//- Private member functions

		void readDict();
				
public:

	TypeName("IBOStream");

	//- Constructor
		IBOStream
		(
			dynamicFvMesh& mesh
		);

	//- Destructor
		virtual ~IBOStream()
		{}

	//- Member functions
		void writeLagrPoints
		(
			const label objectID,
			const pointField& lPoints
		);

		void writeShadPoints
		(
			const label objectID,
			const pointField& shd1Points,
			const pointField& shd2Points
		);

		void writeLagrForces
		(
			const label objectID,
			const pointField& lPoints,
			const vectorField& lForces
		);

		void writeIBForces
		(
			const label objectID,
			const vectorField& lForces,
			const scalar dVolume
		);

		void writeObjVTU
		(
			label objectID,
		    List<point> Nodes, 
		    point pntCenter,
		    label numberFaces_,
		    List<label> nPntsOfFace_,
		    List<List<label> > PntsOfFace_
		);

		void writeObjData
		(
			const label objectID,
			const point CG,
			const vector uTransl,
			const vector uRotate
		);
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //