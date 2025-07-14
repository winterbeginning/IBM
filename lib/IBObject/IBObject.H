/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBObject class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#ifndef IBObject_H
#define IBObject_H

#include "point.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "eulerMesh.H"
#include "IBMotion.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// class IBMotion;
/*---------------------------------------------------------------------------*\
                        Class  IBObjectSelection Declaration
\*---------------------------------------------------------------------------*/

class IBObject
{

public:
	//- Runtime type information
	TypeName("IBObject");

	//- Declare runtime constructor selection table
		declareRunTimeSelectionTable
		(
			autoPtr,
			IBObject,
			dictionary,
			(
				const word& typeName,
				eulerMesh& emesh,
				const dictionary& dict
			),
			(typeName, emesh, dict)
		);

	//- Constructor
		IBObject
		(
			const word& typeName,
			eulerMesh& emesh,
			const dictionary& dict
		);

	//- Selectors
		static autoPtr<IBObject> New 
		(
			const word& typeName,
			eulerMesh& emesh,
			const dictionary& dict
		);

	//- Destructor
		virtual ~IBObject()
		{}

	//- Member functions

			virtual void findNeiCells() = 0;

			virtual void findSolidCells() = 0;

			virtual void findSolidCellsExt() = 0;

		//- Access
			//- Object's geometry and identity
				virtual const word& name() = 0;
				virtual const word& objectType() = 0; 
				virtual const scalar& rho() = 0;

				//- Return object's volume
				virtual const scalar& V() = 0;
				
				//- Return Moment of Inertia
				virtual const vector& Ip() = 0;
				
				//- Return the center of gravity
				virtual point& CG() = 0;
				
				//- Return number of Lagrangian points
				virtual label& nPoints() = 0;

				//- Return Lagrang points
				virtual pointField& lPoints() = 0;

				//- Return the discrete volume of Lagrange points
				virtual const scalar& dVL() = 0;

				//- Return neighbour cells
				virtual const labelListList& neiCells() = 0;

				//- Return solid cells
				virtual const labelList& solidCells() = 0;

				//- Return solid cells extend
				virtual const labelList& solidCellsExt() = 0;

			//- Object's motion
				virtual bool movable() = 0;
				//- Return type of motions
				virtual PtrList<IBMotion>& motions() = 0;

				//- Return translational velocity
				virtual const vector& uTranslate() = 0;

				//- Return rotational velocity
				virtual const vector& uRotate() = 0;

				//- Return velocity at Lagrangian points
				virtual const vectorField& uBoundary() = 0;
		
				//- Return true if wall collision force is considered
				virtual bool calcWallRepulsive() = 0;

				//- Return repulsive force with wall boundary	
				virtual const vector wallRepulsiveForce() = 0;	

				//- Return stiffness coefficient of wall collision
				virtual const scalar& epsilonW() = 0;

				//- Return stiffness coefficient of object mutual collision
				virtual const scalar& epsilonP() = 0;

			//- Data for outputing
				//- Return number of faces to write object
				virtual const label& nFaces() = 0;

				//- Return number of points of each object's face
				virtual const labelList& nPointsOfFaces() = 0;

				//- Return coordinate of points of each object's face
				virtual const labelListList& pointOfFace() = 0;
		
		//- Managing private data
			//- Update motion at each timeStep, Uhlmann (2005)
			virtual void updateObjectMotionUhlmann
			(
				const vectorField& ForceLagrang,
				vector& repulsiveForce,
				const scalar& rhoFluid,
				const vector& g
			) = 0;
			
			//- Update motion at each timeStep, Tobias (2012)
			virtual void updateObjectMotionTobias
			(
				const vectorField& ForceLagrang,
				vector& repulsiveForce,
				const scalar& rhoFluid,
				const vector& g,
				const vector& volIntegralU,
				const vector& volIntegralRxU

			) = 0;

			//- Update object's velocity
			virtual void updateVelocity() = 0;

			//- Move points
			virtual void movePoints() = 0;

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
