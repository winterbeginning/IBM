/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBParticles class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#ifndef IBParticle_H
#define IBParticle_H

#include "IBObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
							Class IBParticle declaration
\*---------------------------------------------------------------------------*/

class IBParticle
:
	public IBObject
{

		const scalar PI;

		eulerMesh& emesh_;

	//- Particle parameters
		word name_;
		
		word objectType_;

		//- Center
		point center_;

		//- Radius
		scalar R_;

		//- Density 
		scalar rho_;

		//- Particle volume
		scalar V_;

		//- Moment of Inertia
		vector Ip_;

		//- Number of Lagrang points
		label nPoints_;

		bool autoPoint_;

		//- Lagrang points
		pointField lPoints_;

		//- Discrete volume associate with L points
		scalar dVL_;
		
		//- Neighbouring cells
		labelListList neiCells_;

		//- Solid cells
		labelList solidCells_;

		//- Solid cells extended
		labelList solidCellsExt_;

	//- Motions
		//- Do object move?
		bool movable_;
		
		//- Pointer to object's motion
		PtrList<IBMotion> motions_;

		//- Translational velocity
		vector uTranslate_;

		//- Rotational velocity
		vector uRotate_;
		
		//- Velocity at Lagrangian points
		vectorField uBoundary_;
		
	//- Points connectivity for writing geometry
		label nFaces_;
		
		labelList nPointsOfFaces_;
		
		labelListList pointOfFace_;

	//- Collision force
		bool calcWallRepulsive_;

		scalar epsilonW_;

		scalar epsilonP_;

	//- Private functions
		
		void initialize(const dictionary& dict);

		//- Create 2D particles automatically
		void createParticle2D();
		
		//- Create 3D particles from input TXT file
		void createParticle3D();

		void calcNLagrangPoints();

		//- Distribute evenly Lagrangian points on sphere surface
		pointField createUnitSphereEqAreaPartition(label nPoints);
		
		scalar capArea(scalar angle);

		//- Create motion
		void addMotions(const dictionary& dict);

public:

	//- Runtime type information
	TypeName("IBParticle");

	//- Constructors
		IBParticle
		(	
			const word& typeName,
			// dynamicFvMesh& mesh,
			eulerMesh& emesh,
			const dictionary& dict
		);
	
	//- Destructor
		virtual ~IBParticle()
		{};

	//- Member functions
		
		//- Manipulate data

			virtual void findNeiCells();

			virtual void findSolidCells();

			virtual void findSolidCellsExt();

		//- Access data

			virtual const word& name();
			virtual const word& objectType();			
			virtual const scalar& rho();
			virtual const scalar& V();
			virtual const vector& Ip();
			virtual const scalar& R();
			virtual point& CG();
			virtual label& nPoints();
			virtual pointField& lPoints();
			virtual const scalar& dVL();

			virtual const labelListList& neiCells();
			virtual const labelList& solidCells();
			virtual const labelList& solidCellsExt();

			virtual bool movable();
			virtual PtrList<IBMotion>& motions();
			virtual const vector& uTranslate();
			virtual const vector& uRotate();
			virtual const vectorField& uBoundary() ;
			virtual bool calcWallRepulsive();
			virtual const vector wallRepulsiveForce();
			virtual const scalar& epsilonW();
			virtual const scalar& epsilonP();

			virtual const label& nFaces();
			virtual const labelList& nPointsOfFaces();
			virtual const labelListList& pointOfFace();

			virtual void updateObjectMotionUhlmann
			(
				const vectorField& ForceLagrang,
				vector& repulsiveForce,
				const scalar& rhoFluid,
				const vector& g
			);

			virtual void updateObjectMotionTobias
			(
				const vectorField& ForceLagrang,
				vector& repulsiveForce,
				const scalar& rhoFluid,
				const vector& g,
				const vector& volIntegralU,
				const vector& volIntegralRxU

			);

			virtual void updateVelocity();
			virtual void movePoints();

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //