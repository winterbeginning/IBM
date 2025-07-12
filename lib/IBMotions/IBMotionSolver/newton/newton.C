/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       newton class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#include "newton.H"
#include "wallPolyPatch.H"
#include "IBParticle.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(newton, 0);
    addToRunTimeSelectionTable(IBMotionSolver, newton, dictionary);
}

//---------------------------------Constructors------------------------------//
Foam::newton::newton
(
    word typeName,
    eulerMesh& emesh,
    PtrList<IBObject>& ibo,
    autoPtr<IBModel>& ibmodelPtr
)
:
    IBMotionSolver(typeName, emesh, ibo, ibmodelPtr),
    emesh_(emesh),
    ibo_(ibo),
    ibmodelPtr_(ibmodelPtr)
{
}

// -------------------------------Member Functions----------------------------//

void  Foam::newton::moveObjects()
{
    for(int i=0; i<ibo_.size(); i++)
    {
        if (ibo_[i].movable())
        {
            Info<< "Moving object "<<ibo_[i].name()<<endl;
            
            vector repulsiveForce(vector::zero);
            
            for(int j=0; j<ibo_.size(); j++)
            {
                if (i == j)
                    continue;
                else
                {
                    repulsiveForce += objMutualRepulsive(i, j);
                }
            }
            
            ibo_[i].updateObjectMotionUhlmann  //只有设置6自由度才会调用
            (
                ibmodelPtr_->FLagrange()[i],
                repulsiveForce,
                emesh_.rhoFluid(),
                emesh_.g()
            );

            ibo_[i].updateVelocity();
            ibo_[i].movePoints();
        }
        calculateAlphaSolid();
        emesh_.mesh().update();
        emesh_.updateEulerMeshInfo();
    }
}

Foam::volScalarField& Foam::newton::alphaSolid()
{
    if (!alphaSolidPtr_.valid())
    {
        alphaSolidPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "alphaSolid",
                    emesh_.mesh().time().timeName(),
                    emesh_.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                emesh_.mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }
    return alphaSolidPtr_();
}

void Foam::newton::calculateAlphaSolid()
{
    volScalarField& alphaSolid = this->alphaSolid();
    alphaSolid = dimensionedScalar("zero", dimless, 0.0);

    forAll(ibo_, objI)
    {
        const labelList& solidCells = ibo_[objI].solidCellsExt();
        forAll(solidCells, cellI)
        {
            const label cellID = solidCells[cellI];
            alphaSolid[cellID] = volFraction(ibo_[objI], cellID);
        }
    }
    alphaSolid.correctBoundaryConditions();
}

Foam::scalar Foam::newton::volFraction(IBObject& ibobj, label cellID)
{
    const faceList& ff = emesh_.mesh().faces();
    const pointField& pp = emesh_.mesh().points();
    const cell& cc = emesh_.mesh().cells()[cellID];
    pointField cellVertices = cc.points(ff, pp);

    IBParticle& ibp = refCast<IBParticle>(ibobj);
    
    scalar alphaIJK(0.0);
    scalar sumPhi(0.0);
    
    forAll(cellVertices, pointI)
    {
        scalar phi_m = LevelSetFunc(ibp.R(), cellVertices[pointI], ibp.CG());
        alphaIJK += -phi_m * HeavisideFunc(-phi_m);
        sumPhi += mag(phi_m);
    }

    return (alphaIJK/sumPhi);
}

Foam::scalar Foam::newton::HeavisideFunc(scalar phi)
{
    if (phi<=0) 
	{
		return 0;
	}
	else
		return 1;
}

Foam::scalar Foam::newton::LevelSetFunc
(
	scalar r, 
	point p, 
	point c
)
{
	if(emesh_.mesh().nGeometricD() == 2)
	{
	    return 
	    (
	    	sqrt
			(
				pow(p.x()-c.x(),2) 
			  + pow(p.y()-c.y(),2)
			)/r - 1.0
	    );
	}
	else 
	{
		return ( mag(p-c)/r - 1.0);
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
