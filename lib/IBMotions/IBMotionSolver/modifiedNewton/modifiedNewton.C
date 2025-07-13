/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       modifiedNewton class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#include "modifiedNewton.H"
#include "wallPolyPatch.H"
#include "IBParticle.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(modifiedNewton, 0);
    addToRunTimeSelectionTable(IBMotionSolver, modifiedNewton, dictionary);
}

//---------------------------------Constructors------------------------------//
Foam::modifiedNewton::modifiedNewton
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
{}

// -------------------------------Member Functions----------------------------//

void  Foam::modifiedNewton::moveObjects()
{
    for(int i=0; i<ibo_.size(); i++)
    {
        if (ibo_[i].movable())
        {
            Info<< "Moving object: "<<ibo_[i].name()<<endl;
            
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
            
            const volVectorField& U = emesh_.mesh().lookupObject<volVectorField>("U");
            const volVectorField& Uold = U.oldTime();
            const scalar dT = emesh_.mesh().time().deltaTValue(); 
            vector volIntegralU_ 
            = 
                (
                    volIntegralU(ibo_[i], U) 
                  - volIntegralU(ibo_[i], Uold)
                ) / dT;

            vector volIntegralRxU_ 
            = 
                (
                    volIntegralRxU(ibo_[i], U) 
                  - volIntegralRxU(ibo_[i], Uold)
                ) / dT;
        
            ibo_[i].updateObjectMotionTobias
            (
                ibmodelPtr_->FLagrange()[i],
                repulsiveForce,
                emesh_.rhoFluid(),
                emesh_.g(),
                volIntegralU_,
                volIntegralRxU_
            );

            ibo_[i].updateVelocity();
            ibo_[i].movePoints();
        }
        //- find new nei cells and solid cells
        ibo_[i].findNeiCells();
        ibo_[i].findSolidCells();
    }
}

Foam::volScalarField Foam::modifiedNewton::calculateAlphaSolid()
{
    Info << "Calculating alphaSolid" << endl;
    volScalarField alphaSolid_
    (
        IOobject
        (
            "alphaSolid",
            emesh_.mesh().time().timeName(),
            emesh_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        emesh_.mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    );

    forAll(ibo_, objI)
    {
        const labelList& solidCells = ibo_[objI].solidCells();
        forAll(solidCells, cellI)
        {
            const label cellID = solidCells[cellI];
            alphaSolid_[cellID] = volFraction(ibo_[objI], cellID);
        }
    }
    return alphaSolid_;
}

Foam::vector Foam::modifiedNewton::volIntegralU
(
    IBObject& ibobj,
    const volVectorField& U
)
{   
    vector IU(vector::zero);

    forAll(ibobj.solidCellsExt(), cellI)
    {
        scalar vf = volFraction(ibobj, ibobj.solidCellsExt()[cellI]);

        IU += U[ibobj.solidCellsExt()[cellI]] * emesh_.dV() * vf;
    }

    return IU;
}

Foam::vector Foam::modifiedNewton::volIntegralRxU
(
    IBObject& ibobj,
    const volVectorField& U
)
{
    const volVectorField& cc = emesh_.mesh().C();

    vector IRxU(vector::zero);

    forAll(ibobj.solidCellsExt(), cellI)
    {
        scalar vf = volFraction(ibobj, ibobj.solidCellsExt()[cellI]);

        vector r = cc[ibobj.solidCellsExt()[cellI]] - ibobj.CG();

        IRxU +=
        ( 
            (r ^ U[ibobj.solidCellsExt()[cellI]]) * emesh_.dV() * vf
        );
    } 

    return IRxU;
}

Foam::scalar Foam::modifiedNewton::volFraction(IBObject& ibobj, label cellID)
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

Foam::scalar Foam::modifiedNewton::HeavisideFunc(scalar phi)
{
    if (phi<=0) 
	{
		return 0;
	}
	else
		return 1;
}

Foam::scalar Foam::modifiedNewton::LevelSetFunc
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
