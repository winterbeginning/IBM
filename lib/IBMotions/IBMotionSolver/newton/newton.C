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

        //- find new nei cells and solid cells
        emesh_.updateEulerMeshInfo();
        if(emesh_.mesh().nGeometricD() == 2)
        {
            ibo_[i].createParticle2D();
        }
        else
        {
            ibo_[i].createParticle3D();
        }
        ibo_[i].findNeiCells();
        ibo_[i].findSolidCells();
        ibo_[i].findSolidCellsExt();
    }
}

Foam::volScalarField Foam::newton::calculateAlphaSolid()
{
    Info << "Calculating meshRefine" << endl;
    volScalarField meshRefine_
    (
        IOobject
        (
            "meshRefine",
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
        forAll(emesh_.mesh().C(), cellI)
        {
            meshRefine_[cellI] = volFraction(ibo_[objI], cellI);
        }
    }
    return meshRefine_;
}

// Foam::scalar Foam::newton::volFraction(IBObject& ibobj, label cellID)
// {
//     const faceList& ff = emesh_.mesh().faces();
//     const pointField& pp = emesh_.mesh().points();
//     const cell& cc = emesh_.mesh().cells()[cellID];
//     pointField cellVertices = cc.points(ff, pp);

//     IBParticle& ibp = refCast<IBParticle>(ibobj);
    
//     List<scalar> phi_m(cellVertices.size(), 0.0);

//     forAll(cellVertices, pointI)
//     {
//         phi_m[pointI] = LevelSetFunc(ibp.R(), cellVertices[pointI], ibp.CG());
//         if (phi_m[pointI] < 0)
//         {
//             return 0;
//         }
//     }

//     return (alphaIJK/sumPhi);
// }

// Foam::scalar Foam::newton::volFraction(IBObject& ibobj, label cellID)
// {
//     const faceList& ff = emesh_.mesh().faces();
//     const pointField& pp = emesh_.mesh().points();
//     const cell& cc = emesh_.mesh().cells()[cellID];
//     pointField cellVertices = cc.points(ff, pp);

//     IBParticle& ibp = refCast<IBParticle>(ibobj);
    
//     bool hasInside = false;
//     bool hasOutside = false;

//     // 检查顶点是否在物体内部/外部
//     forAll(cellVertices, pointI)
//     {
//         scalar phi_m = LevelSetFunc(ibp.R(), cellVertices[pointI], ibp.CG());
        
//         if (phi_m <= 0) 
//         {
//             hasInside = true;
//         }
//         else 
//         {
//             hasOutside = true;
//         }
//     }
    
//     // 判断边界情况
//     if (hasInside && hasOutside)
//     {
//         // 边界单元：根据距离比例估算过渡值
//         return 0.5; // 0.3-0.7范围
//     }
//     else if (hasInside)
//     {
//         return 1.0; // 完全在内部
//     }
//     else
//     {
//         return 0.0; // 完全在外部
//     }
// }

Foam::scalar Foam::newton::volFraction(IBObject& ibobj, label cellID)
{
    // 获取单元中心坐标
    const point& cellCenter = emesh_.mesh().C()[cellID];
    
    IBParticle& ibp = refCast<IBParticle>(ibobj);
    
    // 计算到球体表面的距离（带符号）
    scalar signedDist = LevelSetFunc(ibp.R(), cellCenter, ibp.CG()) * ibp.R();  //偏向外侧加密
    
    // 获取网格特征尺寸
    scalar h = emesh_.h();
    
    // 如果单元在边界2h范围内
    if (signedDist <= 5.0 * h && signedDist >= -5.0 * h)
    {
        return 0.5;
    }
    // 如果单元在物体内部
    else if (signedDist <= 0)
    {
        return 1.0;
    }
    // 外部单元
    else
    {
        return 0.0;
    }
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
