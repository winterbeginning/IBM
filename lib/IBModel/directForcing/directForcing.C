/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       directForcing class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#include "directForcing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(directForcing, 0);
	addToRunTimeSelectionTable(IBModel, directForcing, dictionary);
}

//---------------------------------Constructors------------------------------//

Foam::directForcing::directForcing
(
    eulerMesh& emesh,
    const dictionary& dict,
    PtrList<IBObject>& ibo
)
:	
	IBModel(emesh, dict, ibo),
    ibOStream_(emesh.mesh()),
    emesh_(emesh),
    mesh_(emesh.mesh()),
    objects_(ibo),
	nMDF_(readScalar(dict.lookup("multiDirForcingIter")))
{
    uLagrang_.setSize(objects_.size());
    fLagrang_.setSize(objects_.size());
    forAll(objects_, objI)
    {
        uLagrang_[objI].setSize(objects_[objI].nPoints(), vector::zero);
        fLagrang_[objI].setSize(objects_[objI].nPoints(), vector::zero);
    }
}

// -------------------------------Member Functions----------------------------//

Foam::volVectorField Foam::directForcing::ibForce(const volVectorField& U)
{
    Info<< "IBM: Calculating ibForce ..."<<endl;
    
    volVectorField ibForce_
    (
        IOobject
        (
            "ibForce",
            mesh_.time().timeName(),
            mesh_,
            // IOobject::MUST_READ,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("ibForce",
                dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0))
    );   

    if (objects_.size() != 0)
    {
        const scalar dT = mesh_.time().deltaTValue();
        for(int objI=0; objI<objects_.size(); objI++)
        {
            IBObject& objectI = objects_[objI];
            const vectorField& Ub = objectI.uBoundary();

            vectorField ULagr(objectI.nPoints(), vector::zero);
            vectorField FLagr(objectI.nPoints(), vector::zero);
            
            forAll(objectI.lPoints(), pointI)
            {

                //- Interpolating Velocity at Lagragian points
                forAll(objectI.neiCells()[pointI], cellI)
                {
                    scalar deltaFuncValue 
                    = 
                        deltaFunc
                        (
                            mesh_.C()[objectI.neiCells()[pointI][cellI]], 
                            objectI.lPoints()[pointI]
                        );
                    
                    ULagr[pointI] 
                    += 
                        U[objectI.neiCells()[pointI][cellI]] 
                      * deltaFuncValue 
                      * emesh_.dV();
                }

                if (mesh_.nGeometricD() == 2)
                {
                    ULagr[pointI].z() = 0;
                }

                reduce(ULagr[pointI], sumOp<vector>());

            //- Calculate Force at Lagrang points
                FLagr[pointI] = ( Ub[pointI] - ULagr[pointI] ) / dT;

            //- Spread Force from Lagrang points to Euler cells.
                forAll(objectI.neiCells()[pointI], cellI)
                {
                    scalar deltaFuncValue 
                    = 
                        deltaFunc
                        (
                            mesh_.C()[objectI.neiCells()[pointI][cellI]], 
                            objectI.lPoints()[pointI]
                        );
                    
                    ibForce_[objectI.neiCells()[pointI][cellI]] 
                    += FLagr[pointI] * deltaFuncValue * objectI.dVL();
                }
            }
            uLagrang_[objI] = ULagr;
            fLagrang_[objI] = FLagr;
        }
    }   
	return ibForce_;
}

Foam::volVectorField Foam::directForcing::ibForceInt()
{
    Info<< "IBM: Calculating ibForceInt..."<<endl;
    
    // // 修正1: 正确创建并注册 rho 场
    // if (!mesh_.foundObject<volScalarField>("rho"))
    // {
    //     const scalar rhoValue = objects_[0].rho();
        
    //     volScalarField* rhoPtr = new volScalarField
    //     (
    //         IOobject
    //         (
    //             "rho",
    //             mesh_.time().timeName(),
    //             mesh_,
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //         ),
    //         mesh_,
    //         dimensionedScalar("rho", dimDensity, rhoValue)
    //     );
        
    //     rhoPtr->store();  // 关键：注册到数据库
    // }
       
    volVectorField ibForceInt_
    (
        IOobject
        (
            "ibForceInt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("ibForceInt",
                //dimensionSet(1,-2,-2,0,0,0,0), vector(0,0,0))
                dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0))
    );   
    volScalarField TSolidCells
    (
        IOobject
        (
            "TSolidCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("TSolidCells", dimless, 0)
    );

    const scalar dT = mesh_.time().deltaTValue();
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    //const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    for(int objI=0; objI< objects_.size(); objI++)
    {
        const labelList& solidCells_ = objects_[objI].solidCells();
        const vector US = objects_[objI].uTranslate();

        forAll(solidCells_, cellI)
        {
            //ibForceInt_[solidCells_[cellI]] = rho[solidCells_[cellI]]*(US - U[solidCells_[cellI]])/dT;
            ibForceInt_[solidCells_[cellI]] = (US - U[solidCells_[cellI]])/dT;
            TSolidCells[solidCells_[cellI]] = 1;
        }
    }
    
    if (mesh_.time().outputTime())
    {
        TSolidCells.write();
    }

    return ibForceInt_;
}

void Foam::directForcing::multiDirectForcing
(
    volVectorField& u,
    volVectorField& ibForce_
)
{
    if(nMDF_ > 0)
    {
        dimensionedScalar dT("dT",dimTime,mesh_.time().deltaTValue());;

        for(int i=0; i<nMDF_; i++)
        {
            Info<< "IBM: Multi-direct Forcing Iteration "<<i+1<<endl;
            volVectorField f = ibForce(u);

            u += dT*f;
        }
    }
}


// void Foam::directForcing::update()
// {
//     //- Calculate velocity of Lagrang points, move points accordingly, 
//     //  and update neighbour cells
//    if (mesh_.time().timeIndex() > 0)
//    {
//        moveObjects(fLagrang_);
//    }

//     writeNeighbourCells();
//     // //- Refine mesh following neighbour cells indicator
//     // mesh_.update(cellsToRefine());

//     // //- Update neighbour cells after mesh refinement
//     // for(int i=0; i<size(); i++)
//     // {
//     //     neiCells()[i] = findNeiCells(LPoints()[i]);
//     // }

//     // //- Update h and dV (only do 1 time at timeIndex = 1)
//     // updateCartesianGridSize();

//     // //- Write neighbour cells indicator
//     // writeNeighbourCells();

// }

void Foam::directForcing::write()
{
    for(int objI=0; objI<objects_.size(); objI++)
    {
        ibOStream_.writeLagrPoints(objI, objects_[objI].lPoints());
        ibOStream_.writeLagrForces(objI, objects_[objI].lPoints(), fLagrang_[objI]);
        ibOStream_.writeIBForces(objI, fLagrang_[objI], emesh_.dV());
        ibOStream_.writeObjData
        (
            objI, 
            objects_[objI].CG(), 
            objects_[objI].uTranslate(),
            objects_[objI].uRotate()
        );
        ibOStream_.writeObjVTU
        (
            objI, 
            objects_[objI].lPoints(), 
            objects_[objI].CG(), 
            objects_[objI].nFaces(), 
            objects_[objI].nPointsOfFaces(), 
            objects_[objI].pointOfFace()
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
