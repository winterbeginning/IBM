/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       eulerMesh class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#include "eulerMesh.H"
#include "uniformDimensionedFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(eulerMesh, 0);
}

// ------------------------Private Member Functions--------------------------//
void Foam::eulerMesh::readEnvironmentVariables()
{
    if (transportProperties_.found("phases")) //two-phases simulations
    {
        Pair<word> phases = transportProperties_.lookup("phases");
        scalar rho1 = readScalar(transportProperties_.subDict(phases[0]).lookup("rho"));
        scalar rho2 = readScalar(transportProperties_.subDict(phases[1]).lookup("rho"));

        // const dictionary& phases = transportProperties_.subDict("phases");
        // scalar rho1 = readScalar(phases.subDict("phase1").lookup("rho"));
        // scalar rho2 = readScalar(phases.subDict("phase2").lookup("rho"));

        if (rho1 > rho2)
        {
            rhoFluid_ = rho1;
            rhoAir_   = rho2;
        }
        else 
        {
            rhoFluid_ = rho2;
            rhoAir_   = rho1;
        }

        uniformDimensionedVectorField gravity
        (
            IOobject
            (
                "g",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        g_ = gravity.value();
        
        Info<<"Reading rho1,rho2 and g for the two-phase problem"<<nl
            <<"  rhoFluid = "<<rhoFluid_<<nl
            <<"    rhoAir = "<<rhoAir_
            <<"         g = "<<g_<<nl<<endl;
    }
    else // one phase simulation
    {
        rhoFluid_ = transportProperties_.lookupOrDefault("rhoF",1000.0);
        rhoAir_ = rhoFluid_;
        uniformDimensionedVectorField gravity
        (
            IOobject
            (
                "g",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        g_ = gravity.value();
        Info<<"Reading rho and g for the single-phase problem"<<nl
            <<"  rhoFluid = "<<rhoFluid_<<nl
            <<"         g = "<<g_<<nl<<endl;
    }
}

Foam::scalar Foam::eulerMesh::cellSize(label cellID)
{
    scalar delta = 0.0;
    if (mesh_.nGeometricD() == 3)
    {
        delta = pow(mesh_.V().field()[cellID], 1.0/3.0);
    }
    else
    {
        scalar thickness = 0.0;
        Vector<label> directions = mesh_.geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        delta = Foam::sqrt(mesh_.V().field()[cellID]/thickness);
    }

    return delta;
}

//---------------------------------Constructors------------------------------//

Foam::eulerMesh::eulerMesh
(
	dynamicFvMesh& mesh
)
:	
	mesh_(mesh),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
	h_(),
	dV_()
{
    readEnvironmentVariables();
	updateEulerMeshInfo();
}

// -------------------------------Member Functions----------------------------//
void Foam::eulerMesh::updateEulerMeshInfo()
{
	scalarField cellsSize(mesh_.C().size(), 0.0);

    forAll(mesh_.C(), cellI)
    {
        cellsSize[cellI] = cellSize(cellI);
    }
    
    h_ = min(cellsSize);

    if (mesh_.nGeometricD() == 2)
    {
    	dV_ = pow(h_,2);
    }
    else 
    {
    	dV_ = pow(h_, 3);
    }
    Info<<"Cartesian grid size :  h = "<<h_<<nl
        <<"Cartesian volumetric: dV = "<<dV_<<nl<<endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
