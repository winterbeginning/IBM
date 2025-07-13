/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBParticle class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#include "IBParticle.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "IBMotion.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "IBMotions/IBMotion/sixDoFMotion/sixDoFMotion.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBParticle, 0);
    addToRunTimeSelectionTable(IBObject, IBParticle, dictionary);
}

// * * * * * * * * * * * * Private Member Fuctions * * * * * * * * * * * * * //
void Foam::IBParticle::initialize(const dictionary& dict)
{
    //- Read input data
	center_  = dict.lookup("CG");
	R_       = readScalar(dict.lookup("radius"));
	rho_     = readScalar(dict.lookup("rho"));
    movable_ = dict.lookupOrDefault<Switch>("movable", false);
    calcWallRepulsive_ = dict.lookupOrDefault<Switch>("calcWallRepulsive", false);
    autoPoint_ = dict.lookupOrDefault<Switch>("autoPoint", true);
    nPoints_ = dict.lookupOrDefault<label>("nPoints", 100);

    if (calcWallRepulsive_)
    {
        epsilonW_ = readScalar(dict.lookup("epsilonW"));
        epsilonP_ = readScalar(dict.lookup("epsilonP"));
    }

    if(emesh_.mesh().nGeometricD() == 2)
    {
        V_ = PI*R_*R_;

        Ip_.x() = 0.0;
        Ip_.y() = 0.0;
        Ip_.z() = rho_*V_*R_*R_/2.0;

        createParticle2D();
    }
    else
    {
        V_ = 4.0/3.0*PI*R_*R_*R_;

        Ip_.x() = 2.0/5.0*rho_*V_*R_*R_;    
        Ip_.y() = 2.0/5.0*rho_*V_*R_*R_;    
        Ip_.z() = 2.0/5.0*rho_*V_*R_*R_;

        createParticle3D();
    }

    addMotions(dict);
    findNeiCells();
    findSolidCells();
    findSolidCellsExt();
}

void Foam::IBParticle::createParticle2D()
{ 
    Info<<"------------------"<<endl;
    Info<<"Create 2D particle"<<endl;
    Info<<"------------------"<<endl;
    
    calcNLagrangPoints();

    Info<<"  R = "<<R_<<"; C = "<<center_<<"; dVL = "<<dVL_<<endl;

    lPoints_.setSize(nPoints_, vector::zero);

    scalar dAlpha = 2.0*PI/nPoints_;
    forAll(lPoints_, i)
    {
        lPoints_[i] = point (   center_.x() + R_*Foam::cos(dAlpha*i),
                                center_.y() + R_*Foam::sin(dAlpha*i),
                                center_.z()
                            );
    }

    nFaces_ = nPoints_;
    nPointsOfFaces_.setSize(nFaces_);
    pointOfFace_.setSize(nFaces_);
    for(int i=0;i<nFaces_;i++)
    {
        nPointsOfFaces_[i] = 3;
        pointOfFace_[i].setSize(3);
        pointOfFace_[i][0] = i;
        pointOfFace_[i][1] = (i+1)%nPoints_;
        pointOfFace_[i][2] = nPoints_;
    }
}

void Foam::IBParticle::createParticle3D()
{
    Info<<"------------------"<<endl;
    Info<<"Create 3D particle"<<endl;
    Info<<"------------------"<<endl;

    calcNLagrangPoints();

    Info<<" R = "<<R_<<"; C = "<<center_<<"; dVL = "<<dVL_<<endl;

    //- Distribute nPoints points on surface of unit sphere
    lPoints_ = createUnitSphereEqAreaPartition(nPoints_);
    
    //- Translate and scale unit sphere according to center_ and R_
    //  to obtain desired 3D sphere
    for (int pointI=0;pointI<nPoints_;pointI++)
    {
        lPoints_[pointI].x()=R_*lPoints_[pointI].x()+center_.x();
        lPoints_[pointI].y()=R_*lPoints_[pointI].y()+center_.y();
        lPoints_[pointI].z()=R_*lPoints_[pointI].z()+center_.z();
    }

    nFaces_ = 2*nPoints_ - 4;
    nPointsOfFaces_.setSize(nFaces_, 3);
    pointOfFace_.setSize(nFaces_);
    for(int i=0; i<nFaces_; i++)
    {
        pointOfFace_[i].setSize(nPointsOfFaces_[i]);

    }
}

void Foam::IBParticle::calcNLagrangPoints()
{
    scalar idealNPoints(0.0);

    if (emesh_.mesh().nGeometricD() == 2)
    {
        idealNPoints = 2.0*PI*R_/emesh_.h();
        if(autoPoint_)
            nPoints_ = static_cast<int> (idealNPoints + 0.5);
        
        dVL_ = 2.0*PI*R_ / nPoints_ * emesh_.h();
    }
    else
    {
        idealNPoints = PI/3.0*(12.0*pow(R_/emesh_.h(), 2.0) + 1);
        if(autoPoint_)
            nPoints_ = static_cast<int> (idealNPoints + 0.5);
        dVL_ = 4.0*PI*R_*R_ / nPoints_ * emesh_.h();
    }
    

    Info<<"  Ideal number of Lagrangian points : "<<idealNPoints<<nl
        <<"  Actual number of Lagrangian points: "<<nPoints_<<endl;;

}

Foam::pointField Foam::IBParticle::createUnitSphereEqAreaPartition(label nPoints)
{
    //- Paul Leopardi algorithm 2007
    //  Evenly distribute nPoints points on a unit sphere

    //- 1. Determine colatitudes of polar caps

        //- Area of a region
        scalar Sr = 4.0*PI/nPoints;
        
        //- Collatitude of the North collar cap.
        //  Calculate form area of the cap: Sr = 4*PI*[sin(phiNC/2)]^2
        scalar phiNC = 2.0*asin(sqrt(Sr/(4.0*PI)));
    
    //- 2. Determine an ideal collar angle
    //  (when number of points reach infinity)
       
        scalar deltaI = pow(Sr, 0.5);
    //- 3. Determine an ideal number of collars

        scalar nI = (PI-2.0*phiNC) / deltaI;

    //- 4. Determine the actual number of collars
       
        label n = static_cast<int>(nI+0.5);

    //- 5. Create a list of the ideal number of regions in each collar
    
        //- Fitting collar angle
        scalar deltaF = (PI - 2.0*phiNC)/n;

        //- List of colatitudes of caps
        List<scalar> colatitudes(n+2, 0.0);
        
        forAll(colatitudes, i)
        {
            if (i == 0)
                colatitudes[i] = phiNC;
            else if( i != 0 && i < n+1 )
                colatitudes[i] = phiNC + i*deltaF;
            else
                colatitudes[i] = PI;    
        }

        //- Ideal number of region in each collar
        List<scalar> idealRegionsPerCollar(n+2, 0.0);
        
        forAll(idealRegionsPerCollar, i)
        {
            if (i==0 || i==n+1)
            {
                idealRegionsPerCollar[i] = 1;
            }
            else 
            {
                idealRegionsPerCollar[i] 
                = 
                    (capArea(colatitudes[i]) - capArea(colatitudes[i-1]))/Sr;
            }
        }

    //- 6. Createa a list of actial number of regions in each collar
        
        List<label> actualRegionsPerCollar(n+2,0);
        
        forAll(actualRegionsPerCollar, i)
        {
            if (i==0 || i==n+1)
            {
                actualRegionsPerCollar[i] = 1;
            }
            else if (i == 1)
            {
                actualRegionsPerCollar[i] = static_cast<int>(idealRegionsPerCollar[i]);
            }
            else 
            {
                scalar ai(0.0);
                for(int j=0; j<i; j++)
                {
                    ai += (idealRegionsPerCollar[j] - actualRegionsPerCollar[j]) ;
                } 

                actualRegionsPerCollar[i] = static_cast<int>(idealRegionsPerCollar[i]+ai);
            }
        }

    //- 7. Create a list of colatitudes of points on each collar
        List<List<scalar>> pointLatitudes(n+2);
        List<List<scalar>> pointLongtitudes(n+2);

        for(int i=0; i<n+2; i++)
        {
            if (i==0)
            {
                pointLatitudes[i].setSize(1);
                pointLongtitudes[i].setSize(1);
                forAll(pointLatitudes[i], j)
                {
                    pointLatitudes[i][j] = 0;
                    pointLongtitudes[i][j] = 0;
                }
            }
            else if ( i!=0 && i<n+1)
            {
                pointLatitudes[i].setSize(actualRegionsPerCollar[i]);
                pointLongtitudes[i].setSize(actualRegionsPerCollar[i]);

                forAll(pointLatitudes[i], j)
                {
                    pointLatitudes[i][j] = colatitudes[i-1] + 0.5*(colatitudes[i] - colatitudes[i-1]);
                    pointLongtitudes[i][j] = j*2.0*PI/actualRegionsPerCollar[i];
                }

            }
            else 
            {
                pointLatitudes[i].setSize(1);
                pointLongtitudes[i].setSize(1);
                forAll(pointLatitudes[i], j)
                {
                    pointLatitudes[i][j] = PI;
                    pointLongtitudes[i][j] = 0;
                }
            }
        }

    //8. Create list of points
        pointField distributedPoints;

        for(int i=0; i<n+2; i++)
        {
            forAll(pointLatitudes[i], j)
            {
                point p;
                p.x() = sin(pointLatitudes[i][j])*cos(pointLongtitudes[i][j]);
                p.y() = sin(pointLatitudes[i][j])*sin(pointLongtitudes[i][j]); 
                p.z() = cos(pointLatitudes[i][j]);
                distributedPoints.append(p);
            }
        }

    //9. prepare to save sphere
        nFaces_ = 2*nPoints - 4;
        nPointsOfFaces_.setSize(nFaces_, 3);
        pointOfFace_.setSize(nFaces_);
        for(int i=0; i<nFaces_; i++)
        {
            pointOfFace_[i].setSize(nPointsOfFaces_[i]);
        }

    return distributedPoints;
}

Foam::scalar Foam::IBParticle::capArea(scalar angle)
{
    return (4.0*PI*pow(sin(angle/2.0), 2));
}

void Foam::IBParticle::addMotions(const dictionary& dict)
{
    if(movable_)
    {
        const dictionary& motionDict = dict.subDict("motions");
        
        label i=0;
        motions_.setSize(motionDict.size());
        forAllConstIter(IDLList<entry>, motionDict, iter)
        {
            if(iter().isDict())
            {
                motions_.set
                (
                    i++,
                    IBMotion::New 
                    (
                        iter().keyword(),
                        *this,
                        iter().dict()
                    )
                );
                Info<<"  Motion "<<i<<": "<<iter().keyword()<<endl;
            }
        }
        motions_.setSize(i);
    }
    else
    {
        Info<< setw(11) <<"  Motions: N/A"<<nl<<endl;
    }
}
//--------------------------------Constructors-------------------------------//

Foam::IBParticle::IBParticle
(	
	const word& typeName,
    eulerMesh& emesh,
    const dictionary& dict
)
:
	IBObject(typeName, emesh, dict),
	PI(3.14159265359),
    emesh_(emesh),
    name_(dict.lookup("name")),
    objectType_(typeName),
    center_(point::zero),
    R_(0.0),
    rho_(0.0),
    V_(0.0),
    Ip_(vector::zero),
    nPoints_(0),   
	lPoints_(),
    dVL_(),
    neiCells_(),
    solidCells_(),
    solidCellsExt_(),
    movable_(),
    motions_(),
    uTranslate_(vector::zero),
    uRotate_(vector::zero),
    uBoundary_(),
	nFaces_(0),
	nPointsOfFaces_(),
	pointOfFace_(),
    calcWallRepulsive_(),
    epsilonW_(0.0),
    epsilonP_(0.0)
{
    initialize(dict);
    updateVelocity();
}

//--------------------------------Member Functions---------------------------//
void Foam::IBParticle::findNeiCells()
{
    neiCells_.clear();

    neiCells_.setSize(nPoints_);

    volScalarField neighbourCells
    (
        IOobject
        (
            "neighbourCells",
            emesh_.mesh().time().timeName(),
            emesh_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        emesh_.mesh(),
        dimensionedScalar("neighbourCells", dimless, 0)
    );

    forAll(lPoints_, pointI)
    {
        forAll(emesh_.mesh().C(), cellI)
        { 
            scalar cellToLpoint = mag(emesh_.mesh().C()[cellI] - lPoints_[pointI]);
            
            scalar span = 2.0*emesh_.h();
            
            if (cellToLpoint <= span)
            {
                neiCells_[pointI].append(cellI);
                neighbourCells[cellI] = 1;
            }
        }
    }

    if (emesh_.mesh().time().timeIndex() == 0 || emesh_.mesh().time().outputTime())
    {
        neighbourCells.write();
    }
}

void Foam::IBParticle::findSolidCells()
{
    solidCells_.clear();

    volScalarField sldCells
    (
        IOobject
        (
            "sldCells",
            emesh_.mesh().time().timeName(),
            emesh_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        emesh_.mesh(),
        dimensionedScalar("sldCells", dimless, 0)
    );

    forAll(emesh_.mesh().C(), cellI)
    { 
        scalar dR = mag(emesh_.mesh().C()[cellI] - center_);
        //- not very effective!
        
        if (dR <= R_)
        {
            solidCells_.append(cellI);
            sldCells[cellI] = 1;
        }
    }

    if (emesh_.mesh().time().timeIndex() == 0 || emesh_.mesh().time().outputTime())
    {
        sldCells.write();
    }

}

void Foam::IBParticle::findSolidCellsExt()
{
    solidCellsExt_.clear();

    volScalarField sldCellsExt
    (
        IOobject
        (
            "sldCellsExt",
            emesh_.mesh().time().timeName(),
            emesh_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        emesh_.mesh(),
        dimensionedScalar("sldCellsExt", dimless, 0)
    );

    forAll(emesh_.mesh().C(), cellI)
    { 
        scalar dR = mag(emesh_.mesh().C()[cellI] - center_);
        //- not very effective!
        
        if (dR <= (R_ + 2.0*emesh_.h()))
        {
            solidCellsExt_.append(cellI);
            sldCellsExt[cellI] = 1;
        }
    }

    if (emesh_.mesh().time().timeIndex() == 0 || emesh_.mesh().time().outputTime())
    {
        sldCellsExt.write();
    }

}

const Foam::word& Foam::IBParticle::name()
{
    return name_;
}

const Foam::word& Foam::IBParticle::objectType() 
{
    return objectType_;
}

const Foam::scalar& Foam::IBParticle::rho() 
{
    return rho_;
}

const Foam::scalar& Foam::IBParticle::V() 
{
    return V_;
}

const Foam::vector& Foam::IBParticle::Ip() 
{
    return Ip_;
}

const Foam::scalar& Foam::IBParticle::R() 
{
    return R_;
}

Foam::point& Foam::IBParticle::CG() 
{
    return center_;
}

Foam::label& Foam::IBParticle::nPoints() 
{
    return nPoints_;
}

Foam::pointField& Foam::IBParticle::lPoints() 
{
    return lPoints_;
}

const Foam::scalar& Foam::IBParticle::dVL()
{
    return dVL_;
}

const Foam::labelListList& Foam::IBParticle::neiCells() 
{
    return neiCells_;
}

const Foam::labelList& Foam::IBParticle::solidCells() 
{
    return solidCells_;
}

const Foam::labelList& Foam::IBParticle::solidCellsExt() 
{
    return solidCellsExt_;
}

bool Foam::IBParticle::movable() 
{
    return movable_;
}

Foam::PtrList<Foam::IBMotion>& Foam::IBParticle::motions() 
{
    return motions_;
}

const Foam::vector& Foam::IBParticle::uTranslate()
{
    return uTranslate_;
}

const Foam::vector& Foam::IBParticle::uRotate()
{
    return uRotate_;
}

const Foam::vectorField& Foam::IBParticle::uBoundary()
{
    return uBoundary_;
}

const Foam::label& Foam::IBParticle::nFaces()
{
    return nFaces_;
}

const Foam::labelList& Foam::IBParticle::nPointsOfFaces()
{
    return nPointsOfFaces_;
}

const Foam::labelListList& Foam::IBParticle::pointOfFace()
{
    return pointOfFace_;
}

bool Foam::IBParticle::calcWallRepulsive()
{
    return calcWallRepulsive_;
}

const Foam::scalar& Foam::IBParticle::epsilonW()
{
    return epsilonW_;
}

const Foam::scalar& Foam::IBParticle::epsilonP()
{
    return epsilonP_;
}

void Foam::IBParticle::updateObjectMotionUhlmann
(
    const vectorField& ForceLagrang,
    vector& repulsiveForce,
    const scalar& rhoFluid,
    const vector& g

)
{
    repulsiveForce += wallRepulsiveForce();

    forAll(motions_, i)
    {
        if (motions_[i].motionType() == "sixDoFMotion")
        {
            sixDoFMotion& sDoF = refCast<sixDoFMotion>(motions_[i]);
            sDoF.updateMotion
            (
                emesh_.mesh(),
                ForceLagrang,
                repulsiveForce,
                rhoFluid,
                g
            );
        }
    }
}

void Foam::IBParticle::updateObjectMotionTobias
(
    const vectorField& ForceLagrang,
    vector& repulsiveForce,
    const scalar& rhoFluid,
    const vector& g,
    const vector& volIntegralU,
    const vector& volIntegralRxU

)
{
    repulsiveForce += wallRepulsiveForce();

    forAll(motions_, i)
    {
        if (motions_[i].motionType() == "sixDoFMotion")
        {
            sixDoFMotion& sDoF = refCast<sixDoFMotion>(motions_[i]);
            sDoF.updateMotion
            (
                emesh_.mesh(),
                ForceLagrang,
                repulsiveForce,
                rhoFluid,
                g,
                volIntegralU,
                volIntegralRxU
            );
        }
    }
}

void Foam::IBParticle::updateVelocity()
{
    uBoundary_.setSize(nPoints_, vector::zero);
    if (movable_)
    {
        //- Reset velocities
        uTranslate_ = vector::zero;
        uRotate_ = vector::zero;

        forAll(motions_, i)
        {
            uTranslate_ += motions_[i].translationalVelocity();
            uRotate_ += motions_[i].rotationalVelocity();
        }
        forAll(uBoundary_, i)
        {
            uBoundary_[i] 
            = 
                uTranslate_ 
              + (uRotate_^(lPoints_[i] - center_));
        }
        Info<<"  Update motion for IB object "<<name_<<nl
            <<"    - Translate: "<<uTranslate_<<nl
            <<"    - Rotate   : "<<uRotate_<<endl;
    }
}

void Foam::IBParticle::movePoints()
{
    if (movable_)
    {
        
        const scalar dT = emesh_.mesh().time().deltaTValue();

        center_ += dT*uTranslate_;

        forAll(lPoints_, i)
        {
            lPoints_[i] += dT*uBoundary_[i];
        }

        Info<<"    Move IB object "<<name_
            <<" from "<<(center_-dT*uTranslate_)
            <<" to "<<center_<<nl<<endl;
    }

}

const Foam::vector Foam::IBParticle::wallRepulsiveForce()
{
    //- Calculate particle-wall interaction
    //  (Wan, Turek 2007)

	vector Fwr(vector::zero);

    if(calcWallRepulsive_)
    {
        forAll(emesh_.mesh().boundary(), patchI)
        {
            const polyPatch& pp = emesh_.mesh().boundaryMesh()[patchI];
            if(emesh_.mesh().boundary()[patchI].Cf().size() == 0)
            {
                continue;
            }
            else if(isA<wallPolyPatch>(pp))
            {
                const vectorField fc = pp.faceCentres();
                    
                //- Find face closest to C
                scalar d(GREAT);
                label nearestFace(0);
                forAll(fc, facei)
                {
                    scalar centerToFace = mag(center_ - fc[facei]);
                    if (centerToFace <= d)
                    {
                        d = centerToFace;
                        nearestFace = facei;
                    }
                }

                //- Find projection of C on the nearest face
                //  **projection = p - ((p-c)*n)*n
                    const vector fcn = fc[nearestFace];

                //- Surface area vector of the nearest face           
                    vector Sf_ = emesh_.mesh().boundary()[patchI].Sf()[nearestFace];

                //- Normal vector of the nearest face
                    vector nf_ = Sf_/mag(Sf_);

                //- Vector from particle's center to nearest face's center
                    vector p_fcn = center_ - fcn;

                //- The nearest distance
                    scalar distance = p_fcn & nf_;
   
                //- the projection point
                point c_proj = center_ - distance*nf_;

                scalar di = 2.0*mag(center_ - c_proj);
                scalar xi = 2.0*emesh_.h();
                scalar epsilonWW = epsilonW_;
                //- Theoretical value of stiffness values:

                if (di <= 2.0*R_)
                {
                    Fwr += 2.0*(center_ - c_proj)*(2.0*R_ - di)/epsilonW_;
                }
                else if ((di > 2.0*R_) && (di <= 2.0*R_ + xi))
                {
                    Fwr 
                    += 
                        2.0*(center_ - c_proj)
                       *Foam::pow(2.0*R_ + xi - di, 2.0)/epsilonWW;
                }
                else 
                {
                    Fwr += vector::zero;
                }
            }
            else
                continue;
        }
    }

    return Fwr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //