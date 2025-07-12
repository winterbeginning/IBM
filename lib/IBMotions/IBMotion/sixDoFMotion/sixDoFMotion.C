/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - sixDoFMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "sixDoFMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFMotion, 0);
    addToRunTimeSelectionTable(IBMotion, sixDoFMotion, dictionary);
}
//--------------------------------Private Members----------------------------//
void Foam::sixDoFMotion::readConstraint(const dictionary& dict)
{
	if (dict.found("constraint"))
	{
		constraintDir_ = dict.subDict("constraint").lookup("axis");
		Info<<"FOUND CONSTRAINT MOTION: "<<constraintDir_<<endl;
	}
}
//---------------------------------Constructors------------------------------//
Foam::sixDoFMotion::sixDoFMotion
(
	const word& typeName,
	IBObject& obj,
	const dictionary& dict
)
:
	IBMotion(typeName, obj, dict),
	obj_(obj),
	uTrans_(vector::zero),
	uRot_(vector::zero),
	constraintDir_(vector::zero)
{
	readConstraint(dict);
}
//-------------------------------Member functions----------------------------//
void Foam::sixDoFMotion::updateMotion
(
	dynamicFvMesh& mesh,
	const vectorField& Fk, 
	vector repulsiveForce,
	scalar rhoF,
	vector g
)
{
	//- M. Uhlmann 2005
	const scalar dT = mesh.time().deltaTValue();  

	vector force_ = vector::zero;
	forAll(Fk, i)
	    force_ += Fk[i]*obj_.dVL();

	uTrans_
	= 
		uTrans_ 
	  + dT*( 
		    - rhoF*force_/(obj_.V()*(obj_.rho()-rhoF)) 
			+ repulsiveForce /(obj_.V()*(obj_.rho() - rhoF))
			+ g
			);

	vector fxR = vector::zero;
	forAll(obj_.lPoints(), pointI)
	{
		vector dr = obj_.lPoints()[pointI] - obj_.CG();
	    fxR += (dr ^ Fk[pointI])*obj_.dVL();  
	}

	if (mesh.nGeometricD() == 2) 
	{
		if (mesh.geometricD().x() < 0)
		{
			uTrans_.x() = 0.0;
			uRot_.y() = 0.0;
			uRot_.z() = 0.0;
			uRot_.x() 
			= 
				uRot_.x() 
			  - dT*obj_.rho()*rhoF*fxR.x()/(obj_.Ip().x()*(obj_.rho()-rhoF));
		}
		else if (mesh.geometricD().y() < 0)
		{
			uTrans_.y() = 0.0;
			uRot_.z() = 0.0;
			uRot_.x() = 0.0;
			uRot_.y() 
			= 
				uRot_.y() 
			  - dT*obj_.rho()*rhoF*fxR.y()/(obj_.Ip().y()*(obj_.rho()-rhoF));
		}
		else 
		{
			uTrans_.z() = 0.0;
			uRot_.x() = 0.0;
			uRot_.y() = 0.0;
			uRot_.z() 
			= 
				uRot_.z() 
			  - dT*obj_.rho()*rhoF*fxR.z()/(obj_.Ip().z()*(obj_.rho()-rhoF));
		}
	}
	else if (mesh.nGeometricD() == 3)
	{
		uRot_.x() 
		= 
			uRot_.x() 
		  - dT*obj_.rho()*rhoF*fxR.x()/(obj_.Ip().x()*(obj_.rho()-rhoF));
		uRot_.y() 
		= 
			uRot_.y() 
		  - dT*obj_.rho()*rhoF*fxR.y()/(obj_.Ip().y()*(obj_.rho()-rhoF));
		uRot_.z() 
		= 
			uRot_.z() 
		  - dT*obj_.rho()*rhoF*fxR.z()/(obj_.Ip().z()*(obj_.rho()-rhoF));
	}
    else
    {
    	FatalErrorIn
    	(
    		"sixDoFMotion::updateMotion(vector& uTransl, "
			"vector& uRotate, const fvMesh& mesh," 
			"const vectorField& Fk, vector repulsiveForce," 
			"point center, pointField lPoints, scalar rhoF, "
			"scalar dV, vector g)"
    	)	<< "The mesh is neither 2D nor 3D " << endl
    		<< "Please check again"<<endl
    		<<exit(FatalError);
    }

	if (constraintDir_ != vector::zero)
	{
		uTrans_ = (uTrans_ & constraintDir_)*constraintDir_;
	}
	
}

void Foam::sixDoFMotion::updateMotion
(
	dynamicFvMesh& mesh,
	const vectorField& Fk, 
	vector repulsiveForce,
	scalar rhoF,
	vector g,
	vector volIntegralU,
	vector volIntegralRxU
)
{
	//- Tobias Kempe & J.Frohlich 2012
	const scalar dT = mesh.time().deltaTValue();  

	vector force_ = vector::zero;
	forAll(Fk, i)
	    force_ += Fk[i]*obj_.dVL();

	uTrans_ 
	=
	  	uTrans_ 
	  + dT/(obj_.V()*obj_.rho())
	  * ( 
	  		- rhoF*force_ 
	  		+ rhoF*volIntegralU 
	  		+ obj_.V()*(obj_.rho()-rhoF)*g 
	  		+ repulsiveForce 
	  	);
	
	vector fxR = vector::zero;
	forAll(obj_.lPoints(), pointI)
	{
		vector r = obj_.lPoints()[pointI] - obj_.CG();
	    fxR += (r ^ Fk[pointI])*obj_.dVL();  
	}

    if (mesh.nGeometricD() == 2) 
    {
    	if (mesh.geometricD().x() < 0)
    	{
    		uTrans_.x() = 0.0;
    		uRot_.y() = 0.0;
    		uRot_.z() = 0.0;
    		uRot_.x() 
    		= 
    			uRot_.x() 
    		  + dT*rhoF/obj_.Ip().x()*( -fxR.x() + volIntegralRxU.x());
    	}
    	else if (mesh.geometricD().y() < 0)
    	{
    		uTrans_.y() = 0.0;
    		uRot_.z() = 0.0;
    		uRot_.x() = 0.0;
    		uRot_.y() 
    		= 
    			uRot_.y() 
    		  + dT*rhoF/obj_.Ip().y()*( -fxR.y() + volIntegralRxU.y());
    	}
    	else 
    	{
    		uTrans_.z() = 0.0;
    		uRot_.x() = 0.0;
    		uRot_.y() = 0.0;
    		uRot_.z() 
    		= 
    			uRot_.z() 
    		  + dT*rhoF/obj_.Ip().z()*( -fxR.z() + volIntegralRxU.z());
    	}
    }
    else if (mesh.nGeometricD() == 3)
    {
    	uRot_.x()
    	= 
    		uRot_.x() 
    	  + dT*rhoF/obj_.Ip().x()*( -fxR.x() + volIntegralRxU.x() );
    	uRot_.y() 
    	= 
    		uRot_.y() 
    	  + dT*rhoF/obj_.Ip().y()*( -fxR.y() + volIntegralRxU.y() );
    	uRot_.z() 
    	= 
    		uRot_.z() 
    	  + dT*rhoF/obj_.Ip().z()*( -fxR.z() + volIntegralRxU.z() );
    }
    else
    {
    	FatalErrorIn
    	(
    		"sixDoFMotion::updateMotion(vector& uTransl, "
			"vector& uRotate, const fvMesh& mesh," 
			"const vectorField& Fk, vector repulsiveForce," 
			"point center, pointField lPoints, scalar rhoF, "
			"scalar dV, vector g, vector volIntegralU, "
			"vector volIntegralRxU )"
    	)	<< "The mesh is neither 2D nor 3D " << endl
    		<< "Please check again"<<endl
    		<<exit(FatalError);
    }

	if (constraintDir_ != vector::zero)
	{
		uTrans_ = (uTrans_ & constraintDir_)*constraintDir_;
	}
}

Foam::vector& Foam::sixDoFMotion::translationalVelocity()
{
	return uTrans_;
}

Foam::vector& Foam::sixDoFMotion::rotationalVelocity()
{
	return uRot_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
