/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBOStream class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBOStream.H"
#include "OFstream.H"
#include "masterUncollatedFileOperation.H"

namespace Foam
{
	defineTypeNameAndDebug(IBOStream, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::IBOStream::readDict()
{
	dictionary writeDict
    (
        IOdictionary
        (
            IOobject
            (
                "IBMDict",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).optionalSubDict("writeData")
    );
    writeInterval_ = writeDict.lookupOrDefault<scalar>("writeInterval", 1);
    writeLagrPoints_ = writeDict.lookupOrDefault<Switch>("writeLagrPoints", false);
    writeShadPoints_ = writeDict.lookupOrDefault<Switch>("writeShadPoints", false);
    writeLagrForces_ = writeDict.lookupOrDefault<Switch>("writeLagrForces", false);
    writeIBForce_ = writeDict.lookupOrDefault<Switch>("writeIBForce", false);
    writeObjVTU_ = writeDict.lookupOrDefault<Switch>("writeObjVTU", false);
    writeObjData_ = writeDict.lookupOrDefault<Switch>("writeObjData", false);
}

Foam::IBOStream::IBOStream
(
	dynamicFvMesh& mesh
)
:
	mesh_(mesh),
	writeInterval_(),
	writeLagrPoints_(NULL),
	writeLagrForces_(NULL),
	writeIBForce_(NULL),
	writeObjVTU_(NULL)
{
	readDict();
}

void Foam::IBOStream::writeLagrPoints
(
	const label objectID,
	const pointField& lPoints
)
{
	if (writeLagrPoints_)
	{
		if (    mesh_.time().timeIndex() == 0 
			 || mesh_.time().outputTime()
		   )
		{
			if (Pstream::master())
			{
				Info<< "IBM: Writing Lagrangian points of object " <<objectID+1<<endl;

				fileName outputDir(fileName::null); 

				if(Pstream::parRun())
				{
					outputDir = mesh_.time().path()/".."/"outputData"/"LagrPoints";
				}
				else
				{
					outputDir = mesh_.time().path()/"outputData"/"LagrPoints";
				}

				mkDir (outputDir);
				
				autoPtr<OFstream> outputFilePtr;
				
				word outputFile = "LagrPoints_Object-"+Foam::name(objectID+1)
									+"_Time-"+Foam::name(mesh_.time().timeIndex())+".csv";
				
				outputFilePtr.reset(new OFstream(outputDir/outputFile));
		    
			
				outputFilePtr() << "X, Y, Z" << endl;
		    
		    	forAll(lPoints, i)
		    	{	
					outputFilePtr() <<lPoints[i].x()<<", "<<lPoints[i].y()<<", "<<lPoints[i].z()<< endl;
				}
			}
		}
	}
}

void Foam::IBOStream::writeShadPoints
(
	const label objectID,
	const pointField& shd1Points,
	const pointField& shd2Points
)
{
	if (writeShadPoints_)
	{
		if (    mesh_.time().timeIndex() == 0 
			 || mesh_.time().outputTime()
		   )
		{
			if (Pstream::master())
		    {
				Info<< "IBM: Writing Shadow points of object " <<objectID+1<< endl;

				fileName outputDir(fileName::null); 

				if(Pstream::parRun())
				{
					outputDir = mesh_.time().path()/".."/"outputData"/"ShadPoints";
				}
				else
				{
					outputDir = mesh_.time().path()/"outputData"/"ShadPoints";
				}

				mkDir (outputDir);
				
				autoPtr<OFstream> outputFilePtr1;
				autoPtr<OFstream> outputFilePtr2;
				
				word outputFile1 = "Shad1Points_Object-"+Foam::name(objectID+1)
									+"_Time-"+Foam::name(mesh_.time().timeIndex())+".csv";
				word outputFile2 = "Shad2Points_Object-"+Foam::name(objectID+1)
									+"_Time-"+Foam::name(mesh_.time().timeIndex())+".csv";
				
				outputFilePtr1.reset(new OFstream(outputDir/outputFile1));
				outputFilePtr2.reset(new OFstream(outputDir/outputFile2));
				
				outputFilePtr1() << "X, Y, Z" << endl;
				outputFilePtr2() << "X, Y, Z" << endl;
				
				forAll(shd1Points, i)
				{
					outputFilePtr1() <<shd1Points[i].x()<<", "<<shd1Points[i].y()<<", "<<shd1Points[i].z()<< endl;
					outputFilePtr2() <<shd2Points[i].x()<<", "<<shd2Points[i].y()<<", "<<shd2Points[i].z()<< endl;
				}
			}
		}
	}
}

void Foam::IBOStream::writeLagrForces
(
	const label objectID,
	const pointField& lPoints,
	const vectorField& lForces
)
{
	if (writeLagrForces_)
	{
		if (    mesh_.time().timeIndex() == 0 
			 || mesh_.time().outputTime()
		   )
		{  
			if (Pstream::master())
			{
				Info<< "IBM: Writing Lagrangian Force of object "<<objectID+1<< endl;
				
				fileName outputDir(fileName::null); 
				
				if(Pstream::parRun())
				{
					outputDir = mesh_.time().path()/".."/"outputData"/"LagrForces";
				}
				else
				{
					outputDir = mesh_.time().path()/"outputData"/"LagrForces";
				}
				
				mkDir (outputDir);
				
				autoPtr<OFstream> outputFilePtr;
				
				word outputFile = "LagrForces_Object-"+Foam::name(objectID+1)
									+"_Time-"+Foam::name(mesh_.time().timeIndex())+".csv";
				
				outputFilePtr.reset(new OFstream(outputDir/outputFile));

				outputFilePtr() << "X, Y, Z, Fx, Fy, Fz" << endl;
				
				forAll(lPoints, i)
				{
					outputFilePtr() <<lPoints[i].x()<<", "<<lPoints[i].y()<<", "<<lPoints[i].z() <<", "
									<<lForces[i].x()<<", "<<lForces[i].y()<<", "<<lForces[i].z()<< endl;
				}
			}
		}
	}
}
void Foam::IBOStream::writeIBForces
(
	const label objectID,
	const vectorField& lForces,
	const scalar dVolume
)
{
	if(writeIBForce_)
	{
	    if(Pstream::master())
	    {
			//- Calculate force acting upon object's surface

			vector ibForce(vector::zero);
			
			forAll(lForces, i)
			{
			    ibForce -= lForces[i]*dVolume;
			}

			Info<< "IBM: Writing ibForces"<< nl
				<< "    Object-" << objectID+1 <<" : "<<ibForce<<endl;

			//- Write
			//- Create dir
			fileName outputDir(fileName::null); 
			
			if(Pstream::parRun())
			{
				outputDir = mesh_.time().path()/".."/"outputData"/"Forces";
			}
			else
			{
				outputDir = mesh_.time().path()/"outputData"/"Forces";
			}
		
			mkDir (outputDir);
		
			//-Create file
			const fileName outputFile(outputDir/"Forces_Object-"
										+Foam::name(objectID+1)+".csv");
			IOstream::streamFormat format=IOstream::ASCII;
			IOstream::versionNumber version=IOstream::currentVersion;
			IOstream::compressionType compression=IOstream::UNCOMPRESSED;

			const bool append(true);

			OFstream os(outputFile, format, version, compression, append);

			//- Write file
			scalar magForce = mag(ibForce);
			if (mesh_.time().timeIndex() == 0 && Pstream::master())
			{
				os 	<< "Time, Fx, Fy, Fz, magF" << endl;
				os  << mesh_.time().timeName() << ", "
					<< ibForce.x()             << ", "
					<< ibForce.y()             << ", "
					<< ibForce.z()             << ", "
					<< magForce				   << endl; 
			}
			if (mesh_.time().timeIndex() % writeInterval_ == 0 && Pstream::master())
			{
				os  << mesh_.time().timeName() << ", "
					<< ibForce.x()             << ", "
					<< ibForce.y()             << ", "
					<< ibForce.z()             << ", "
					<< magForce				   << endl; 
			}
		}
	}
}
void Foam::IBOStream::writeObjVTU
(
	label objectID,
    List<point> Nodes, 
    point pntCenter,
    label numberFaces_,
    List<label> nPntsOfFace_,
    List<List<label> > PntsOfFace_
)
{
	label nPoints = Nodes.size();

	if (writeObjVTU_)
	{
		if (    mesh_.time().timeIndex() == 0 
			 || mesh_.time().outputTime()
		)
		{
		    //- Start writing
			if (Pstream::master())
			{
				//- Create output file
				fileName outputDir(fileName::null);

				if(Pstream::parRun())
				{
				    outputDir = mesh_.time().path()/".."/"outputData"/"ObjectVTU";
				}
				else
				    outputDir = mesh_.time().path()/"outputData"/"ObjectVTU";

				mkDir (outputDir);
				
				autoPtr<OFstream> outputFilePtr;

				word outputFile = "Object-"+Foam::name(objectID+1)+"_Time-"
									+Foam::name(mesh_.time().timeIndex())+".vtu";

				outputFilePtr.reset(new OFstream(outputDir/outputFile));
		    
				outputFilePtr() << "<?xml version=\"1.0\"?>" << "\n";
				outputFilePtr() << "<VTKFile type=\"UnstructuredGrid\">" << "\n";
				outputFilePtr() << "<UnstructuredGrid>" << "\n";

				if(mesh_.nGeometricD()==2)
				{
					outputFilePtr()	<< "<Piece NumberOfPoints=\"" << Nodes.size()+1
									<< "\" NumberOfCells=\"" 	  << numberFaces_ << "\">" << "\n";
				}
				if(mesh_.nGeometricD()==3)
				{
					outputFilePtr() << "<Piece NumberOfPoints=\"" << Nodes.size()
									<< "\" NumberOfCells=\""      << numberFaces_ << "\">" << "\n";
				}
				
				outputFilePtr() << "<Points>" << "\n";
				outputFilePtr() << "<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >" << "\n";

				forAll(Nodes, i)
				{
					outputFilePtr() << Nodes[i].x() << " "
									<< Nodes[i].y() << " "
									<< Nodes[i].z() << "\n";
				}

				if(mesh_.nGeometricD()==2)
				{
					outputFilePtr() << pntCenter.x() << " "
									<< pntCenter.y() << " "
									<< pntCenter.z() << "\n";
				}

				outputFilePtr() << "</DataArray>" << "\n";
				outputFilePtr() << "</Points>"    << "\n";
				outputFilePtr() << "<Cells>"      << "\n";
				outputFilePtr() << "<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\" >" << "\n";

				for (int i=0;i<numberFaces_;i++)
				{
					for (int j=0;j<nPntsOfFace_[i];j++)
					{
						outputFilePtr() << PntsOfFace_[i][j] << "\t";
					}

					outputFilePtr() << "\n";
				}

				outputFilePtr() << "</DataArray>" << "\n";
				outputFilePtr() << "<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\" >" << "\n";

				scalar offset=0;
				
				for (int i=0;i<numberFaces_;i++)
				{
					offset+=nPntsOfFace_[i];/// Triangular face
					outputFilePtr() << offset << " ";
					if((i+1)%6==0)
					{
						outputFilePtr() << "\n";
					} 
				}
				outputFilePtr() << "\n</DataArray>" << "\n";
				outputFilePtr() << "<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" >" << "\n";

				for (int i=0;i<numberFaces_;i++)
				{
					if(nPntsOfFace_[i]==3)
					{
						outputFilePtr() << "5 "; /// Triangular face
					}
					if(nPntsOfFace_[i]==4)
					{
						outputFilePtr() << "9 "; /// Quard
					}
					if(nPntsOfFace_[i] > 4)
					{
						outputFilePtr() << "6 "; /// Triangular strip for top and bottom faces of equal sphere
					}
					if((i+1)%6==0)
					{
						outputFilePtr() << "\n";
					} 
				}

				outputFilePtr() << "\n";
				outputFilePtr() << "</DataArray>" << "\n";
				outputFilePtr() << "</Cells>"     << "\n";
				/// Vectors at nodes:
				outputFilePtr() << "<CellData Vectors=\"" << "Fk" << "\">" << "\n";
				outputFilePtr() << "<DataArray Name=\""   << "Fk" 
					<< "\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\"3\">"<<"\n";
				if(mesh_.nGeometricD()==2)
				{
					for (int i=0;i<numberFaces_;i++)
					{
						if(i<numberFaces_/2.0)
							outputFilePtr()<<0<<" "<<0<<" "<<0<<" ";
						else
							outputFilePtr()<<1<<" "<<1<<" "<<1<<" ";
						// outputFilePtr()<<(Fk[PntsOfFace_[i][0]].x()+Fk[PntsOfFace_[i][1]].x())/2.0<<" "
						// 				  <<(Fk[PntsOfFace_[i][0]].y()+Fk[PntsOfFace_[i][1]].y())/2.0<<" "
						//  			  <<0<<" ";
						if((i+1)%2==0) outputFilePtr()<<"\n";
					}
				}
				if(mesh_.nGeometricD()==3)
				{
					if(numberFaces_ == nPoints)
					{
						for (int i=0;i<numberFaces_;i++)
						{
							if(i<numberFaces_/2.0)
								outputFilePtr()<<0<<" "<<0<<" "<<0<<" ";
							else
								outputFilePtr()<<1<<" "<<1<<" "<<1<<" ";
							// outputFilePtr()<<Fk[i].x()<<" "
							//     <<Fk[i].y()<<" "
							//     <<Fk[i].z()<<" ";
							if((i+1)%2==0) outputFilePtr()<<"\n";
						}
					}
					else
					{
						for (int i=0;i<numberFaces_;i++)
						{
							if(i<numberFaces_/2.0)
								outputFilePtr()<<0<<" "<<0<<" "<<0<<" ";
							else
								outputFilePtr()<<1<<" "<<1<<" "<<1<<" ";
							// outputFilePtr()<<(Fk[PntsOfFace_[i][0]].x()+Fk[PntsOfFace_[i][1]].x()+Fk[PntsOfFace_[i][2]].x())/3.0<<" "
							//     <<(Fk[PntsOfFace_[i][0]].y()+Fk[PntsOfFace_[i][1]].y()+Fk[PntsOfFace_[i][2]].y())/3.0<<" "
							//     <<(Fk[PntsOfFace_[i][0]].z()+Fk[PntsOfFace_[i][1]].z()+Fk[PntsOfFace_[i][2]].z())/3.0<<" ";
							if((i+1)%2==0) outputFilePtr()<<"\n";
						}
					}
				}
				outputFilePtr()<<"\n</DataArray>"<<"\n";
				outputFilePtr()<<"</CellData>"<<"\n";
				outputFilePtr()<<"</Piece>"<<"\n";
				outputFilePtr()<<"</UnstructuredGrid>"<<"\n";
				outputFilePtr()<<"</VTKFile>"<<"\n";
			}
		}
	}
}

void Foam::IBOStream::writeObjData
(
	const label objectID,
	const point CG,
	const vector uTransl,
	const vector uRotate
)
{
	if(writeObjData_)
	{
				    
		if (Pstream::master())
		{
			Info<< "IBM: Writing data of object "<<objectID+1<< endl;

			//- Create dir
			fileName outputDir(fileName::null); 
				
			if(Pstream::parRun())
			{
				outputDir = mesh_.time().path()/".."/"outputData"/"ObjectsData";
			}
			else
			{
				outputDir = mesh_.time().path()/"outputData"/"ObjectsData";
			}
			
			mkDir (outputDir);
			
			//-Create file
			const fileName outputFile(outputDir/"Info_Object-"
										+Foam::name(objectID+1)+".csv");

			IOstream::streamFormat format=IOstream::ASCII;
			IOstream::versionNumber version=IOstream::currentVersion;
			IOstream::compressionType compression=IOstream::UNCOMPRESSED;

			const bool append(true);

			OFstream os(outputFile, format, version, compression, append);

			//- Write file
			if (mesh_.time().timeIndex() == 0)
			{
				os 	<< "Time, "
					<< "CG_x, CG_y, CG_z," 
					<< "Ux, Uy, Uz, "
					<< "Omega_x, Omega_y, Omega_z" << endl;
			}
			if (mesh_.time().timeIndex() % writeInterval_ == 0)
			{
				os  << mesh_.time().timeName()	<< ", "
					<< CG.x()             		<< ", "
					<< CG.y()             		<< ", "
					<< CG.z()             		<< ", "
					<< uTransl.x()             	<< ", "
					<< uTransl.y()             	<< ", "
					<< uTransl.z()             	<< ", "
					<< uRotate.x()             	<< ", "
					<< uRotate.y()             	<< ", "
					<< uRotate.z()             	<< ", "<< endl; 
			}
		}
	}
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
