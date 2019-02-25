/*
Copyright (c) 2016, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/
#define NEW_CODE

#ifdef NEW_CODE
#define BIG_DATA								// Supports processing requiring more than 32-bit integers for indexing
												// Note: enabling BIG_DATA can generate .ply files using "longlong" for face indices instead of "int".
												// These are not standardly supported by .ply reading/writing applications.
#define NEW_THREADS								// Enabling this flag replaces the OpenMP implementation of parallelism with C++11's
#define FORCE_PARALLEL							// Forces parallel methods to pass in a thread pool
#endif // NEW_CODE


#undef ARRAY_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "PPolynomial.h"
#include "FEMTree.h"
#include "Ply.h"
#include "PointStreamData.h"
#include "Image.h"

cmdLineParameter< char* >
	In( "in" ) ,
	OutMesh( "mesh" ) ,
	OutGrid( "grid" );

cmdLineReadable
	PolygonMesh( "polygonMesh" ) ,
	NonManifold( "nonManifold" ) ,
	FlipOrientation( "flip" ) ,
	ASCII( "ascii" ) ,
	NonLinearFit( "nonLinearFit" ) ,
	PrimalGrid( "primalGrid" ) ,
	Verbose( "verbose" );

cmdLineParameter< int >
#ifdef NEW_THREADS
	ParallelType( "parallelType" , (int)ThreadPool::OPEN_MP ) ,
	ScheduleType( "scheduleType" , (int)ThreadPool::DefaultSchedule ) ,
	ThreadChunkSize( "tChunkSize" , (int)ThreadPool::DefaultChunkSize ) ,
	Threads( "threads" , (int)std::thread::hardware_concurrency() );
#else // !NEW_THREADS
	Threads( "threads" , omp_get_num_procs() );
#endif // NEW_THREADS

cmdLineParameter< float >
	IsoValue( "iso" , 0.f );

cmdLineReadable* params[] =
{
	&In , 
	&OutMesh , &NonManifold , &PolygonMesh , &FlipOrientation , &ASCII , &NonLinearFit , &IsoValue ,
	&OutGrid , &PrimalGrid ,
	&Threads ,
	&Verbose , 
#ifdef NEW_THREADS
	&ParallelType ,
	&ScheduleType ,
	&ThreadChunkSize ,
#endif // NEW_THREADS
	NULL
};


void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input tree>\n" , In.name );
	printf( "\t[--%s <ouput triangle mesh>]\n" , OutMesh.name );
	printf( "\t[--%s <ouput grid>]\n" , OutGrid.name );
#ifdef NEW_THREADS
	printf( "\t[--%s <num threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <parallel type>=%d]\n" , ParallelType.name , ParallelType.value );
	for( size_t i=0 ; i<ThreadPool::ParallelNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ParallelNames[i].c_str() );
	printf( "\t[--%s <schedue type>=%d]\n" , ScheduleType.name , ScheduleType.value );
	for( size_t i=0 ; i<ThreadPool::ScheduleNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ScheduleNames[i].c_str() );
	printf( "\t[--%s <thread chunk size>=%d]\n" , ThreadChunkSize.name , ThreadChunkSize.value );
#else // !NEW_THREADS
#ifdef _OPENMP
	printf( "\t[--%s <num threads>=%d]\n" , Threads.name , Threads.value );
#endif // _OPENMP
#endif // NEW_THREADS
	printf( "\t[--%s <iso-value for extraction>=%f]\n" , IsoValue.name , IsoValue.value );
	printf( "\t[--%s]\n" , NonManifold.name );
	printf( "\t[--%s]\n" , PolygonMesh.name );
	printf( "\t[--%s]\n" , NonLinearFit.name );
	printf( "\t[--%s]\n" , FlipOrientation.name );
	printf( "\t[--%s]\n" , PrimalGrid.name );
	printf( "\t[--%s]\n" , ASCII.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

template< typename Real , unsigned int Dim >
#ifdef NEW_THREADS
void WriteGrid( ThreadPool &tp , ConstPointer( Real ) values , int res , const char *fileName )
#else // !NEW_THREADS
void WriteGrid( ConstPointer( Real ) values , int res , const char *fileName )
#endif // NEW_THREADS
{
	int resolution = 1;
	for( int d=0 ; d<Dim ; d++ ) resolution *= res;

	char *ext = GetFileExtension( fileName );

	if( Dim==2 && ImageWriter::ValidExtension( ext ) )
	{
		Real avg = 0;
#ifdef NEW_THREADS
		std::vector< Real > avgs( tp.threadNum() , 0 );
		tp.parallel_for( 0 , resolution , [&]( const ThreadPool::ThreadNum &thread , size_t i ){ avgs[thread()] += values[i]; } );
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) avg += avgs[t];
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : avg )
		for( int i=0 ; i<resolution ; i++ ) avg += values[i];
#endif // NEW_THREADS
		avg /= (Real)resolution;

		Real std = 0;
#ifdef NEW_THREADS
		std::vector< Real > stds( tp.threadNum() , 0 );
		tp.parallel_for( 0 , resolution , [&]( const ThreadPool::ThreadNum &thread , size_t i ){ stds[thread()] += ( values[i] - avg ) * ( values[i] - avg ); } );
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) std += stds[t];
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : std )
		for( int i=0 ; i<resolution ; i++ ) std += ( values[i] - avg ) * ( values[i] - avg );
#endif // NEW_THREADS
		std = (Real)sqrt( std / resolution );

		if( Verbose.set ) printf( "Grid to image: [%.2f,%.2f] -> [0,255]\n" , avg - 2*std , avg + 2*std );

		unsigned char *pixels = new unsigned char[ resolution*3 ];
#ifdef NEW_THREADS
		tp.parallel_for( 0 , resolution , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
		for( int i=0 ; i<resolution ; i++ )
#endif // NEW_THREADS
		{
			Real v = (Real)std::min< Real >( (Real)1. , std::max< Real >( (Real)-1. , ( values[i] - avg ) / (2*std ) ) );
			v = (Real)( ( v + 1. ) / 2. * 256. );
			unsigned char color = (unsigned char )std::min< Real >( (Real)255. , std::max< Real >( (Real)0. , v ) );
			for( int c=0 ; c<3 ; c++ ) pixels[i*3+c ] = color;
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
		ImageWriter::Write( fileName , pixels , res , res , 3 );
		delete[] pixels;
	}
	else
	{

		FILE *fp = fopen( fileName , "wb" );
		if( !fp ) ERROR_OUT( "Failed to open grid file for writing: " , fileName );
		else
		{
			fwrite( &res , sizeof(int) , 1 , fp );
			if( typeid(Real)==typeid(float) ) fwrite( values , sizeof(float) , resolution , fp );
			else
			{
				float *fValues = new float[resolution];
				for( int i=0 ; i<resolution ; i++ ) fValues[i] = float( values[i] );
				fwrite( fValues , sizeof(float) , resolution , fp );
				delete[] fValues;
			}
			fclose( fp );
		}
	}
	delete[] ext;
}

template< unsigned int Dim , class Real , unsigned int FEMSig >
void _Execute( const FEMTree< Dim , Real >* tree , FILE* fp )
{
#ifdef NEW_THREADS
	ThreadPool tp( (ThreadPool::ParallelType)ParallelType.value , Threads.value );
#endif // NEW_THREADS
	static const unsigned int Degree = FEMSignature< FEMSig >::Degree;
	DenseNodeData< Real , IsotropicUIntPack< Dim , FEMSig > > coefficients;

	coefficients.read( fp );

	// Output the grid
	if( OutGrid.set )
	{
		int res = 0;
		double t = Time();
#ifdef NEW_THREADS
		Pointer( Real ) values = tree->template regularGridEvaluate< true >( tp , coefficients , res , -1 , PrimalGrid.set );
#else // !NEW_THREADS
		Pointer( Real ) values = tree->template regularGridEvaluate< true >( coefficients , res , -1 , PrimalGrid.set );
#endif // NEW_THREADS
		if( Verbose.set ) printf( "Got grid: %.2f(s)\n" , Time()-t );
#ifdef NEW_THREADS
		WriteGrid< Real , Dim >( tp , values , res , OutGrid.value );
#else // !NEW_THREADS
		WriteGrid< Real , Dim >( values , res , OutGrid.value );
#endif // NEW_THREADS
		DeletePointer( values );
	}

	// Output the mesh
	if( OutMesh.set )
	{
		double t = Time();
		typedef PlyVertex< Real , Dim > Vertex;
#ifdef NEW_CODE
		CoredFileMeshData< Vertex , node_index_type > mesh;
#else // !NEW_CODE
		CoredFileMeshData< Vertex > mesh;
#endif // NEW_CODE
		std::function< void ( Vertex& , Point< Real , Dim > , Real , Real ) > SetVertex = []( Vertex& v , Point< Real , Dim > p , Real , Real ){ v.point = p; };
#ifdef NEW_THREADS
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
		static const unsigned int DataSig = FEMDegreeAndBType< 0 , BOUNDARY_FREE >::Signature;
		IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< Real >( IsotropicUIntPack< Dim , FEMSig >() , UIntPack< 0 >() , UIntPack< FEMTrivialSignature >() , tp , *tree , ( typename FEMTree< Dim , Real >::template DensityEstimator< 0 >* )NULL , ( SparseNodeData< ProjectiveData< Real , Real > , IsotropicUIntPack< Dim , DataSig > > * )NULL , coefficients , IsoValue.value , mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , FlipOrientation.set );
#else // !__GNUC__ || __GNUC__ >=5
		IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< Real >( IsotropicUIntPack< Dim , FEMSig >() , UIntPack< 0 >() , UIntPack< FEMTrivialSignature >() , tp , *tree , ( typename FEMTree< Dim , Real >::template DensityEstimator< 0 >* )NULL , NULL , coefficients , IsoValue.value , mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , FlipOrientation.set );
#endif // __GNUC__ || __GNUC__ < 4
#else // !NEW_THREADS
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
		static const unsigned int DataSig = FEMDegreeAndBType< 0 , BOUNDARY_FREE >::Signature;
		IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< Real >( IsotropicUIntPack< Dim , FEMSig >() , UIntPack< 0 >() , UIntPack< FEMTrivialSignature >() , *tree , ( typename FEMTree< Dim , Real >::template DensityEstimator< 0 >* )NULL , ( SparseNodeData< ProjectiveData< Real , Real > , IsotropicUIntPack< Dim , DataSig > > * )NULL , coefficients , IsoValue.value , mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , FlipOrientation.set );
#else // !__GNUC__ || __GNUC__ >=5
		IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< Real >( IsotropicUIntPack< Dim , FEMSig >() , UIntPack< 0 >() , UIntPack< FEMTrivialSignature >() , *tree , ( typename FEMTree< Dim , Real >::template DensityEstimator< 0 >* )NULL , NULL , coefficients , IsoValue.value , mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , FlipOrientation.set );
#endif // __GNUC__ || __GNUC__ < 4
#endif // NEW_THREADS

		if( Verbose.set ) printf( "Got iso-surface: %.2f(s)\n" , Time()-t );
#ifdef NEW_CODE
		if( Verbose.set ) printf( "Vertices / Polygons: %llu / %llu\n" , (unsigned long long)( mesh.outOfCorePointCount()+mesh.inCorePoints.size() ) , (unsigned long long)mesh.polygonCount() );
#else // !NEW_CODE
		if( Verbose.set ) printf( "Vertices / Polygons: %d / %d\n" , (int)( mesh.outOfCorePointCount()+mesh.inCorePoints.size() ) , (int)mesh.polygonCount() );
#endif // NEW_CODE

		std::vector< std::string > comments;
#ifdef NEW_CODE
		if( !PlyWritePolygons< Vertex , node_index_type , Real , Dim >( OutMesh.value , &mesh , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , comments , XForm< Real , Dim+1 >::Identity() ) )
#else // !NEW_CODE
		if( !PlyWritePolygons< Vertex , Real , Dim >( OutMesh.value , &mesh , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , comments , XForm< Real , Dim+1 >::Identity() ) )
#endif // NEW_CODE
			ERROR_OUT( "Could not write mesh to: " , OutMesh.value );
	}
}


template< unsigned int Dim , class Real >
void Execute( FILE* fp , int degree , BoundaryType bType )
{
	FEMTree< Dim , Real > tree( fp , MEMORY_ALLOCATOR_BLOCK_SIZE );

#ifdef NEW_CODE
	if( Verbose.set ) printf( "Leaf Nodes / Active Nodes / Ghost Nodes: %llu / %llu / %llu\n" , (unsigned long long)tree.leaves() , (unsigned long long)tree.nodes() , (unsigned long long)tree.ghostNodes() );
#else // !NEW_CODE
	if( Verbose.set ) printf( "Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n" , (int)tree.leaves() , (int)tree.nodes() , (int)tree.ghostNodes() );
#endif // NEW_CODE

	switch( bType )
	{
	case BOUNDARY_FREE:
	{
		switch( degree )
		{
			case 1: _Execute< Dim , Real , FEMDegreeAndBType< 1 , BOUNDARY_FREE >::Signature >( &tree , fp ) ; break;
			case 2: _Execute< Dim , Real , FEMDegreeAndBType< 2 , BOUNDARY_FREE >::Signature >( &tree , fp ) ; break;
			case 3: _Execute< Dim , Real , FEMDegreeAndBType< 3 , BOUNDARY_FREE >::Signature >( &tree , fp ) ; break;
			case 4: _Execute< Dim , Real , FEMDegreeAndBType< 4 , BOUNDARY_FREE >::Signature >( &tree , fp ) ; break;
			default: ERROR_OUT( "Only B-Splines of degree 1 - 4 are supported" );
		}
	}
	break;
	case BOUNDARY_NEUMANN:
	{
		switch( degree )
		{
			case 1: _Execute< Dim , Real , FEMDegreeAndBType< 1 , BOUNDARY_NEUMANN >::Signature >( &tree , fp ) ; break;
			case 2: _Execute< Dim , Real , FEMDegreeAndBType< 2 , BOUNDARY_NEUMANN >::Signature >( &tree , fp ) ; break;
			case 3: _Execute< Dim , Real , FEMDegreeAndBType< 3 , BOUNDARY_NEUMANN >::Signature >( &tree , fp ) ; break;
			case 4: _Execute< Dim , Real , FEMDegreeAndBType< 4 , BOUNDARY_NEUMANN >::Signature >( &tree , fp ) ; break;
			default: ERROR_OUT( "Only B-Splines of degree 1 - 4 are supported" );
		}
	}
	break;
	case BOUNDARY_DIRICHLET:
	{
		switch( degree )
		{
			case 1: _Execute< Dim , Real , FEMDegreeAndBType< 1 , BOUNDARY_DIRICHLET >::Signature >( &tree , fp ) ; break;
			case 2: _Execute< Dim , Real , FEMDegreeAndBType< 2 , BOUNDARY_DIRICHLET >::Signature >( &tree , fp ) ; break;
			case 3: _Execute< Dim , Real , FEMDegreeAndBType< 3 , BOUNDARY_DIRICHLET >::Signature >( &tree , fp ) ; break;
			case 4: _Execute< Dim , Real , FEMDegreeAndBType< 4 , BOUNDARY_DIRICHLET >::Signature >( &tree , fp ) ; break;
			default: ERROR_OUT( "Only B-Splines of degree 1 - 4 are supported" );
		}
	}
	break;
	default: ERROR_OUT( "Not a valid boundary type: " , bType );
	}
}

int main( int argc , char* argv[] )
{
#ifdef ARRAY_DEBUG
	WARN( "Array debugging enabled" );
#endif // ARRAY_DEBUG
	cmdLineParse( argc-1 , &argv[1] , params );
#ifdef NEW_THREADS
	ThreadPool::DefaultChunkSize = ThreadChunkSize.value;
	ThreadPool::DefaultSchedule = (ThreadPool::ScheduleType)ScheduleType.value;
#else // !NEW_THREADS
	omp_set_num_threads( Threads.value > 1 ? Threads.value : 1 );
#endif // NEW_THREADS
	if( Verbose.set )
	{
		printf( "**************************************************\n" );
		printf( "**************************************************\n" );
		printf( "** Running Octree Visualization (Version %s) **\n" , VERSION );
		printf( "**************************************************\n" );
		printf( "**************************************************\n" );
	}

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	FILE* fp = fopen( In.value , "rb" );
	if( !fp ) ERROR_OUT( "Failed to open file for reading: " , In.value );
	FEMTreeRealType realType ; int degree ; BoundaryType bType;
#ifdef NEW_CODE
	unsigned int dimension;
#else // !NEW_CODE
	int dimension;
#endif // NEW_CODE
	ReadFEMTreeParameter( fp , realType , dimension );
	{
		unsigned int dim = dimension;
		unsigned int* sigs = ReadDenseNodeDataSignatures( fp , dim );
		if( dimension!=dim ) ERROR_OUT( "Octree and node data dimensions don't math: " , dimension , " != " , dim );
		for( unsigned int d=1 ; d<dim ; d++ ) if( sigs[0]!=sigs[d] ) ERROR_OUT( "Anisotropic signatures" );
		degree = FEMSignatureDegree( sigs[0] );
		bType = FEMSignatureBType( sigs[0] );
		delete[] sigs;
	}
	if( Verbose.set ) printf( "%d-dimension , %s-precision , degree-%d , %s-boundary\n" , dimension , FEMTreeRealNames[ realType ] , degree , BoundaryNames[ bType ] );

	switch( dimension )
	{
	case 2:
		switch( realType )
		{
			case FEM_TREE_REAL_FLOAT:  Execute< 2 , float  >( fp , degree , bType ) ; break;
			case FEM_TREE_REAL_DOUBLE: Execute< 2 , double >( fp , degree , bType ) ; break;
			default: ERROR_OUT( "Unrecognized real type: " , realType );
		}
		break;
	case 3:
		switch( realType )
		{
			case FEM_TREE_REAL_FLOAT:  Execute< 3 , float  >( fp , degree , bType ) ; break;
			case FEM_TREE_REAL_DOUBLE: Execute< 3 , double >( fp , degree , bType ) ; break;
			default: ERROR_OUT( "Unrecognized real type: " , realType );
		}
		break;
	default: ERROR_OUT( "Only dimensions 1-4 supported" );
	}

	fclose( fp );
	return EXIT_SUCCESS;
}
