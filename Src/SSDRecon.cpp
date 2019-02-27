/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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
#include "PreProcessor.h"
#endif // NEW_CODE

#undef USE_DOUBLE								// If enabled, double-precesion is used
#define DATA_DEGREE 0							// The order of the B-Spline used to splat in data for color interpolation
												// This can be changed to zero if more interpolatory performance is desired.
#define WEIGHT_DEGREE 2							// The order of the B-Spline used to splat in the weights for density estimation
#define NORMAL_DEGREE 2							// The order of the B-Spline used to splat int the normals for constructing the Laplacian constraints
#define DEFAULT_FEM_DEGREE 2					// The default finite-element degree
#define DEFAULT_FEM_BOUNDARY BOUNDARY_NEUMANN	// The default finite-element boundary type
#define DIMENSION 3								// The dimension of the system

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

MessageWriter messageWriter;

double BaseSSDWeights[] = { 5e+1f , 5e-4f , 1e-5f }; 

cmdLineParameter< char* >
	In( "in" ) ,
	Out( "out" ) ,
	TempDir( "tempDir" ) ,
	Grid( "grid" ) ,	
	Tree( "tree" ) ,
	Transform( "xForm" );

cmdLineReadable
	Performance( "performance" ) ,
	ShowResidual( "showResidual" ) ,
	NoComments( "noComments" ) ,
	PolygonMesh( "polygonMesh" ) ,
	NonManifold( "nonManifold" ) ,
	ASCII( "ascii" ) ,
	Density( "density" ) ,
	NonLinearFit( "nonLinearFit" ) ,
	PrimalGrid( "primalGrid" ) ,
	ExactInterpolation( "exact" ) ,
	Normals( "normals" ) ,
	Colors( "colors" ) ,
	InCore( "inCore" ) ,
	Verbose( "verbose" );

cmdLineParameter< int >
#ifndef FAST_COMPILE
	Degree( "degree" , DEFAULT_FEM_DEGREE ) ,
#endif // !FAST_COMPILE
	Depth( "depth" , 8 ) ,
	KernelDepth( "kernelDepth" ) ,
	Iters( "iters" , 8 ) ,
	FullDepth( "fullDepth" , 5 ) ,
	BaseDepth( "baseDepth" , 5 ) ,
	BaseVCycles( "baseVCycles" , 4 ) ,
#ifndef FAST_COMPILE
	BType( "bType" , DEFAULT_FEM_BOUNDARY+1 ) ,
#endif // !FAST_COMPILE
	MaxMemoryGB( "maxMemory" , 0 ) ,
#ifdef NEW_THREADS
	ParallelType( "parallel" , (int)ThreadPool::OPEN_MP ) ,
	ScheduleType( "schedule" , (int)ThreadPool::DefaultSchedule ) ,
	ThreadChunkSize( "tChunkSize" , (int)ThreadPool::DefaultChunkSize ) ,
	Threads( "threads" , (int)std::thread::hardware_concurrency() );
#else // !NEW_THREADS
	Threads( "threads" , omp_get_num_procs() );
#endif // NEW_THREADS

cmdLineParameter< float >
	DataX( "data" , 32.f ) ,
	SamplesPerNode( "samplesPerNode" , 1.5f ) ,
	Scale( "scale" , 1.1f ) ,
	Width( "width" , 0.f ) ,
	Confidence( "confidence" , 0.f ) ,
	ConfidenceBias( "confidenceBias" , 0.f ) ,
	CGSolverAccuracy( "cgAccuracy" , 1e-3f ) ,
	ValueWeight   (    "valueWeight" , 1.f ) ,
	GradientWeight( "gradientWeight" , 1.f ) ,
	BiLapWeight   (    "biLapWeight" , 1.f );


cmdLineReadable* params[] =
{
#ifndef FAST_COMPILE
	&Degree , &BType ,
#endif // !FAST_COMPILE
	&In , &Depth , &Out , &Transform ,
	&Width ,
	&Scale , &Verbose , &CGSolverAccuracy , &NoComments ,
	&KernelDepth , &SamplesPerNode , &Confidence , &NonManifold , &PolygonMesh , &ASCII , &ShowResidual ,
	&ConfidenceBias ,
	&ValueWeight , &GradientWeight , &BiLapWeight ,
	&Grid , &Threads ,
	&Tree ,
	&Density ,
	&FullDepth ,
	&BaseDepth , &BaseVCycles ,
	&Iters ,
	&DataX ,
	&Colors ,
	&Normals ,
	&NonLinearFit ,
	&PrimalGrid ,
	&TempDir ,
	&ExactInterpolation ,
	&Performance ,
	&MaxMemoryGB ,
	&InCore ,
#ifdef NEW_THREADS
	&ParallelType ,
	&ScheduleType ,
	&ThreadChunkSize ,
#endif // NEW_THREADS
	NULL
};

void ShowUsage(char* ex)
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input points>\n" , In.name );
	printf( "\t[--%s <ouput triangle mesh>]\n" , Out.name );
	printf( "\t[--%s <ouput grid>]\n" , Grid.name );
	printf( "\t[--%s <ouput fem tree>]\n" , Tree.name );
#ifndef FAST_COMPILE
	printf( "\t[--%s <b-spline degree>=%d]\n" , Degree.name , Degree.value );
	printf( "\t[--%s <boundary type>=%d]\n" , BType.name , BType.value );
	for( int i=0 ; i<BOUNDARY_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i+1 , BoundaryNames[i] );
#endif // !FAST_COMPILE
	printf( "\t[--%s <maximum reconstruction depth>=%d]\n" , Depth.name , Depth.value );
	printf( "\t[--%s <grid width>]\n" , Width.name );
	printf( "\t[--%s <full depth>=%d]\n" , FullDepth.name , FullDepth.value );
	printf( "\t[--%s <coarse MG solver depth>=%d]\n" , BaseDepth.name , BaseDepth.value );
	printf( "\t[--%s <coarse MG solver v-cycles>=%d]\n" , BaseVCycles.name , BaseVCycles.value );
	printf( "\t[--%s <scale factor>=%f]\n" , Scale.name , Scale.value );
	printf( "\t[--%s <minimum number of samples per node>=%f]\n" , SamplesPerNode.name, SamplesPerNode.value );
	printf( "\t[--%s <zero-crossing weight>=%.3e]\n" , ValueWeight.name , ValueWeight.value );
	printf( "\t[--%s <gradient weight>=%.3e]\n" , GradientWeight.name , GradientWeight.value );
	printf( "\t[--%s <bi-laplacian weight>=%.3e]\n" , BiLapWeight.name , BiLapWeight.value );
	printf( "\t[--%s <iterations>=%d]\n" , Iters.name , Iters.value );
	printf( "\t[--%s]\n" , ExactInterpolation.name );
	printf( "\t[--%s <pull factor>=%f]\n" , DataX.name , DataX.value );
	printf( "\t[--%s]\n" , Colors.name );
	printf( "\t[--%s]\n" , Normals.name );
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
	printf( "\t[--%s <normal confidence exponent>=%f]\n" , Confidence.name , Confidence.value );
	printf( "\t[--%s <normal confidence bias exponent>=%f]\n" , ConfidenceBias.name , ConfidenceBias.value );
	printf( "\t[--%s]\n" , NonManifold.name );
	printf( "\t[--%s]\n" , PolygonMesh.name );
	printf( "\t[--%s <cg solver accuracy>=%g]\n" , CGSolverAccuracy.name , CGSolverAccuracy.value );
	printf( "\t[--%s <maximum memory (in GB)>=%d]\n" , MaxMemoryGB.name , MaxMemoryGB.value );
	printf( "\t[--%s]\n" , Performance.name );
	printf( "\t[--%s]\n" , Density.name );
	printf( "\t[--%s]\n" , NonLinearFit.name );
	printf( "\t[--%s]\n" , PrimalGrid.name );
	printf( "\t[--%s]\n" , ASCII.name );
	printf( "\t[--%s]\n" , NoComments.name );
	printf( "\t[--%s]\n" , TempDir.name );
	printf( "\t[--%s]\n" , InCore.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

double Weight( double v , double start , double end )
{
	v = ( v - start ) / ( end - start );
	if     ( v<0 ) return 1.;
	else if( v>1 ) return 0.;
	else
	{
		// P(x) = a x^3 + b x^2 + c x + d
		//		P (0) = 1 , P (1) = 0 , P'(0) = 0 , P'(1) = 0
		// =>	d = 1 , a + b + c + d = 0 , c = 0 , 3a + 2b + c = 0
		// =>	c = 0 , d = 1 , a + b = -1 , 3a + 2b = 0
		// =>	a = 2 , b = -3 , c = 0 , d = 1
		// =>	P(x) = 2 x^3 - 3 x^2 + 1
		return 2. * v * v * v - 3. * v * v + 1.;
	}
}

template< unsigned int Dim , class Real >
struct FEMTreeProfiler
{
	FEMTree< Dim , Real >& tree;
	double t;

	FEMTreeProfiler( FEMTree< Dim , Real >& t ) : tree(t) { ; }
	void start( void ){ t = Time() , FEMTree< Dim , Real >::ResetLocalMemoryUsage(); }
	void print( const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) printf( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
		else         printf(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
	}
	void dumpOutput( const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) messageWriter( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
		else         messageWriter(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
	}
	void dumpOutput2( std::vector< std::string >& comments , const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) messageWriter( comments , "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
		else         messageWriter( comments ,    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
	}
};

template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetBoundingBoxXForm( Point< Real , Dim > min , Point< Real , Dim > max , Real scaleFactor )
{
	Point< Real , Dim > center = ( max + min ) / 2;
	Real scale = max[0] - min[0];
	for( int d=1 ; d<Dim ; d++ ) scale = std::max< Real >( scale , max[d]-min[d] );
	scale *= scaleFactor;
	for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
	XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
	for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
	return sXForm * tXForm;
}
template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetBoundingBoxXForm( Point< Real , Dim > min , Point< Real , Dim > max , Real width , Real scaleFactor , int& depth )
{
	// Get the target resolution (along the largest dimension)
	Real resolution = ( max[0]-min[0] ) / width;
	for( int d=1 ; d<Dim ; d++ ) resolution = std::max< Real >( resolution , ( max[d]-min[d] ) / width );
	resolution *= scaleFactor;
	depth = 0;
	while( (1<<depth)<resolution ) depth++;

	Point< Real , Dim > center = ( max + min ) / 2;
	Real scale = (1<<depth) * width;

	for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
	XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
	for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
	return sXForm * tXForm;
}

template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetPointXForm( InputPointStream< Real , Dim >& stream , Real width , Real scaleFactor , int& depth )
{
	Point< Real , Dim > min , max;
	stream.boundingBox( min , max );
	return GetBoundingBoxXForm( min , max , width , scaleFactor , depth );
}
template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetPointXForm( InputPointStream< Real , Dim >& stream , Real scaleFactor )
{
	Point< Real , Dim > min , max;
	stream.boundingBox( min , max );
	return GetBoundingBoxXForm( min , max , scaleFactor );
}

template< unsigned int Dim , typename Real , typename TotalPointSampleData >
struct ConstraintDual
{
	Real target , vWeight , gWeight;
	ConstraintDual( Real t , Real v , Real g ) : target(t) , vWeight(v) , gWeight(g) { }
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( const Point< Real , Dim >& p , const TotalPointSampleData& data ) const 
	{
		Point< Real , Dim > n = std::get<0>( data.data ).data;
		CumulativeDerivativeValues< Real , Dim , 1 > cdv;
		cdv[0] = target*vWeight;
		for( int d=0 ; d<Dim ; d++ ) cdv[1+d] = -n[d]*gWeight;
		return cdv;
	}
};
template< unsigned int Dim , typename Real , typename TotalPointSampleData >
struct SystemDual
{
	CumulativeDerivativeValues< Real , Dim , 1 > weight;
	SystemDual( Real v , Real g )
	{
		weight[0] = v;
		for( int d=0 ; d<Dim ; d++ ) weight[d+1] = g;
	}
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const TotalPointSampleData& data , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
	{
		return dValues * weight;
	}
	CumulativeDerivativeValues< double , Dim , 1 > operator()( Point< Real , Dim > p , const TotalPointSampleData& data , const CumulativeDerivativeValues< double , Dim , 1 >& dValues ) const
	{
		return dValues * weight;
	};
};
template< unsigned int Dim , class TotalPointSampleData >
struct SystemDual< Dim , double , TotalPointSampleData >
{
	typedef double Real;
	CumulativeDerivativeValues< Real , Dim , 1 > weight;
	SystemDual( Real v , Real g ) : weight( v , g , g , g ) { }
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const TotalPointSampleData& data , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
	{
		return dValues * weight;
	}
};

template< typename Vertex , typename Real , unsigned int ... FEMSigs , typename ... SampleData >
#ifdef NEW_THREADS
void ExtractMesh( UIntPack< FEMSigs ... > , std::tuple< SampleData ... > , ThreadPool &tp , FEMTree< sizeof ... ( FEMSigs ) , Real >& tree , const DenseNodeData< Real , UIntPack< FEMSigs ... > >& solution , Real isoValue , const std::vector< typename FEMTree< sizeof ... ( FEMSigs ) , Real >::PointSample >* samples , std::vector< MultiPointStreamData< Real , PointStreamNormal< Real , DIMENSION > , MultiPointStreamData< Real , SampleData ... > > >* sampleData , const typename FEMTree< sizeof ... ( FEMSigs ) , Real >::template DensityEstimator< WEIGHT_DEGREE >* density , std::function< void ( Vertex& , Point< Real , DIMENSION > , Real , MultiPointStreamData< Real , PointStreamNormal< Real , DIMENSION > , MultiPointStreamData< Real , SampleData ... > > ) > SetVertex , std::vector< std::string > &comments , XForm< Real , sizeof...(FEMSigs)+1 > iXForm )
#else // !NEW_THREADS
void ExtractMesh( UIntPack< FEMSigs ... > , std::tuple< SampleData ... > , FEMTree< sizeof ... ( FEMSigs ) , Real >& tree , const DenseNodeData< Real , UIntPack< FEMSigs ... > >& solution , Real isoValue , const std::vector< typename FEMTree< sizeof ... ( FEMSigs ) , Real >::PointSample >* samples , std::vector< MultiPointStreamData< Real , PointStreamNormal< Real , DIMENSION > , MultiPointStreamData< Real , SampleData ... > > >* sampleData , const typename FEMTree< sizeof ... ( FEMSigs ) , Real >::template DensityEstimator< WEIGHT_DEGREE >* density , std::function< void ( Vertex& , Point< Real , DIMENSION > , Real , MultiPointStreamData< Real , PointStreamNormal< Real , DIMENSION > , MultiPointStreamData< Real , SampleData ... > > ) > SetVertex , std::vector< std::string > &comments , XForm< Real , sizeof...(FEMSigs)+1 > iXForm )
#endif // NEW_THREADS
{
	static const int Dim = sizeof ... ( FEMSigs );
	typedef UIntPack< FEMSigs ... > Sigs;
	typedef PointStreamNormal< Real , Dim > NormalPointSampleData;
	typedef MultiPointStreamData< Real , SampleData ... > AdditionalPointSampleData;
	typedef MultiPointStreamData< Real , NormalPointSampleData , AdditionalPointSampleData > TotalPointSampleData;
	static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
	typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;

	FEMTreeProfiler< Dim , Real > profiler( tree );

	char tempHeader[1024];
	{
		char tempPath[1024];
		tempPath[0] = 0;
		if( TempDir.set ) strcpy( tempPath , TempDir.value );
		else SetTempDirectory( tempPath , sizeof(tempPath) );
		if( strlen(tempPath)==0 ) sprintf( tempPath , ".%c" , FileSeparator );
		if( tempPath[ strlen( tempPath )-1 ]==FileSeparator ) sprintf( tempHeader , "%sPR_" , tempPath );
		else                                                  sprintf( tempHeader , "%s%cPR_" , tempPath , FileSeparator );
	}

#ifdef NEW_CODE
	CoredMeshData< Vertex , node_index_type > *mesh;
	if( InCore.set ) mesh = new CoredVectorMeshData< Vertex , node_index_type >();
	else             mesh = new CoredFileMeshData< Vertex , node_index_type >( tempHeader );
#else // !NEW_CODE
	CoredMeshData< Vertex > *mesh;
	if( InCore.set ) mesh = new CoredVectorMeshData< Vertex >();
	else             mesh = new CoredFileMeshData< Vertex >( tempHeader );
#endif // NEW_CODE
	profiler.start();
	typename IsoSurfaceExtractor< Dim , Real , Vertex >::IsoStats isoStats;
	if( sampleData )
	{
		SparseNodeData< ProjectiveData< TotalPointSampleData , Real > , IsotropicUIntPack< Dim , DataSig > > _sampleData = tree.template setDataField< DataSig , false >( *samples , *sampleData , (DensityEstimator*)NULL );
#ifdef NEW_CODE
		for( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* n = tree.tree().nextNode() ; n ; n=tree.tree().nextNode( n ) )
#else // !NEW_CODE
		for( const RegularTreeNode< Dim , FEMTreeNodeData >* n = tree.tree().nextNode() ; n ; n=tree.tree().nextNode( n ) )
#endif // NEW_CODE
		{
			ProjectiveData< TotalPointSampleData , Real >* clr = _sampleData( n );
			if( clr ) (*clr) *= (Real)pow( DataX.value , tree.depth( n ) );
		}
#ifdef NEW_THREADS
		isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tp , tree , density , &_sampleData , solution , isoValue , *mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#else // !NEW_THREADS
		isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , &_sampleData , solution , isoValue , *mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#endif // NEW_THREADS
	}
#ifdef NEW_THREADS
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
	else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tp , tree , density , (SparseNodeData< ProjectiveData< TotalPointSampleData , Real > , IsotropicUIntPack< Dim , DataSig > > *)NULL , solution , isoValue , *mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#else // !__GNUC__ || __GNUC__ >=5
	else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tp , tree , density , NULL , solution , isoValue , *mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#endif // __GNUC__ || __GNUC__ < 4
#else // !NEW_THREADS
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
	else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , (SparseNodeData< ProjectiveData< TotalPointSampleData , Real > , IsotropicUIntPack< Dim , DataSig > > *)NULL , solution , isoValue , *mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#else // !__GNUC__ || __GNUC__ >=5
	else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , NULL , solution , isoValue , *mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#endif // __GNUC__ || __GNUC__ < 4
#endif // NEW_THREADS
#ifdef NEW_CODE
	messageWriter( "Vertices / Polygons: %llu / %llu\n" , (unsigned long long)( mesh->outOfCorePointCount()+mesh->inCorePoints.size() ) , (unsigned long long)mesh->polygonCount() );
#else // !NEW_CODE
	messageWriter( "Vertices / Polygons: %d / %d\n" , mesh->outOfCorePointCount()+mesh->inCorePoints.size() , mesh->polygonCount() );
#endif // NEW_CODE
	std::string isoStatsString = isoStats.toString() + std::string( "\n" );
	messageWriter( isoStatsString.c_str() );
	if( PolygonMesh.set ) profiler.dumpOutput2( comments , "#         Got polygons:" );
	else                  profiler.dumpOutput2( comments , "#        Got triangles:" );

	std::vector< std::string > noComments;
#ifdef NEW_CODE
	if( !PlyWritePolygons< Vertex , node_index_type , Real , Dim >( Out.value , mesh , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , NoComments.set ? noComments : comments , iXForm ) )
#else // !NEW_CODE
	if( !PlyWritePolygons< Vertex , Real , Dim >( Out.value , mesh , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , NoComments.set ? noComments : comments , iXForm ) )
#endif // NEW_CODE
		ERROR_OUT( "Could not write mesh to: " , Out.value );
	delete mesh;
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
		tp.parallel_for( 0 , resolution , [&]( unsigned int thread , size_t i ){ avgs[thread] += values[i]; } );
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) avg += avgs[t];
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : avg )
		for( int i=0 ; i<resolution ; i++ ) avg += values[i];
#endif // NEW_THREADS
		avg /= (Real)resolution;

		Real std = 0;
#ifdef NEW_THREADS
		std::vector< Real > stds( tp.threadNum() , 0 );
		tp.parallel_for( 0 , resolution , [&]( unsigned int thread , size_t i ){ stds[thread] += ( values[i] - avg ) * ( values[i] - avg ); } );
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) std += stds[t];
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : std )
		for( int i=0 ; i<resolution ; i++ ) std += ( values[i] - avg ) * ( values[i] - avg );
#endif // NEW_THREADS
		std = (Real)sqrt( std / resolution );

		if( Verbose.set ) printf( "Grid to image: [%.2f,%.2f] -> [0,255]\n" , avg - 2*std , avg + 2*std );

		unsigned char *pixels = new unsigned char[ resolution*3 ];
#ifdef NEW_THREADS
		tp.parallel_for( 0 , resolution , [&]( unsigned int , size_t i )
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


template< class Real , typename ... SampleData , unsigned int ... FEMSigs >
void Execute( int argc , char* argv[] , UIntPack< FEMSigs ... > )
{
	static const int Dim = sizeof ... ( FEMSigs );
	typedef UIntPack< FEMSigs ... > Sigs;
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;
	typedef UIntPack< FEMDegreeAndBType< NORMAL_DEGREE , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > NormalSigs;
	static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
	typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
	typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 1 > InterpolationInfo;
	typedef PointStreamNormal< Real , Dim > NormalPointSampleData;
	typedef MultiPointStreamData< Real , SampleData ... > AdditionalPointSampleData;
	typedef MultiPointStreamData< Real , NormalPointSampleData , AdditionalPointSampleData > TotalPointSampleData;
	typedef InputPointStreamWithData< Real , Dim , TotalPointSampleData > InputPointStream;
	typedef TransformedInputPointStreamWithData< Real , Dim , TotalPointSampleData > XInputPointStream;
	std::vector< std::string > comments;
	messageWriter( comments , "************************************************\n" );
	messageWriter( comments , "************************************************\n" );
	messageWriter( comments , "** Running SSD Reconstruction (Version %s) **\n" , VERSION );
	messageWriter( comments , "************************************************\n" );
	messageWriter( comments , "************************************************\n" );


#ifdef NEW_THREADS
	ThreadPool tp( (ThreadPool::ParallelType)ParallelType.value , Threads.value );
#endif // NEW_THREADS

	XForm< Real , Dim+1 > xForm , iXForm;
	if( Transform.set )
	{
		FILE* fp = fopen( Transform.value , "r" );
		if( !fp )
		{
			WARN( "Could not read x-form from: " , Transform.value );
			xForm = XForm< Real , Dim+1 >::Identity();
		}
		else
		{
			for( int i=0 ; i<Dim+1 ; i++ ) for( int j=0 ; j<Dim+1 ; j++ )
			{
				float f;
				if( fscanf( fp , " %f " , &f )!=1 ) ERROR_OUT( "Failed to read xform" );
				xForm(i,j) = (Real)f;
			}
			fclose( fp );
		}
	}
	else xForm = XForm< Real , Dim+1 >::Identity();

	char str[1024];
	for( int i=0 ; params[i] ; i++ )
		if( params[i]->set )
		{
			params[i]->writeValue( str );
			if( strlen( str ) ) messageWriter( comments , "\t--%s %s\n" , params[i]->name , str );
			else                messageWriter( comments , "\t--%s\n" , params[i]->name );
		}

	double startTime = Time();
	Real isoValue = 0;

	FEMTree< Dim , Real > tree( MEMORY_ALLOCATOR_BLOCK_SIZE );
	FEMTreeProfiler< Dim , Real > profiler( tree );

	if( Depth.set && Width.value>0 )
	{
		WARN( "Both --" , Depth.name , " and --" , Width.name , " set, ignoring --" , Width.name );
		Width.value = 0;
	}

#ifdef NEW_CODE
	size_t pointCount;
#else // !NEW_CODE
	int pointCount;
#endif // NEW_CODE

	Real pointWeightSum;
	std::vector< typename FEMTree< Dim , Real >::PointSample >* samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
	std::vector< TotalPointSampleData >* sampleData = NULL;
	DensityEstimator* density = NULL;
	SparseNodeData< Point< Real , Dim > , NormalSigs >* normalInfo = NULL;
	Real targetValue = (Real)0.;

	// Read in the samples (and color data)
	{
		profiler.start();
		InputPointStream* pointStream;
		char* ext = GetFileExtension( In.value );
		sampleData = new std::vector< TotalPointSampleData >();
		std::vector< std::pair< Point< Real , Dim > , TotalPointSampleData > > inCorePoints;
		if( InCore.set )
		{
			InputPointStream *_pointStream;
			if     ( !strcasecmp( ext , "bnpts" ) ) _pointStream = new BinaryInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::ReadBinary );
			else if( !strcasecmp( ext , "ply"   ) ) _pointStream = new    PLYInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::PlyReadProperties() , TotalPointSampleData::PlyReadNum , TotalPointSampleData::ValidPlyReadProperties );
			else                                    _pointStream = new  ASCIIInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::ReadASCII );
			Point< Real , Dim > p;
			TotalPointSampleData d;
			while( _pointStream->nextPoint( p , d ) ) inCorePoints.push_back( std::pair< Point< Real , Dim > , TotalPointSampleData >( p , d ) );
			delete _pointStream;

			pointStream = new MemoryInputPointStreamWithData< Real , Dim , TotalPointSampleData >( inCorePoints.size() , &inCorePoints[0] );
		}
		else
		{
			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::ReadBinary );
			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::PlyReadProperties() , TotalPointSampleData::PlyReadNum , TotalPointSampleData::ValidPlyReadProperties );
			else                                    pointStream = new  ASCIIInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::ReadASCII );
		}
		delete[] ext;
		typename TotalPointSampleData::Transform _xForm( xForm );
		XInputPointStream _pointStream( [&]( Point< Real , Dim >& p , TotalPointSampleData& d ){ p = xForm*p , d = _xForm(d); } , *pointStream );
		if( Width.value>0 ) xForm = GetPointXForm< Real , Dim >( _pointStream , Width.value , (Real)( Scale.value>0 ? Scale.value : 1. ) , Depth.value ) * xForm;
		else                xForm = Scale.value>0 ? GetPointXForm< Real , Dim >( _pointStream , (Real)Scale.value ) * xForm : xForm;
		{
			typename TotalPointSampleData::Transform _xForm( xForm );
			XInputPointStream _pointStream( [&]( Point< Real , Dim >& p , TotalPointSampleData& d ){ p = xForm*p , d = _xForm(d); } , *pointStream );
			auto ProcessDataWithConfidence = [&]( const Point< Real , Dim >& p , TotalPointSampleData& d )
			{
				Real l = (Real)Length( std::get< 0 >( d.data ).data );
				if( !l || l!=l ) return (Real)-1.;
				return (Real)pow( l , Confidence.value );
			};
			auto ProcessData = []( const Point< Real , Dim >& p , TotalPointSampleData& d )
			{
				Real l = (Real)Length( std::get< 0 >( d.data ).data );
				if( !l || l!=l ) return (Real)-1.;
				std::get< 0 >( d.data ).data /= l;
				return (Real)1.;
			};
			if( Confidence.value>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< TotalPointSampleData >( tree.spaceRoot() , _pointStream , Depth.value , *samples , *sampleData , true , tree.nodeAllocators.size() ? tree.nodeAllocators[0] : NULL , tree.initializer() , ProcessDataWithConfidence );
			else                     pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< TotalPointSampleData >( tree.spaceRoot() , _pointStream , Depth.value , *samples , *sampleData , true , tree.nodeAllocators.size() ? tree.nodeAllocators[0] : NULL , tree.initializer() , ProcessData );
		}
		iXForm = xForm.inverse();
		delete pointStream;

#ifdef NEW_CODE
		messageWriter( "Input Points / Samples: %llu / %llu\n" , pointCount , (unsigned long long)samples->size() );
#else // !NEW_CODE
		messageWriter( "Input Points / Samples: %d / %d\n" , pointCount , samples->size() );
#endif // NEW_CODE
		profiler.dumpOutput2( comments , "# Read input into tree:" );
	}
	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value-2;
	if( kernelDepth>Depth.value )
	{
		WARN( KernelDepth.name , " can't be greater than " , Depth.name , ": " , KernelDepth.value , " <= " , Depth.value );
		kernelDepth = Depth.value;
	}

	DenseNodeData< Real , Sigs > solution;
	{
		DenseNodeData< Real , Sigs > constraints;
		InterpolationInfo* iInfo = NULL;
		int solveDepth = Depth.value;

		tree.resetNodeIndices();

		// Get the kernel density estimator
		{
			profiler.start();
#ifdef NEW_THREADS
			density = tree.template setDensityEstimator< WEIGHT_DEGREE >( tp , *samples , kernelDepth , SamplesPerNode.value , 1 );
#else // !NEW_THREADS
			density = tree.template setDensityEstimator< WEIGHT_DEGREE >( *samples , kernelDepth , SamplesPerNode.value , 1 );
#endif // NEW_THREADS
			profiler.dumpOutput2( comments , "#   Got kernel density:" );
		}

		// Transform the Hermite samples into a vector field
		{
			profiler.start();
			normalInfo = new SparseNodeData< Point< Real , Dim > , NormalSigs >();
#ifdef NEW_THREADS
			if( ConfidenceBias.value>0 ) *normalInfo = tree.setNormalField( NormalSigs() , tp , *samples , *sampleData , density , pointWeightSum , [&]( Real conf ){ return (Real)( log( conf ) * ConfidenceBias.value / log( 1<<(Dim-1) ) ); } );
			else                         *normalInfo = tree.setNormalField( NormalSigs() , tp , *samples , *sampleData , density , pointWeightSum );
#else // !NEW_THREADS
			if( ConfidenceBias.value>0 ) *normalInfo = tree.setNormalField( NormalSigs() , *samples , *sampleData , density , pointWeightSum , [&]( Real conf ){ return (Real)( log( conf ) * ConfidenceBias.value / log( 1<<(Dim-1) ) ); } );
			else                         *normalInfo = tree.setNormalField( NormalSigs() , *samples , *sampleData , density , pointWeightSum );
#endif // NEW_THREADS
			profiler.dumpOutput2( comments , "#     Got normal field:" );
			messageWriter( "Point weight / Estimated Area: %g / %g\n" , pointWeightSum , pointCount*pointWeightSum );
		}

		if( !Density.set ) delete density , density = NULL;

		// Trim the tree and prepare for multigrid
		{
			profiler.start();
			constexpr int MAX_DEGREE = NORMAL_DEGREE > Degrees::Max() ? NORMAL_DEGREE : Degrees::Max();
#ifdef NEW_THREADS
			tree.template finalizeForMultigrid< MAX_DEGREE >( tp , FullDepth.value , typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs >( *normalInfo ) , normalInfo , density );
#else // !NEW_THREADS
			tree.template finalizeForMultigrid< MAX_DEGREE >( FullDepth.value , typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs >( *normalInfo ) , normalInfo , density );
#endif // NEW_THREADS
			profiler.dumpOutput2( comments , "#       Finalized tree:" );
		}

		// Free up the normal info [If we don't need it for subsequent iterations.]
		if( normalInfo ) delete normalInfo , normalInfo = NULL;

		// Add the interpolation constraints
		if( ValueWeight.value>0 || GradientWeight.value>0 )
		{
			profiler.start();
#ifdef NEW_THREADS
			if( ExactInterpolation.set ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointAndDataInterpolationInfo< Real , TotalPointSampleData , 1 >( tp , tree , *samples , GetPointer( *sampleData ) , ConstraintDual< Dim , Real , TotalPointSampleData >( targetValue , (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum  ) , SystemDual< Dim , Real , TotalPointSampleData >( (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum ) , true , false );
			else                         iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointAndDataInterpolationInfo< Real , TotalPointSampleData , 1 >( tp , tree , *samples , GetPointer( *sampleData ) , ConstraintDual< Dim , Real , TotalPointSampleData >( targetValue , (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum  ) , SystemDual< Dim , Real , TotalPointSampleData >( (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum ) , true , 1 );
			constraints = tree.initDenseNodeData( Sigs() );
			tree.addInterpolationConstraints( tp , constraints , solveDepth , *iInfo );
#else // !NEW_THREADS
			if( ExactInterpolation.set ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointAndDataInterpolationInfo< Real , TotalPointSampleData , 1 >( tree , *samples , GetPointer( *sampleData ) , ConstraintDual< Dim , Real , TotalPointSampleData >( targetValue , (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum  ) , SystemDual< Dim , Real , TotalPointSampleData >( (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum ) , true , false );
			else                         iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointAndDataInterpolationInfo< Real , TotalPointSampleData , 1 >( tree , *samples , GetPointer( *sampleData ) , ConstraintDual< Dim , Real , TotalPointSampleData >( targetValue , (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum  ) , SystemDual< Dim , Real , TotalPointSampleData >( (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum ) , true , 1 );
			constraints = tree.initDenseNodeData( Sigs() );
			tree.addInterpolationConstraints( constraints , solveDepth , *iInfo );
#endif // NEW_THREADS
			profiler.dumpOutput2( comments , "#Set point constraints:" );
			if( DataX.value<=0 || ( !Colors.set && !Normals.set ) ) delete sampleData , sampleData = NULL;
		}

#ifdef NEW_CODE
		messageWriter( "Leaf Nodes / Active Nodes / Ghost Nodes: %llu / %llu / %llu\n" , (unsigned long long)tree.leaves() , (unsigned long long)tree.nodes() , (unsigned long long)tree.ghostNodes() );
#else // !NEW_CODE
		messageWriter( "Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n" , (int)tree.leaves() , (int)tree.nodes() , (int)tree.ghostNodes() );
#endif // NEW_CODE
		messageWriter( "Memory Usage: %.3f MB\n" , float( MemoryInfo::Usage())/(1<<20) );

		// Solve the linear system
		{
			profiler.start();
			typename FEMTree< Dim , Real >::SolverInfo sInfo;
			sInfo.cgDepth = 0 , sInfo.cascadic = true , sInfo.vCycles = 1 , sInfo.iters = Iters.value , sInfo.cgAccuracy = CGSolverAccuracy.value , sInfo.verbose = Verbose.set , sInfo.showResidual = ShowResidual.set , sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , sInfo.sliceBlockSize = 1;
			sInfo.baseDepth = BaseDepth.value , sInfo.baseVCycles = BaseVCycles.value;
			typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 2 > > F( { 0. , 0. , (double)BiLapWeight.value } );
#ifdef NEW_THREADS
			solution = tree.solveSystem( Sigs() , tp , F , constraints , solveDepth , sInfo , iInfo );
#else // !NEW_THREADS
			solution = tree.solveSystem( Sigs() , F , constraints , solveDepth , sInfo , iInfo );
#endif // NEW_THREADS
			profiler.dumpOutput2( comments , "# Linear system solved:" );
			if( iInfo ) delete iInfo , iInfo = NULL;
		}
	}

	{
		profiler.start();
		double valueSum = 0 , weightSum = 0;
#ifdef NEW_THREADS
		typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( tp , &tree , solution );
		std::vector< double > valueSums( tp.threadNum() , 0 ) , weightSums( tp.threadNum() , 0 );
		tp.parallel_for( 0 , samples->size() , [&]( unsigned int thread , size_t j )
#else // !NEW_THREADS
		typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( &tree , solution );
#pragma omp parallel for reduction( + : valueSum , weightSum )
		for( int j=0 ; j<samples->size() ; j++ )
#endif // NEW_THREADS
		{
			ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
			Real w = sample.weight;
#ifdef NEW_THREADS
			if( w>0 ) weightSums[thread] += w , valueSums[thread] += evaluator.values( sample.data / sample.weight , thread , (*samples)[j].node )[0] * w;
#else // !NEW_THREADS
			if( w>0 ) weightSum += w , valueSum += evaluator.values( sample.data / sample.weight , omp_get_thread_num() , (*samples)[j].node )[0] * w;
#endif // NEW_THREADS
		}
#ifdef NEW_THREADS
		);
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) valueSum += valueSums[t] , weightSum += weightSums[t];
#endif // NEW_THREADS
		isoValue = (Real)( valueSum / weightSum );
		if( DataX.value<=0 || ( !Colors.set && !Normals.set ) ) delete samples , samples = NULL;
		profiler.dumpOutput( "Got average:" );
		messageWriter( "Iso-Value: %e = %g / %g\n" , isoValue , valueSum , weightSum );
	}
	if( Tree.set )
	{
		FILE* fp = fopen( Tree.value , "wb" );
		if( !fp ) ERROR_OUT( "Failed to open file for writing: " , Tree.value );
		FEMTree< Dim , Real >::WriteParameter( fp );
		DenseNodeData< Real , Sigs >::WriteSignatures( fp );
		tree.write( fp );
		solution.write( fp );
		fclose( fp );
	}

	if( Grid.set )
	{
		int res = 0;
		profiler.start();
#ifdef NEW_THREADS
		Pointer( Real ) values = tree.template regularGridEvaluate< true >( tp , solution , res , -1 , PrimalGrid.set );
#else // !NEW_THREADS
		Pointer( Real ) values = tree.template regularGridEvaluate< true >( solution , res , -1 , PrimalGrid.set );
#endif // NEW_THREADS
		size_t resolution = 1;
		for( int d=0 ; d<Dim ; d++ ) resolution *= res;
#ifdef NEW_THREADS
		tp.parallel_for( 0 , resolution , [&]( unsigned int , size_t i ){ values[i] -= isoValue; } );
#else // !NEW_THREADS
#pragma omp parallel for
		for( int i=0 ; i<resolution ; i++ ) values[i] -= isoValue;
#endif // NEW_THREADS
		profiler.dumpOutput( "Got grid:" );
#ifdef NEW_THREADS
		WriteGrid< Real , DIMENSION >( tp , values , res , Grid.value );
#else // !NEW_THREADS
		WriteGrid< Real , DIMENSION >( values , res , Grid.value );
#endif // NEW_THREADS
		DeletePointer( values );
		if( Verbose.set )
		{
			printf( "Transform:\n" );
			for( int i=0 ; i<Dim+1 ; i++ )
			{
				printf( "\t" );
				for( int j=0 ; j<Dim+1 ; j++ ) printf( " %f" , iXForm(j,i) );
				printf( "\n" );
			}
		}
	}

	if( Out.set )
	{
		if( Normals.set )
		{
			if( Density.set )
			{
				typedef PlyVertexWithData< Real , Dim , MultiPointStreamData< Real , PointStreamNormal< Real , Dim > , PointStreamValue< Real > , AdditionalPointSampleData > > Vertex;
				std::function< void ( Vertex& , Point< Real , Dim > , Real , TotalPointSampleData ) > SetVertex = []( Vertex& v , Point< Real , Dim > p , Real w , TotalPointSampleData d ){ v.point = p , std::get< 0 >( v.data.data ) = std::get< 0 >( d.data ) , std::get< 1 >( v.data.data ).data = w , std::get< 2 >( v.data.data ) = std::get< 1 >( d.data ); };
#ifdef NEW_THREADS
				ExtractMesh< Vertex >( UIntPack< FEMSigs ... >() , std::tuple< SampleData ... >() , tp , tree , solution , isoValue , samples , sampleData , density , SetVertex , comments , iXForm );
#else // !NEW_THREADS
				ExtractMesh< Vertex >( UIntPack< FEMSigs ... >() , std::tuple< SampleData ... >() , tree , solution , isoValue , samples , sampleData , density , SetVertex , comments , iXForm );
#endif // NEW_THREADS
			}
			else
			{
				typedef PlyVertexWithData< Real , Dim , MultiPointStreamData< Real , PointStreamNormal< Real , Dim > , AdditionalPointSampleData > > Vertex;
				std::function< void ( Vertex& , Point< Real , Dim > , Real , TotalPointSampleData ) > SetVertex = []( Vertex& v , Point< Real , Dim > p , Real w , TotalPointSampleData d ){ v.point = p , std::get< 0 >( v.data.data ) = std::get< 0 >( d.data ) , std::get< 1 >( v.data.data ) = std::get< 1 >( d.data ); };
#ifdef NEW_THREADS
				ExtractMesh< Vertex >( UIntPack< FEMSigs ... >() , std::tuple< SampleData ... >() , tp , tree , solution , isoValue , samples , sampleData , density , SetVertex , comments , iXForm );
#else // !NEW_THREADS
				ExtractMesh< Vertex >( UIntPack< FEMSigs ... >() , std::tuple< SampleData ... >() , tree , solution , isoValue , samples , sampleData , density , SetVertex , comments , iXForm );
#endif // NEW_THREADS
			}
		}
		else
		{
			if( Density.set )
			{
				typedef PlyVertexWithData< Real , Dim , MultiPointStreamData< Real , PointStreamValue< Real > , AdditionalPointSampleData > > Vertex;
				std::function< void ( Vertex& , Point< Real , Dim > , Real , TotalPointSampleData ) > SetVertex = []( Vertex& v , Point< Real , Dim > p , Real w , TotalPointSampleData d ){ v.point = p , std::get< 0 >( v.data.data ).data = w , std::get< 1 >( v.data.data ) = std::get< 1 >( d.data ); };
#ifdef NEW_THREADS
				ExtractMesh< Vertex >( UIntPack< FEMSigs ... >() , std::tuple< SampleData ... >() , tp , tree , solution , isoValue , samples , sampleData , density , SetVertex , comments , iXForm );
#else // !NEW_THREADS
				ExtractMesh< Vertex >( UIntPack< FEMSigs ... >() , std::tuple< SampleData ... >() , tree , solution , isoValue , samples , sampleData , density , SetVertex , comments , iXForm );
#endif // NEW_THREADS
			}
			else
			{
				typedef PlyVertexWithData< Real , Dim , MultiPointStreamData< Real , AdditionalPointSampleData > > Vertex;
				std::function< void ( Vertex& , Point< Real , Dim > , Real , TotalPointSampleData ) > SetVertex = []( Vertex& v , Point< Real , Dim > p , Real w , TotalPointSampleData d ){ v.point = p , std::get< 0 >( v.data.data ) = std::get< 1 >( d.data ); };
#ifdef NEW_THREADS
				ExtractMesh< Vertex >( UIntPack< FEMSigs ... >() , std::tuple< SampleData ... >() , tp , tree , solution , isoValue , samples , sampleData , density , SetVertex , comments , iXForm );
#else // !NEW_THREADS
				ExtractMesh< Vertex >( UIntPack< FEMSigs ... >() , std::tuple< SampleData ... >() , tree , solution , isoValue , samples , sampleData , density , SetVertex , comments , iXForm );
#endif // NEW_THREADS
			}
		}
		if( sampleData ){ delete sampleData ; sampleData = NULL; }
	}
	if( density ) delete density , density = NULL;
	messageWriter( comments , "#          Total Solve: %9.1f (s), %9.1f (MB)\n" , Time()-startTime , FEMTree< Dim , Real >::MaxMemoryUsage() );
}

#ifndef FAST_COMPILE
template< unsigned int Dim , class Real , typename ... SampleData >
void Execute( int argc , char* argv[] )
{
	switch( BType.value )
	{
	case BOUNDARY_FREE+1:
	{
		switch( Degree.value )
		{
			case 2: return Execute< Real , SampleData ... >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BOUNDARY_FREE >::Signature >() );
			case 3: return Execute< Real , SampleData ... >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BOUNDARY_FREE >::Signature >() );
//			case 4: return Execute< Real , SampleData ... >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BOUNDARY_FREE >::Signature >() );
			default: ERROR_OUT( "Only B-Splines of degree 2 - 3 are supported" );
		}
	}
	case BOUNDARY_NEUMANN+1:
	{
		switch( Degree.value )
		{
			case 2: return Execute< Real , SampleData ... >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BOUNDARY_NEUMANN >::Signature >() );
			case 3: return Execute< Real , SampleData ... >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BOUNDARY_NEUMANN >::Signature >() );
//			case 4: return Execute< Real , SampleData ... >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BOUNDARY_NEUMANN >::Signature >() );
			default: ERROR_OUT( "Only B-Splines of degree 2 - 3 are supported" );
		}
	}
	case BOUNDARY_DIRICHLET+1:
	{
		switch( Degree.value )
		{
			case 2: return Execute< Real , SampleData ... >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BOUNDARY_DIRICHLET >::Signature >() );
			case 3: return Execute< Real , SampleData ... >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BOUNDARY_DIRICHLET >::Signature >() );
//			case 4: return Execute< Real , SampleData ... >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BOUNDARY_DIRICHLET >::Signature >() );
			default: ERROR_OUT( "Only B-Splines of degree 2 - 3 are supported" );
		}
	}
	default: ERROR_OUT( "Not a valid boundary type: " , BType.value );
	}
}
#endif // !FAST_COMPILE

int main( int argc , char* argv[] )
{
	Timer timer;
#ifdef ARRAY_DEBUG
	WARN( "Array debugging enabled" );
#endif // ARRAY_DEBUG

	cmdLineParse( argc-1 , &argv[1] , params );
	if( MaxMemoryGB.value>0 ) SetPeakMemoryMB( MaxMemoryGB.value<<10 );
#ifdef NEW_THREADS
	ThreadPool::DefaultChunkSize = ThreadChunkSize.value;
	ThreadPool::DefaultSchedule = (ThreadPool::ScheduleType)ScheduleType.value;
#else // !NEW_THREADS
	omp_set_num_threads( Threads.value > 1 ? Threads.value : 1 );
#endif // NEW_THREADS
	messageWriter.echoSTDOUT = Verbose.set;

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return 0;
	}
	if( GradientWeight.value<=0 ) ERROR_OUT( "Gradient weight must be positive: " , GradientWeight.value , "> 0" );
	if( BiLapWeight.value<=0 ) ERROR_OUT( "Bi-Laplacian weight must be positive: " , BiLapWeight.value , " > 0" );
	if( DataX.value<=0 ) Normals.set = Colors.set = false;
	if( BaseDepth.value>FullDepth.value )
	{
		if( BaseDepth.set ) WARN( "Base depth must be smaller than full depth: " , BaseDepth.value , " <= " , FullDepth.value );
		BaseDepth.value = FullDepth.value;
	}
	ValueWeight.value    *= (float)BaseSSDWeights[0];
	GradientWeight.value *= (float)BaseSSDWeights[1];
	BiLapWeight.value    *= (float)BaseSSDWeights[2];

#ifdef USE_DOUBLE
	typedef double Real;
#else // !USE_DOUBLE
	typedef float  Real;
#endif // USE_DOUBLE

#ifdef FAST_COMPILE
	static const int Degree = DEFAULT_FEM_DEGREE;
	static const BoundaryType BType = DEFAULT_FEM_BOUNDARY;
	typedef IsotropicUIntPack< DIMENSION , FEMDegreeAndBType< Degree , BType >::Signature > FEMSigs;
#ifdef NEW_CODE
	WARN( "Compiled for degree-" , Degree , ", boundary-" , BoundaryNames[ BType ] , ", " , sizeof(Real)==4 ? "single" : "double" , "-precision _only_" );
	if( Colors.set ) Execute< Real , PointStreamColor< Real > >( argc , argv , FEMSigs() );
#else // !NEW_CODE
	WARN( "Compiled for degree-" , Degree , ", boundary-" , BoundaryNames[ BType ] , ", " , sizeof(DefaultFloatType)==4 ? "single" : "double" , "-precision _only_" );
	if( Colors.set ) Execute< Real , PointStreamColor< DefaultFloatType > >( argc , argv , FEMSigs() );
#endif // NEW_CODE
	else             Execute< Real >( argc , argv , FEMSigs() );
#else // !FAST_COMPILE
	if( Colors.set ) Execute< DIMENSION , Real , PointStreamColor< float > >( argc , argv );
	else             Execute< DIMENSION , Real >( argc , argv );
#endif // FAST_COMPILE
	if( Performance.set )
	{
		printf( "Time (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
		printf( "Peak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );
	}
	return EXIT_SUCCESS;
}
