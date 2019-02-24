/*
Copyright (c) 2013, Michael Kazhdan
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

#undef ARRAY_DEBUG
#undef FAST_COMPILE
#undef USE_DOUBLE
#define DIMENSION 2
#define USE_DEEP_TREE_NODES
#define ROW_BLOCK_SIZE 16
#define DEFAULT_FEM_DEGREE 1

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include "Image.h"
#include "MyMiscellany.h"
#include "Array.h"
#include "CmdLineParser.h"
#include "Geometry.h"
#include "FEMTree.h"

cmdLineParameterArray< char* , 2 >
	In( "in" );
cmdLineParameter< char* >
	Out( "out" );
cmdLineParameter< int >
#ifdef FAST_COMPILE
#else // !FAST_COMPILE
	Degree( "degree" , DEFAULT_FEM_DEGREE ) ,
#endif // FAST_COMPILE
#ifdef NEW_THREADS
	Threads( "threads" , ThreadPool::DefaultThreadNum ) ,
#else // !NEW_THREADS
	Threads( "threads" , omp_get_num_procs() ) ,
#endif // NEW_THREADS
	MaxMemoryGB( "maxMemory" , 0 ) ,
	GSIterations( "iters" , 8 ) ,
	FullDepth( "fullDepth" , 6 ) ,
	BaseDepth( "baseDepth" , 6 ) ,
	BaseVCycles( "baseVCycles" , 4 );
cmdLineReadable
	Verbose( "verbose" ) ,
	ShowResidual( "residual" ) ,
	Performance( "performance" );
cmdLineParameter< float >
	WeightScale   ( "wScl", 0.125f ) ,
	WeightExponent( "wExp" , 6.f );

cmdLineReadable* params[] =
{
	&In , &Out , &Threads , &Verbose , &ShowResidual , &GSIterations , &FullDepth ,
	&BaseDepth , &BaseVCycles ,
	&WeightScale , &WeightExponent ,
	&Performance ,
	&MaxMemoryGB ,
#if !defined( FAST_COMPILE )
	&Degree , 
#endif // !FAST_COMPILE
	NULL
};

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input color / labels>\n" , In.name );
	printf( "\t[--%s <ouput stitched image>]\n" , Out.name );
#if !defined( FAST_COMPILE )
	printf( "\t[--%s <b-spline degree>=%d]\n" , Degree.name , Degree.value );
#endif // !FAST_COMPILE
	printf( "\t[--%s <GS iterations>=%d]\n" , GSIterations.name , GSIterations.value );
	printf( "\t[--%s <full depth>=%d]\n" , FullDepth.name , FullDepth.value );
	printf( "\t[--%s <parallelization threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <successive under-relaxation scale>=%f]\n", WeightScale.name , WeightScale.value );
	printf( "\t[--%s <successive under-relaxation exponent>=%f]\n", WeightExponent.name , WeightExponent.value );
	printf( "\t[--%s <maximum memory (in GB)>=%d]\n" , MaxMemoryGB.name , MaxMemoryGB.value );
	printf( "\t[--%s]\n" , Performance.name );
	printf( "\t[--%s <coarse MG solver depth>=%d]\n" , BaseDepth.name , BaseDepth.value );
	printf( "\t[--%s <coarse MG solver v-cycles>=%d]\n" , BaseVCycles.name , BaseVCycles.value );
	printf( "\t[--%s]\n" , ShowResidual.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

struct RGBPixel
{
	unsigned char rgb[3];
	unsigned char& operator[]( int idx ){ return rgb[idx]; }
	const unsigned char& operator[]( int idx ) const { return rgb[idx]; }
	int mask( void ) const
	{
		if( rgb[0]==255 && rgb[1]==255 && rgb[2]==255 ) return -1;
		else return ( (int)rgb[0] )<<16 | ( (int)rgb[1] )<<8 | ( (int)rgb[2] );
	}

	RGBPixel( void ){ rgb[0] = rgb[1] = rgb[2] = 0; }
	RGBPixel( double r , double g , double b )
	{
		rgb[0] = (unsigned char)( std::max< int >( 0 , std::min< int >( 255 , (int)( r*255 ) ) ) );
		rgb[1] = (unsigned char)( std::max< int >( 0 , std::min< int >( 255 , (int)( g*255 ) ) ) );
		rgb[2] = (unsigned char)( std::max< int >( 0 , std::min< int >( 255 , (int)( b*255 ) ) ) );
	}
	RGBPixel( float r , float g , float b )
	{
		rgb[0] = (unsigned char)( std::max< int >( 0 , std::min< int >( 255 , (int)( r*255 ) ) ) );
		rgb[1] = (unsigned char)( std::max< int >( 0 , std::min< int >( 255 , (int)( g*255 ) ) ) );
		rgb[2] = (unsigned char)( std::max< int >( 0 , std::min< int >( 255 , (int)( b*255 ) ) ) );
	}

	template< class Real >
	static Point< Real , 3 > ToPoint( RGBPixel rgb )
	{
		Point< Real , 3 > p;
		for( int c=0 ; c<3 ; c++ ) p[c] = (Real)( ( (double)rgb[c] ) / 255. );
		return p;
	}
};

void WriteImage( char* fileName , RGBPixel* pixels , int w , int h )
{
	unsigned int _w = w , _h = h , _c = 3;
	ImageWriter::Write( fileName , (const unsigned char*)pixels , _w , _h , _c );
}

struct Profiler
{
	double t;
	Profiler( void ){ t = Time(); }
	void print( bool newLine=false ) const
	{
		printf( "%.2f (s) ; %d (MB)" , Time()-t , MemoryInfo::PeakMemoryUsageMB() );
		if( newLine ) printf( "\n" );
	}
};

template< unsigned int Colors >
void ReadAndWrite( ImageReader* pixels , ImageReader* labels , ImageWriter* output )
{
	RGBPixel* pixelRow = new RGBPixel[ pixels->width() ];
	RGBPixel* labelRow = new RGBPixel[ labels->width() ];
	for( unsigned int r=0 ; r<pixels->height() ; r++ )
	{
		if( Verbose.set ) printf( "%d / %d       \r" , r ,pixels->height() );
		pixels->nextRow( (unsigned char*)pixelRow );
		labels->nextRow( (unsigned char*)labelRow );
		output->nextRow( (unsigned char*)pixelRow );
	}
	if( Verbose.set ) printf( "\n" );
	delete[] pixelRow;
	delete[] labelRow;
}

template< class Real , unsigned int Colors >
struct BufferedImageDerivativeStream : public FEMTreeInitializer< DIMENSION , Real >::template DerivativeStream< Point< Real , Colors > >
{
	BufferedImageDerivativeStream( const unsigned int resolution[] , ImageReader* pixels , ImageReader* labels ) : _pixels( pixels ) , _labels( labels )
	{
		memcpy( _resolution , resolution , sizeof( unsigned int ) * DIMENSION );
		for( int i=0 ; i<3 ; i++ )
		{
			_pixelRows[i] = new RGBPixel[ _resolution[0] ];
			_labelRows[i] = new RGBPixel[ _resolution[0] ];
			_maskRows [i] = new      int[ _resolution[0] ];
		}
		if( pixels->channels()!=3 && pixels->channels()!=1 ) ERROR_OUT( "Pixel input must have 1 or 3 channels: " , pixels->channels() );
		if( labels->channels()!=3 && labels->channels()!=1 ) ERROR_OUT( "Label input must have 1 or 3 channels: " , labels->channels() );
		__pixelRow = pixels->channels()==3 ? NULL : new unsigned char[ _resolution[0] ];
		__labelRow = labels->channels()==3 ? NULL : new unsigned char[ _resolution[0] ];
		_r = -2 ; prefetch();
		_r = -1 ; prefetch();
		_c = _r = _dir = 0;
	}
	~BufferedImageDerivativeStream( void )
	{
		for( int i=0 ; i<3 ; i++ ) delete[] _pixelRows[i] , delete[] _labelRows[i] , delete[] _maskRows[i];
		if( __pixelRow ) delete[] __pixelRow;
		if( __labelRow ) delete[] __labelRow;
	}
	void resolution( unsigned int res[] ) const { memcpy( res , _resolution , sizeof(_resolution) ); }

	void advance( void ){ _c = _dir = 0 , _r++; }
	void prefetch( void )
	{
#ifdef NEW_THREADS
		static ThreadPool tp;
#endif // NEW_THREADS
		if( _r+2<(int)_resolution[1] )
		{
			int j = (_r+2)%3;
			RGBPixel *pixelRow = _pixelRows[j] , *labelRow = _labelRows[j];
			int *maskRow = _maskRows[j];
			if( _pixels->channels()==3 ) _pixels->nextRow( (unsigned char*)pixelRow );
			else
			{
				_pixels->nextRow( __pixelRow );
				for( int i=0 ; i<(int)_resolution[0] ; i++ ) pixelRow[i][0] = pixelRow[i][1] = pixelRow[i][2] = __pixelRow[i];
			}
			if( _labels->channels()==3 ) _labels->nextRow( (unsigned char*)labelRow );
			else
			{
				_labels->nextRow( __labelRow );
				for( int i=0 ; i<(int)_resolution[0] ; i++ ) labelRow[i][0] = labelRow[i][1] = labelRow[i][2] = __labelRow[i];
			}
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			tp.parallel_for( 0 , _resolution[0] , [&]( const ThreadPool::ThreadNum & , size_ i ){ maskRow[i] = labelRow[i].mask(); } );
#else // !NEW_THREAD_NUM
			tp.parallel_for( 0 , _resolution[0] , [&]( unsigned int , size_ i ){ maskRow[i] = labelRow[i].mask(); } );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
			for( int i=0 ; i<(int)_resolution[0] ; i++ ) maskRow[i] = labelRow[i].mask();
#endif // NEW_THREADS
		}
	}

	bool nextDerivative( unsigned int idx[] , unsigned int& dir , Point< Real , Colors >& dValue )
	{
		const RGBPixel *pixelRow1 = _pixelRows[_r%3] , *pixelRow2 = _pixelRows[(_r+1)%3];
		const int *maskRow1 = _maskRows[_r%3] , *maskRow2 = _maskRows[(_r+1)%3];
		if( _dir==0 )
		{
			for( ; _c<(int)_resolution[0]-1 ; _c++ )
			{
				if( maskRow1[_c]!=maskRow1[_c+1] && maskRow1[_c]>=0 && maskRow1[_c+1]>=0 )
				{
					idx[0] = _c , idx[1] = _r , dir = _dir;
					dValue = RGBPixel::ToPoint< Real >( pixelRow1[_c+1] ) - RGBPixel::ToPoint< Real >( pixelRow1[_c] );
					_c++;
					return true;
				}
			}
			_dir = 1 , _c = 0;
		}
		if( _dir==1 )
		{
			if( _r+1<(int)_resolution[1] )
			{
				for( ; _c<(int)_resolution[0] ; _c++ )
				{
					if( maskRow1[_c]!=maskRow2[_c] && maskRow1[_c]>=0 && maskRow2[_c]>=0 )
					{
						idx[0] = _c , idx[1] = _r , dir = _dir;
						dValue = RGBPixel::ToPoint< Real >( pixelRow2[_c] ) - RGBPixel::ToPoint< Real >( pixelRow1[_c] );
						_c++;
						return true;
					}
				}
			}
		}
		return false;
	}
protected:
	int _r , _c , _dir;
	unsigned int _resolution[DIMENSION];
	ImageReader *_pixels , *_labels;
	RGBPixel *_pixelRows[3] , *_labelRows[3];
	unsigned char *__pixelRow , *__labelRow;
	int* _maskRows[3];
};

template< typename Real , unsigned int Degree >
void _Execute( void )
{
#ifdef NEW_THREADS
	ThreadPool tp;
#endif // NEW_THREADS
	int w , h;
	{
		unsigned int _w , _h , _c;
		ImageReader::GetInfo( In.values[0] , _w , _h , _c );
		w = _w , h = _h;
		ImageReader::GetInfo( In.values[1] , _w , _h , _c );
		if( w!=_w || h!=_h ) ERROR_OUT( "Pixel and label dimensions don't match: " , _w , " x " , _h , " != " , w , " x " , h );
	}
	if( Verbose.set ) printf( "Resolution: %d x %d\n" , w , h );

	static const unsigned int Dim = DIMENSION;
	static const unsigned int Colors = 3;
	static const unsigned int FEMSig = FEMDegreeAndBType< Degree , BOUNDARY_NEUMANN >::Signature;

	FEMTree< Dim , Real > tree( MEMORY_ALLOCATOR_BLOCK_SIZE );
	std::vector< NodeSample< Dim , Point< Real , Colors > > > derivatives[Dim];
	int maxDepth;
	DenseNodeData< Point< Real , Colors > , IsotropicUIntPack< Dim , FEMSig > > constraints;
	DenseNodeData< Point< Real , Colors > , IsotropicUIntPack< Dim , FEMSig > > solution;
	{
		Profiler p;
		ImageReader* pixels = ImageReader::Get( In.values[0] );
		ImageReader* labels = ImageReader::Get( In.values[1] );
		unsigned int resolution[] = { (unsigned int )w , (unsigned int )h };
		BufferedImageDerivativeStream< Real , Colors > dStream( resolution , pixels , labels );
		for( int j=0 ; j<h ; j++ )
		{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			tp.parallel_for( 0 , 2 , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( 0 , 2 , [&]( unsigned int , size_t i )
#endif // NEW_THREAD_NUM
			{
				if( i==0 ) dStream.prefetch();
				else maxDepth = FEMTreeInitializer< Dim , Real >::template Initialize< (Degree&1)==0 , Point< Real , Colors > >( tree.spaceRoot() , dStream , derivatives , tree.nodeAllocator , tree.initializer() );
#ifdef NEW_THREAD_LOOP
			} , 1 , 1
#else // !NEW_THREAD_LOOP
			} , 1
#endif // NEW_THREAD_LOOP
			);
#else // !NEW_THREADS
#pragma omp parallel sections
			{
#pragma omp section
				{
					dStream.prefetch();
				}
#pragma omp section
				{
					maxDepth = FEMTreeInitializer< Dim , Real >::template Initialize< (Degree&1)==0 , Point< Real , Colors > >( tree.spaceRoot() , dStream , derivatives , tree.nodeAllocator , tree.initializer() );
				}
			}
#endif // NEW_THREADS
			dStream.advance();
		}
		delete pixels;
		delete labels;
		{
			std::vector< typename FEMTree< Dim , Real >::FEMTreeNode* > nodes;
			nodes.reserve( derivatives[0].size() + derivatives[1].size() );
			for( int i=0 ; i<derivatives[0].size() ; i++ ) nodes.push_back( derivatives[0][i].node );
			for( int i=0 ; i<derivatives[1].size() ; i++ ) nodes.push_back( derivatives[1][i].node );
			tree.template thicken< 1 , 0 >( &nodes[0] , (int)nodes.size() );
		}
#ifdef NEW_CODE
#ifdef NEW_THREADS
		tree.template finalizeForMultigrid< Degree >( tp , FullDepth.value , []( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* ){ return true; } );
#else // !NEW_THREADS
		tree.template finalizeForMultigrid< Degree >( FullDepth.value , []( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* ){ return true; } );
#endif // NEW_THREADS
#else // !NEW_CODE
		tree.template finalizeForMultigrid< Degree >( FullDepth.value , []( const RegularTreeNode< Dim , FEMTreeNodeData >* ){ return true; } );
#endif // NEW_CODE
		if( Verbose.set )
		{
#ifdef NEW_CODE
			printf( "Valid FEM Nodes / Edges: %llu %llu\n" , (unsigned long long)tree.validFEMNodes( IsotropicUIntPack< Dim , FEMSig >() ) , (unsigned long long)( derivatives[0].size() + derivatives[1].size() ) );
#else // !NEW_CODE
			printf( "Valid FEM Nodes / Edges: %d %d\n" , (int)tree.validFEMNodes( IsotropicUIntPack< Dim , FEMSig >() ) , (int)( derivatives[0].size() + derivatives[1].size() ) );
#endif // NEW_CODE
			printf( "Set tree [%d]: " , maxDepth ) , p.print( true );
		}
	}

	{
		Profiler p;
		constraints = tree.template initDenseNodeData< Point< Real , Colors > >( IsotropicUIntPack< Dim , FEMSig >() );
		static const unsigned int DFEMSig = FEMSignature< FEMSig >::DSignature();
		// Generate the partial-x constraints
		{
			typedef UIntPack< DFEMSig , FEMSig > CSignature;
			typedef IsotropicUIntPack< 2 , 0 > CDerivative;
			typedef UIntPack< FEMSig , FEMSig > FEMSignature;
			typedef UIntPack< 1 , 0 > FEMDerivative;
			SparseNodeData< Point< Real , Colors > , CSignature > partialX;
			for( int i=0 ; i<derivatives[0].size() ; i++ ) partialX[ derivatives[0][i].node ] = -derivatives[0][i].data * (1<<maxDepth);

			unsigned int derivatives1[] = { 1 , 0 } , derivatives2[] = { 0 , 0 };
			typename FEMIntegrator::template Constraint< FEMSignature , FEMDerivative , CSignature , CDerivative , 1 > F;
			F.weights[0][ TensorDerivatives< FEMDerivative >::Index( derivatives1 ) ][ TensorDerivatives< CDerivative >::Index( derivatives2 ) ] = 1;
#ifdef NEW_THREADS
			tree.addFEMConstraints( tp , F , partialX , constraints , maxDepth );
#else // !NEW_THREADS
			tree.addFEMConstraints( F , partialX , constraints , maxDepth );
#endif // NEW_THREADS
		}
		// Generate the partial-y constraints
		{
			typedef UIntPack< FEMSig , DFEMSig > CSignature;
			typedef IsotropicUIntPack< 2 , 0 > CDerivative;
			typedef UIntPack< FEMSig , FEMSig > FEMSignature;
			typedef UIntPack< 0 , 1 > FEMDerivative;
			SparseNodeData< Point< Real , Colors > , CSignature > partialY;
			for( int i=0 ; i<derivatives[1].size() ; i++ ) partialY[ derivatives[1][i].node ] = -derivatives[1][i].data * (1<<maxDepth);

			unsigned int derivatives1[] = { 0 , 1 } , derivatives2[] = { 0 , 0 };
			typename FEMIntegrator::template Constraint< FEMSignature , FEMDerivative , CSignature , CDerivative , 1 > F;
			F.weights[0][ TensorDerivatives< FEMDerivative >::Index( derivatives1 ) ][ TensorDerivatives< CDerivative >::Index( derivatives2 ) ] = 1;
#ifdef NEW_THREADS
			tree.addFEMConstraints( tp , F , partialY , constraints , maxDepth );
#else // !NEW_THREADS
			tree.addFEMConstraints( F , partialY , constraints , maxDepth );
#endif // NEW_THREADS
		}
		if( Verbose.set ) printf( "Set constraints: " ) , p.print( true );
	}
	// Solve the system
	{
		Profiler p;
		double t = Time();
		solution = tree.template initDenseNodeData< Point< Real , Colors > >( IsotropicUIntPack< Dim , FEMSig >() );
		typename FEMTree< Dim , Real >::SolverInfo sInfo;
		sInfo.cgDepth = 0 , sInfo.cascadic = false , sInfo.vCycles = 1 , sInfo.cgAccuracy = 0 , sInfo.verbose = Verbose.set , sInfo.showResidual = ShowResidual.set , sInfo.showGlobalResidual = false , sInfo.sliceBlockSize = ROW_BLOCK_SIZE;
		sInfo.baseDepth = BaseDepth.value , sInfo.baseVCycles = BaseVCycles.value;
		sInfo.iters = GSIterations.value;

		sInfo.useSupportWeights = true;
		sInfo.sorRestrictionFunction = [&] ( Real w , Real ){ return (Real)( WeightScale.value * pow( w , WeightExponent.value ) ); };
		sInfo.wCycle = false;
		typename FEMIntegrator::template System< IsotropicUIntPack< Dim , FEMSig > , IsotropicUIntPack< Dim , 1 > > F( { 0. , 1. } );
		DenseNodeData< Point< Real , Colors > , IsotropicUIntPack< Dim , FEMSig > > _constraints = tree.template initDenseNodeData< Point< Real , Colors > >( IsotropicUIntPack< Dim , FEMSig >() );
#ifdef NEW_THREADS
		tree.solveSystem( IsotropicUIntPack< Dim , FEMSig >() , tp , F , constraints , solution , Point< Real , Colors >::Dot , maxDepth , sInfo );
#else // !NEW_THREADS
		tree.solveSystem( IsotropicUIntPack< Dim , FEMSig >() , F , constraints , solution , Point< Real , Colors >::Dot , maxDepth , sInfo );
#endif // NEW_THREADS
		if( Verbose.set ) printf( "Solved system: " ) , p.print( true );
	}

	Point< Real , Colors > average;
	{
		Profiler p;
		Real begin[] = { 0 , 0 } , end[] = { (Real)w/(1<<maxDepth) , (Real)h/(1<<maxDepth) };
#ifdef NEW_THREADS
		average = tree.average( tp , solution , begin , end );
#else // !NEW_THREADS
		average = tree.average( solution , begin , end );
#endif // NEW_THREADS
		if( Verbose.set ) printf( "Got average: " ) , p.print( true );
	}
	// Stitch the image
	if( Out.set )
	{
		Profiler p;
		int begin[2] , end[2];
		ImageReader* in = ImageReader::Get( In.values[0] );
		ImageWriter* out = ImageWriter::Get( Out.value , w , h , 3 );

		RGBPixel *inRows[2] , *outRows[2];
		unsigned char* inRow = NULL;
		for( int i=0 ; i<2 ; i++ ) inRows[i] = new RGBPixel[w*ROW_BLOCK_SIZE] , outRows[i] = new RGBPixel[w*ROW_BLOCK_SIZE];
		if( in->channels()==1 ) inRow = new unsigned char[w];

		auto FetchInput = [&]( unsigned int block )
		{
#ifdef NEW_THREADS
			static ThreadPool tp;
#endif // NEW_THREADS
			int rStart = block*ROW_BLOCK_SIZE;
			int rEnd = rStart + ROW_BLOCK_SIZE < h ? rStart + ROW_BLOCK_SIZE : h;
			for( int r=rStart , rr=0 ; r<rEnd ; r++ , rr++ )
			{
				if( in->channels()==3 ) in->nextRow( (unsigned char*)( inRows[block&1] + rr*w ) );
				else
				{
					in->nextRow( inRow );
					RGBPixel *_inRow = inRows[block&1] + rr*w;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
					tp.parallel_for( 0 , w , [&]( const ThreadPool::ThreadNum & , size_t i ){ _inRow[i][0] = _inRow[i][1] = _inRow[i][2] = inRow[i]; } );
#else // !NEW_THREAD_NUM
					tp.parallel_for( 0 , w , [&]( unsigned int , size_t i ){ _inRow[i][0] = _inRow[i][1] = _inRow[i][2] = inRow[i]; } );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
					for( int i=0 ; i<w ; i++ ) _inRow[i][0] = _inRow[i][1] = _inRow[i][2] = inRow[i];
#endif // NEW_THREADS
				}
			}
		};
		auto SetOutput = [&]( unsigned int block )
		{
			int rStart = block*ROW_BLOCK_SIZE;
			int rEnd = rStart + ROW_BLOCK_SIZE < h ? rStart + ROW_BLOCK_SIZE : h;
			out->nextRows( (unsigned char*)outRows[block&1] , rEnd-rStart );
		};
		int blockNum = ( h + ROW_BLOCK_SIZE - 1 ) / ROW_BLOCK_SIZE;

		// Prefetch the first block
		FetchInput( 0 );
#ifdef NEW_THREADS
#else // !NEW_THREADS
		omp_set_nested( true );
#endif // NEW_THREADS
		for( int rStart=0 , block=0 ; rStart<h ; rStart+=ROW_BLOCK_SIZE , block++ )
		{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			tp.parallel_for( 0 , 3 , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( 0 , 3 , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
			{
				switch( i )
				{
				case 0: if( block<blockNum ) FetchInput( block+1 ) ; break;
				case 1: if( block>0 ) SetOutput( block-1 ) ; break;
				case 2:
				{
					static ThreadPool tp;
					RGBPixel *_inRows = inRows[block&1] , *_outRows = outRows[block&1];
					int rEnd = rStart + ROW_BLOCK_SIZE < h ? rStart + ROW_BLOCK_SIZE : h;

					// Expand the next block of values
					begin[0] = 0 , begin[1] = rStart , end[0] = w , end[1] = rEnd;
					Pointer( Point< Real , Colors > ) outBlock = tree.template regularGridUpSample< true >( tp , solution , begin , end );
					int size = (rEnd-rStart)*w;
					tp.parallel_for( 0 , size , [&]( unsigned int , size_t ii )
					{
						Point< Real , Colors > c = Point< Real , Colors >( _inRows[ii][0] , _inRows[ii][1] , _inRows[ii][2] ) / 255;
						c += outBlock[ii] - average;
						_outRows[ii] = RGBPixel( c[0] , c[1] , c[2] );
					}
					);
					DeletePointer( outBlock );
					break;
				}
				}
			} , 1
			);
#else // !NEW_THREADS
#pragma omp parallel sections
			{
#pragma omp section
				{
					double t = Time();
					if( block<blockNum ) FetchInput( block+1 );
				}
#pragma omp section
				{
					double t = Time();
					if( block>0 ) SetOutput( block-1 );
				}
#pragma omp section
				{
					double t = Time();
					RGBPixel *_inRows = inRows[block&1] , *_outRows = outRows[block&1];
					int rEnd = rStart + ROW_BLOCK_SIZE < h ? rStart + ROW_BLOCK_SIZE : h;

					// Expand the next block of values
					begin[0] = 0 , begin[1] = rStart , end[0] = w , end[1] = rEnd;
#ifdef NEW_THREADS
					Pointer( Point< Real , Colors > ) outBlock = tree.template regularGridUpSample< true >( tp , solution , begin , end );
#else // !NEW_THREADS
					Pointer( Point< Real , Colors > ) outBlock = tree.template regularGridUpSample< true >( solution , begin , end );
#endif // NEW_THREADS
					int size = (rEnd-rStart)*w;
#pragma omp parallel for
					for( int ii=0 ; ii<size ; ii++ )
					{
						Point< Real , Colors > c = Point< Real , Colors >( _inRows[ii][0] , _inRows[ii][1] , _inRows[ii][2] ) / 255;
						c += outBlock[ii] - average;
						_outRows[ii] = RGBPixel( c[0] , c[1] , c[2] );
					}
					DeletePointer( outBlock );
				}
			}
#endif // NEW_THREADS
		}
		// Write out the last block
		SetOutput( blockNum-1 );
		if( Verbose.set ) printf( "Wrote output: " ) , p.print( true );
		delete[] inRows[0];
		delete[] outRows[0];
		delete[] inRows[1];
		delete[] outRows[1];
		if( inRow ) delete[] inRow;
		delete in;
		delete out;
	}
}

#ifdef FAST_COMPILE
#else // !FAST_COMPILE
template< typename Real >
void _Execute( void )
{
	switch( Degree.value )
	{
	case 1: _Execute< Real , 1 >() ; break;
	case 2: _Execute< Real , 2 >() ; break;
//	case 3: _Execute< Real , 3 >() ; break;
//	case 4: _Execute< Real , 4 >() ; break;
	default: ERROR_OUT( "Only B-Splines of degree 1 - 2 are supported" );
	}
}
#endif // FAST_COMPILE

int main( int argc , char* argv[] )
{
	Timer timer;
	cmdLineParse( argc-1 , &argv[1] , params );
	if( MaxMemoryGB.value>0 ) SetPeakMemoryMB( MaxMemoryGB.value<<10 );
#ifdef NEW_THREADS
	ThreadPool::DefaultThreadNum = Threads.value > 1 ? Threads.value : 0;
#else // !NEW_THREADS
	omp_set_num_threads( Threads.value > 1 ? Threads.value : 1 );
#endif // NEW_THREADS
	if( Verbose.set )
	{
		printf( "*********************************************\n" );
		printf( "*********************************************\n" );
		printf( "** Running Image Stitching (Version %s) **\n" , VERSION );
		printf( "*********************************************\n" );
		printf( "*********************************************\n" );
	}

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	if( BaseDepth.value>FullDepth.value )
	{
		if( BaseDepth.set ) WARN( "Base depth must be smaller than full depth: " , BaseDepth.value , " <= " , FullDepth.value );
		BaseDepth.value = FullDepth.value;
	}

#ifdef USE_DOUBLE
	typedef double Real;
#else // !USE_DOUBLE
	typedef float  Real;
#endif // USE_DOUBLE

#ifdef FAST_COMPILE
	static const int Degree = DEFAULT_FEM_DEGREE;
#ifdef NEW_CODE
	WARN( "Compiled for degree-" , Degree , ", " , sizeof(Real)==4 ? "single" : "double" , "-precision _only_" );
#else // !NEW_CODE
	WARN( "Compiled for degree-" , Degree , ", " , sizeof(DefaultFloatType)==4 ? "single" : "double" , "-precision _only_" );
#endif // NEW_CODE
	_Execute< Real , Degree >();
#else // !FAST_COMPILE
	_Execute< Real >();
#endif // FAST_COMPILE

	if( Performance.set )
	{
		printf( "Time (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
		printf( "Peak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );
	}
	return EXIT_SUCCESS;
}
