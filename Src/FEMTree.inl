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

#include <functional>
#include <cmath>
#include <climits>
#include "MyMiscellany.h"

/////////////////////
// FEMTreeNodeData //
/////////////////////
FEMTreeNodeData::FEMTreeNodeData( void ){ flags = 0; }
FEMTreeNodeData::~FEMTreeNodeData( void ) { }


/////////////
// FEMTree //
/////////////
template< unsigned int Dim , class Real >
double FEMTree< Dim , Real >::MemoryUsage( void )
{
	double mem = double( MemoryInfo::Usage() ) / (1<<20);
	_MaxMemoryUsage = std::max< double >( mem , _MaxMemoryUsage );
	_LocalMemoryUsage = std::max< double >( mem , _LocalMemoryUsage );
	return mem;
}

#ifdef NEW_CODE
template< unsigned int Dim , class Real > FEMTree< Dim , Real >::FEMTree( size_t blockSize )
#else // !NEW_CODE
template< unsigned int Dim , class Real > FEMTree< Dim , Real >::FEMTree( int blockSize )
#endif // NEW_CODE
{
#ifdef NEW_CODE
	if( blockSize )
#else // !NEW_CODE
	if( blockSize>0 )
#endif // NEW_CODE
	{
#ifdef NEW_CODE
#ifdef NEW_THREADS
		nodeAllocators.resize( std::thread::hardware_concurrency() );
#else // !NEW_THREADS
		nodeAllocators.resize( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<nodeAllocators.size() ; i++ )
		{
			nodeAllocators[i] = new Allocator< FEMTreeNode >();
			nodeAllocators[i]->set( blockSize );
		}
#else // !NEW_CODE
		nodeAllocator = new Allocator< FEMTreeNode >();
		nodeAllocator->set( blockSize );
#endif // NEW_CODE
	}
#ifdef NEW_CODE
#else // !NEW_CODE
	else nodeAllocator = NULL;
#endif // NEW_CODE
	_nodeCount = 0;
#ifdef NEW_CODE
	_tree = FEMTreeNode::NewBrood( nodeAllocators.size() ? nodeAllocators[0] : NULL , _NodeInitializer( *this ) );
	_tree->initChildren( nodeAllocators.size() ? nodeAllocators[0] : NULL , _NodeInitializer( *this ) ) , _spaceRoot = _tree->children;
#else // !NEW_CODE
	_tree = FEMTreeNode::NewBrood( nodeAllocator , _NodeInitializer( *this ) );
	_tree->initChildren( nodeAllocator , _NodeInitializer( *this ) ) , _spaceRoot = _tree->children;
#endif // NEW_CODE
	int offset[Dim];
	for( int d=0 ; d<Dim ; d++ ) offset[d] = 0;
#ifdef NEW_CODE
	RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::ResetDepthAndOffset( _spaceRoot , 0 , offset );
#else // !NEW_CODE
	RegularTreeNode< Dim , FEMTreeNodeData >::ResetDepthAndOffset( _spaceRoot , 0 , offset );
#endif // NEW_CODE
	_depthOffset = 0;
	memset( _femSigs1 , -1 , sizeof( _femSigs1 ) );
	memset( _femSigs2 , -1 , sizeof( _femSigs2 ) );
	memset( _refinableSigs , -1 , sizeof( _refinableSigs ) );
}
template< unsigned int Dim , class Real >
#ifdef NEW_CODE
FEMTree< Dim , Real >::FEMTree( FILE* fp , size_t blockSize )
#else // !NEW_CODE
FEMTree< Dim , Real >::FEMTree( FILE* fp , int blockSize )
#endif // NEW_CODE
{
#ifdef NEW_CODE
	if( blockSize )
#else // !NEW_CODE
	if( blockSize>0 )
#endif // NEW_CODE
	{
#ifdef NEW_CODE
#ifdef NEW_THREADS
		nodeAllocators.resize( std::thread::hardware_concurrency() );
#else // !NEW_THREADS
		nodeAllocators.resize( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<nodeAllocators.size() ; i++ )
		{
			nodeAllocators[i] = new Allocator< FEMTreeNode >();
			nodeAllocators[i]->set( blockSize );
		}
#else // !NEW_CODE
		nodeAllocator = new Allocator< FEMTreeNode >();
		nodeAllocator->set( blockSize );
#endif // NEW_CODE
	}
#ifdef NEW_CODE
	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
#else // !NEW_CODE
	else nodeAllocator = NULL;
#endif // NEW_CODE
	if( fp )
	{
		if( fread( &_depthOffset , sizeof( int ) , 1 , fp )!=1 ) ERROR_OUT( "Failed to read depth offset" );
		_tree = FEMTreeNode::NewBrood( nodeAllocator , _NodeInitializer( *this ) );
		_tree->read( fp , nodeAllocator , _NodeInitializer( *this ) );
		_maxDepth = _tree->maxDepth() - _depthOffset;

		_spaceRoot = _tree->children;

		if( _depthOffset>1 )
		{
			_spaceRoot = _tree->children + (1<<Dim)-1;
			for( int d=1 ; d<_depthOffset ; d++ )
				if( !_spaceRoot->children ) ERROR_OUT( "Expected children" );
				else _spaceRoot = _spaceRoot->children;
		}
		_sNodes.set( *_tree , NULL );
	}
	else
	{
		_tree = FEMTreeNode::NewBrood( nodeAllocator , _NodeInitializer( *this ) );
		_tree->initChildren( nodeAllocator , _NodeInitializer( *this ) ) , _spaceRoot = _tree->children;
		int offset[Dim];
		for( int d=0 ; d<Dim ; d++ ) offset[d] = 0;
#ifdef NEW_CODE
		RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::ResetDepthAndOffset( _spaceRoot , 0 , offset );
#else // !NEW_CODE
		RegularTreeNode< Dim , FEMTreeNodeData >::ResetDepthAndOffset( _spaceRoot , 0 , offset );
#endif // NEW_CODE
		_depthOffset = 0;
	}
}
template< unsigned int Dim , class Real > void FEMTree< Dim , Real >::write( FILE* fp ) const
{
	fwrite( &_depthOffset , sizeof( int ) , 1 , fp );
	_tree->write( fp );
}

template< unsigned int Dim , class Real >
#ifdef NEW_CODE
const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* FEMTree< Dim , Real >::leaf( Point< Real , Dim > p ) const
#else // !NEW_CODE
const RegularTreeNode< Dim , FEMTreeNodeData >* FEMTree< Dim , Real >::leaf( Point< Real , Dim > p ) const
#endif // NEW_CODE
{
	if( !_InBounds( p ) ) return NULL;
	Point< Real , Dim > center;
	for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
	Real width = Real(1.0);
	FEMTreeNode* node = _spaceRoot;
	while( node->children )
	{
		int cIndex = FEMTreeNode::ChildIndex( center , p );
		node = node->children + cIndex;
		width /= 2;
		for( int d=0 ; d<Dim ; d++ )
			if( (cIndex>>d) & 1 ) center[d] += width/2;
			else                  center[d] -= width/2;
	}
	return node;
}
template< unsigned int Dim , class Real >
#ifdef NEW_CODE
RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* FEMTree< Dim , Real >::_leaf( Allocator< FEMTreeNode > *nodeAllocator , Point< Real , Dim > p , LocalDepth maxDepth )
#else // !NEW_CODE
RegularTreeNode< Dim , FEMTreeNodeData >* FEMTree< Dim , Real >::leaf( Point< Real , Dim > p , LocalDepth maxDepth )
#endif // NEW_CODE
{
	if( !_InBounds( p ) ) return NULL;
	Point< Real , Dim > center;
	for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
	Real width = Real(1.0);
	FEMTreeNode* node = _spaceRoot;
	LocalDepth d = _localDepth( node );
	while( ( d<0 && node->children ) || ( d>=0 && d<maxDepth ) )
	{
		if( !node->children ) node->initChildren( nodeAllocator , _NodeInitializer( *this ) );
		int cIndex = FEMTreeNode::ChildIndex( center , p );
		node = node->children + cIndex;
		d++;
		width /= 2;
		for( int d=0 ; d<Dim ; d++ )
			if( (cIndex>>d) & 1 ) center[d] += width/2;
			else                  center[d] -= width/2;
	}
	return node;
}
#ifdef NEW_CODE
template< unsigned int Dim , class Real >
RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* FEMTree< Dim , Real >::_leaf_s( Allocator< FEMTreeNode > *nodeAllocator , Point< Real , Dim > p , LocalDepth maxDepth )
{
	if( !_InBounds( p ) ) return NULL;
	Point< Real , Dim > center;
	for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
	Real width = Real(1.0);
	FEMTreeNode* node = _spaceRoot;
	LocalDepth d = _localDepth( node );
	while( ( d<0 && node->children ) || ( d>=0 && d<maxDepth ) )
	{
		if( !node->children ) node->initChildren_s( nodeAllocator , _NodeInitializer( *this ) );
		int cIndex = FEMTreeNode::ChildIndex( center , p );
		node = node->children + cIndex;
		d++;
		width /= 2;
		for( int d=0 ; d<Dim ; d++ )
			if( (cIndex>>d) & 1 ) center[d] += width/2;
			else                  center[d] -= width/2;
	}
	return node;
}
#endif // NEW_CODE

template< unsigned int Dim , class Real > bool FEMTree< Dim , Real >::_InBounds( Point< Real , Dim > p ){ for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return false ; return true; }
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSignatures >
bool FEMTree< Dim , Real >::isValidFEMNode( UIntPack< FEMSignatures ... > , const FEMTreeNode* node ) const
{
	if( GetGhostFlag< Dim >( node ) ) return false;
	LocalDepth d ; LocalOffset off ; _localDepthAndOffset( node , d , off );
	if( d<0 ) return false;
	return FEMIntegrator::IsValidFEMNode( UIntPack< FEMSignatures ... >() , d , off );
}
template< unsigned int Dim , class Real >
bool FEMTree< Dim , Real >::isValidSpaceNode( const FEMTreeNode* node ) const
{
	if( !node ) return false;
	LocalDepth d ; LocalOffset off ; _localDepthAndOffset( node , d , off );
	if( d<0 ) return false;
	int res = 1<<d;
	for( int dd=0 ; dd<Dim ; dd++ ) if( off[dd]<0 || off[dd]>=res ) return false;
	return true;
}

template< unsigned int Dim , class Real >
template< unsigned int ... Degrees >
#ifdef NEW_CODE
void FEMTree< Dim , Real >::_setFullDepth( UIntPack< Degrees ... > , Allocator< FEMTreeNode > *nodeAllocator , FEMTreeNode* node , LocalDepth depth )
#else // !NEW_CODE
void FEMTree< Dim , Real >::_setFullDepth( UIntPack< Degrees ... > , FEMTreeNode* node , LocalDepth depth )
#endif // NEW_CODE
{
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	bool refine = d<depth && ( d<0 || !FEMIntegrator::IsOutOfBounds( UIntPack< FEMDegreeAndBType< Degrees , BOUNDARY_FREE >::Signature ... >() , d , off ) );
	if( refine )
	{
		if( !node->children ) node->initChildren( nodeAllocator , _NodeInitializer( *this ) );
#ifdef NEW_CODE
		for( int c=0 ; c<(1<<Dim) ; c++ ) _setFullDepth( UIntPack< Degrees ... >() , nodeAllocator , node->children+c , depth );
#else // !NEW_CODE
		for( int c=0 ; c<(1<<Dim) ; c++ ) _setFullDepth( UIntPack< Degrees ... >() , node->children+c , depth );
#endif // NEW_CODE
	}
}
template< unsigned int Dim , class Real >
template< unsigned int ... Degrees >
#ifdef NEW_CODE
void FEMTree< Dim , Real >::_setFullDepth( UIntPack< Degrees ... > , Allocator< FEMTreeNode > *nodeAllocator , LocalDepth depth )
#else // !NEW_CODE
void FEMTree< Dim , Real >::_setFullDepth( UIntPack< Degrees ... > , LocalDepth depth )
#endif // NEW_CODE
{
	if( !_tree->children ) _tree->initChildren( nodeAllocator , _NodeInitializer( *this ) );
#ifdef NEW_CODE
	for( int c=0 ; c<(1<<Dim) ; c++ ) _setFullDepth( UIntPack< Degrees ... >() , nodeAllocator , _tree->children+c , depth );
#else // !NEW_CODE
	for( int c=0 ; c<(1<<Dim) ; c++ ) _setFullDepth( UIntPack< Degrees ... >() , _tree->children+c , depth );
#endif // NEW_CODE
}
template< unsigned int Dim , class Real >
template< unsigned int ... Degrees >
typename FEMTree< Dim , Real >::LocalDepth FEMTree< Dim , Real >::_getFullDepth( UIntPack< Degrees ... > , const FEMTreeNode* node ) const
{
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	bool refine = d<0 || !FEMIntegrator::IsOutOfBounds( UIntPack< FEMDegreeAndBType< Degrees , BOUNDARY_FREE >::Signature ... >() , d , off );

	if( refine )
	{
		if( !node->children ) return d;
		else
		{
			LocalDepth depth = INT_MAX;
			for( int c=0 ; c<(1<<Dim) ; c++ )
			{
				LocalDepth d = _getFullDepth( UIntPack< Degrees ... >() , node->children+c );
				if( d<depth ) depth = d;
			}
			return depth;
		}
	}
	else return INT_MAX;
}
template< unsigned int Dim , class Real >
template< unsigned int ... Degrees >
typename FEMTree< Dim , Real >::LocalDepth FEMTree< Dim , Real >::getFullDepth( UIntPack< Degrees ... > ) const
{
	if( !_tree->children ) return -1;
	LocalDepth depth = INT_MAX;
	for( int c=0 ; c<(1<<Dim) ; c++ )
	{
		LocalDepth d = _getFullDepth( UIntPack< Degrees ... >() , _tree->children+c );
		if( d<depth ) depth = d;
	}
	return depth;
}

template< unsigned int Dim , class Real >
template< unsigned int LeftRadius , unsigned int RightRadius , class ... DenseOrSparseNodeData > 
void FEMTree< Dim , Real >::thicken( FEMTreeNode **nodes , size_t nodeCount, DenseOrSparseNodeData* ... data )
{
#ifdef NEW_CODE
	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
#endif // NEW_CODE
#ifdef NEW_CODE
	std::vector< node_index_type > map( _nodeCount );
	for( node_index_type i=0 ; i<_nodeCount ; i++ ) map[i] = i;
#else // !NEW_CODE
	std::vector< int > map( _nodeCount );
	for( int i=0 ; i<_nodeCount ; i++ ) map[i] = i;
#endif // NEW_CODE
	{
		int d=0 , off[Dim];
		for( int d=0 ; d<Dim ; d++ ) off[d] = 0;
		FEMTreeNode::ResetDepthAndOffset( _tree , d , off );
	}
#ifdef NEW_CODE
	typename RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::template NeighborKey< IsotropicUIntPack< Dim , LeftRadius > , IsotropicUIntPack< Dim , RightRadius > > neighborKey;
#else // !NEW_CODE
	typename RegularTreeNode< Dim , FEMTreeNodeData >::template NeighborKey< IsotropicUIntPack< Dim , LeftRadius > , IsotropicUIntPack< Dim , RightRadius > > neighborKey;
#endif // NEW_CODE
	neighborKey.set( _tree->maxDepth() );
#ifdef NEW_CODE
	for( size_t i=0 ; i<nodeCount ; i++ ) neighborKey.template getNeighbors< true >( nodes[i] , nodeAllocator , _NodeInitializer( *this ) );
#else // !NEW_CODE
	for( int i=0 ; i<nodeCount ; i++ ) neighborKey.template getNeighbors< true >( nodes[i] , nodeAllocator , _NodeInitializer( *this ) );
#endif // NEW_CODE
	{
		int d=0 , off[Dim];
		for( int d=0 ; d<Dim ; d++ ) off[d] = 0;
		FEMTreeNode::ResetDepthAndOffset( _spaceRoot , d , off );
	}

	_reorderDenseOrSparseNodeData( &map[0] , _nodeCount , data ... );
}
template< unsigned int Dim , class Real >
template< unsigned int LeftRadius , unsigned int RightRadius , class IsThickenNode , class ... DenseOrSparseNodeData > 
void FEMTree< Dim , Real >::thicken( IsThickenNode F , DenseOrSparseNodeData* ... data )
{
	std::vector< FEMTreeNode* > nodes;
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( IsActiveNode( node ) && F( node ) ) nodes.push_back( node );
	thicken< LeftRadius , RightRadius >( &nodes[0] , nodes.size() , data ... );
}

template< unsigned int Dim , class Real >
template< unsigned int DensityDegree >
#ifdef NEW_THREADS
typename FEMTree< Dim , Real >::template DensityEstimator< DensityDegree >* FEMTree< Dim , Real >::setDensityEstimator( ThreadPool &tp , const std::vector< PointSample >& samples , LocalDepth splatDepth , Real samplesPerNode , int coDimension )
#else // !NEW_THREADS
typename FEMTree< Dim , Real >::template DensityEstimator< DensityDegree >* FEMTree< Dim , Real >::setDensityEstimator( const std::vector< PointSample >& samples , LocalDepth splatDepth , Real samplesPerNode , int coDimension )
#endif // NEW_THREADS
{
#ifdef NEW_CODE
	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
#endif // NEW_CODE
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	splatDepth = std::max< LocalDepth >( 0 , std::min< LocalDepth >( splatDepth , maxDepth ) );
	DensityEstimator< DensityDegree >* _density = new DensityEstimator< DensityDegree >( splatDepth , coDimension );
	DensityEstimator< DensityDegree >& density = *_density;
	PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > densityKey;
	densityKey.set( _localToGlobal( splatDepth ) );

#ifdef NEW_CODE
	std::vector< node_index_type > sampleMap( nodeCount() , -1 );
#else // !NEW_CODE
	std::vector< int > sampleMap( nodeCount() , -1 );
#endif // NEW_CODE
#ifdef NEW_THREADS
	tp.parallel_for( 0 , samples.size() , [&]( unsigned int , size_t i ){ if( samples[i].sample.weight>0 ) sampleMap[ samples[i].node->nodeData.nodeIndex ] = (node_index_type)i; } );
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ ) if( samples[i].sample.weight>0 ) sampleMap[ samples[i].node->nodeData.nodeIndex ] = i;
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ ) if( samples[i].sample.weight>0 ) sampleMap[ samples[i].node->nodeData.nodeIndex ] = i;
#endif // NEW_CODE
#endif // NEW_THREADS
	std::function< ProjectiveData< Point< Real , Dim > , Real > ( FEMTreeNode* ) > SetDensity = [&] ( FEMTreeNode* node )
	{
		ProjectiveData< Point< Real , Dim > , Real > sample;
		LocalDepth d = _localDepth( node );
#ifdef NEW_CODE
		node_index_type idx = node->nodeData.nodeIndex;
#else // !NEW_CODE
		int idx = node->nodeData.nodeIndex;
#endif // NEW_CODE
		if( node->children )
			for( int c=0 ; c<(1<<Dim) ; c++ )
			{
				ProjectiveData< Point< Real , Dim > , Real > s = SetDensity( node->children + c );
				if( d<=splatDepth && s.weight>0 )
				{
					Point< Real , Dim > p = s.data / s.weight;
					Real w = s.weight / samplesPerNode;
#ifdef NEW_CODE
					_addWeightContribution( nodeAllocator , density , node , p , densityKey , w );
#else // !NEW_CODE
					_addWeightContribution( density , node , p , densityKey , w );
#endif // NEW_CODE
				}
				sample += s;
			}
#ifdef NEW_CODE
		else if( idx<(node_index_type)sampleMap.size() && sampleMap[idx]!=-1 )
#else // !NEW_CODE
		else if( idx<sampleMap.size() && sampleMap[idx]!=-1 )
#endif // NEW_CODE
		{
			sample = samples[ sampleMap[ idx ] ].sample;
			if( d<=splatDepth && sample.weight>0 )
			{
				Point< Real , Dim > p = sample.data / sample.weight;
				Real w = sample.weight / samplesPerNode;
#ifdef NEW_CODE
				_addWeightContribution( nodeAllocator , density , node , p , densityKey , w );
#else // !NEW_CODE
				_addWeightContribution( density , node , p , densityKey , w );
#endif // NEW_CODE
			}
		}
		return sample;
	};
	SetDensity( _spaceRoot );

	MemoryUsage();
	return _density;
}
template< unsigned int Dim , class Real >
template< unsigned int ... NormalSigs , unsigned int DensityDegree , class Data >
#ifdef NEW_THREADS
SparseNodeData< Point< Real , Dim > , UIntPack< NormalSigs ... > > FEMTree< Dim , Real >::setNormalField( UIntPack< NormalSigs ... > , ThreadPool &tp , const std::vector< PointSample >& samples , const std::vector< Data >& normalData , const DensityEstimator< DensityDegree >* density , Real& pointWeightSum , std::function< Real ( Real ) > BiasFunction )
#else // !NEW_THREADS
SparseNodeData< Point< Real , Dim > , UIntPack< NormalSigs ... > > FEMTree< Dim , Real >::setNormalField( UIntPack< NormalSigs ... > , const std::vector< PointSample >& samples , const std::vector< Data >& normalData , const DensityEstimator< DensityDegree >* density , Real& pointWeightSum , std::function< Real ( Real ) > BiasFunction )
#endif // NEW_THREADS
{
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	typedef PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > DensityKey;
	typedef UIntPack< FEMSignature< NormalSigs >::Degree ... > NormalDegrees;
	typedef PointSupportKey< UIntPack< FEMSignature< NormalSigs >::Degree ... > > NormalKey;
#ifdef NEW_THREADS
	std::vector< DensityKey > densityKeys( tp.threadNum() );
	std::vector<  NormalKey >  normalKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< DensityKey > densityKeys( omp_get_max_threads() );
	std::vector<  NormalKey >  normalKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	bool oneKey = DensityDegree==NormalDegrees::Min() && DensityDegree==NormalDegrees::Max();
#ifdef NEW_CODE
	for( size_t i=0 ; i<densityKeys.size() ; i++ ) densityKeys[i].set( _localToGlobal( maxDepth ) );
	if( !oneKey ) for( size_t i=0 ; i<normalKeys.size() ; i++ ) normalKeys[i].set( _localToGlobal( maxDepth ) );
#else // !NEW_CODE
	for( int i=0 ; i<densityKeys.size() ; i++ ) densityKeys[i].set( _localToGlobal( maxDepth ) );
	if( !oneKey ) for( int i=0 ; i<normalKeys.size() ; i++ ) normalKeys[i].set( _localToGlobal( maxDepth ) );
#endif // NEW_CODE
	Real weightSum = 0;
	pointWeightSum = 0;
	SparseNodeData< Point< Real , Dim > , UIntPack< NormalSigs ... > > normalField;
	Real _pointWeightSum = 0;
#ifdef NEW_THREADS
	tp.parallel_for( 0 , samples.size() , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : weightSum , _pointWeightSum )
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
#ifdef NEW_THREADS
		DensityKey& densityKey = densityKeys[ thread ];
		NormalKey& normalKey = normalKeys[ thread ];
#else // !NEW_THREADS
		DensityKey& densityKey = densityKeys[ omp_get_thread_num() ];
		NormalKey& normalKey = normalKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
		const ProjectiveData< Point< Real , Dim > , Real >& sample = samples[i].sample;
		if( sample.weight>0 )
		{
			Point< Real , Dim > p = sample.data / sample.weight , n = std::get< 0 >( normalData[i].data ).data;
			Real l = (Real)Length( n );
			// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
#ifdef NEW_THREADS
			if( !l ) return;
#else // !NEW_THREADS
			if( !l ) continue;
#endif // NEW_THREADS
			Real confidence = l / sample.weight;
			n *= sample.weight / l;
			Real depthBias = BiasFunction( confidence );
#ifdef NEW_THREADS
			AddAtomic( weightSum , sample.weight );
#else // !NEW_THREADS
			weightSum += sample.weight;
#endif // NEW_THREADS
			if( !_InBounds(p) )
			{
				WARN( "Point sample is out of bounds" );
#ifdef NEW_THREADS
				return;
#else // !NEW_THREADS
				continue;
#endif // NEW_THREADS
			}
#ifdef NEW_CODE
#ifdef NEW_THREADS
			Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[ thread ] : NULL;
#else // !NEW_THREADS
			Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[ omp_get_thread_num() ] : NULL;
#endif // NEW_THREADS
#endif // NEW_CODE

#ifdef NEW_CODE
#ifdef NEW_THREADS
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
			if( density ) AddAtomic( _pointWeightSum , _splatPointData< true , DensityDegree , Point< Real , Dim > >( nodeAllocator , *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight );
#else // !__GNUC__ || __GNUC__ >=5
			if( density ) AddAtomic( _pointWeightSum , _splatPointData< true , DensityDegree , Point< Real , Dim > , NormalSigs ... >( nodeAllocator , *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight );
#endif // __GNUC__ || __GNUC__ < 4
			else
			{
				Real width = (Real)( 1.0 / ( 1<<maxDepth ) );
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
				_splatPointData< true , Point< Real , Dim > >( nodeAllocator , _leaf_s( nodeAllocator , p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#else // !__GNUC__ || __GNUC__ >=5
				_splatPointData< true , Point< Real , Dim > , NormalSigs ... >( nodeAllocator , _leaf_s( nodeAllocator , p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#endif // __GNUC__ || __GNUC__ < 4
				AddAtomic( _pointWeightSum , sample.weight );
			}
#else // !NEW_THREADS
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
			if( density ) _pointWeightSum += _splatPointData< true , DensityDegree , Point< Real , Dim > >( nodeAllocator , *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight;
#else // !__GNUC__ || __GNUC__ >=5
			if( density ) _pointWeightSum += _splatPointData< true , DensityDegree , Point< Real , Dim > , NormalSigs ... >( nodeAllocator , *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight;
#endif // __GNUC__ || __GNUC__ < 4
			else
			{
				Real width = (Real)( 1.0 / ( 1<<maxDepth ) );
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
				_splatPointData< true , Point< Real , Dim > >( nodeAllocator , _leaf( nodeAllocator , p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#else // !__GNUC__ || __GNUC__ >=5
				_splatPointData< true , Point< Real , Dim > , NormalSigs ... >( nodeAllocator , _leaf( nodeAllocator , p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#endif // __GNUC__ || __GNUC__ < 4
				_pointWeightSum += sample.weight;
			}
#endif // NEW_THREADS
#else // !NEW_CODE
#ifdef NEW_THREADS
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
			if( density ) AddAtomic( _pointWeightSum , _splatPointData< true , DensityDegree , Point< Real , Dim > >( *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight );
#else // !__GNUC__ || __GNUC__ >=5
			if( density ) AddAtomic( _pointWeightSum , _splatPointData< true , DensityDegree , Point< Real , Dim > , NormalSigs ... >( *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight );
#endif // __GNUC__ || __GNUC__ < 4
			else
			{
				Real width = (Real)( 1.0 / ( 1<<maxDepth ) );
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
				_splatPointData< true , Point< Real , Dim > >( leaf( p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#else // !__GNUC__ || __GNUC__ >=5
				_splatPointData< true , Point< Real , Dim > , NormalSigs ... >( leaf( p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#endif // __GNUC__ || __GNUC__ < 4
				AddAtomic( _pointWeightSum , sample.weight );
			}
#else // !NEW_THREADS
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
			if( density ) _pointWeightSum += _splatPointData< true , DensityDegree , Point< Real , Dim > >( *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight;
#else // !__GNUC__ || __GNUC__ >=5
			if( density ) _pointWeightSum += _splatPointData< true , DensityDegree , Point< Real , Dim > , NormalSigs ... >( *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight;
#endif // __GNUC__ || __GNUC__ < 4
			else
			{
				Real width = (Real)( 1.0 / ( 1<<maxDepth ) );
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
				_splatPointData< true , Point< Real , Dim > >( leaf( p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#else // !__GNUC__ || __GNUC__ >=5
				_splatPointData< true , Point< Real , Dim > , NormalSigs ... >( leaf( p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#endif // __GNUC__ || __GNUC__ < 4
				_pointWeightSum += sample.weight;
			}
#endif // NEW_THREADS
#endif // NEW_CODE
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	pointWeightSum = _pointWeightSum / weightSum;
	MemoryUsage();
	return normalField;
}
template< unsigned int Dim , class Real >
template< unsigned int DataSig , bool CreateNodes , unsigned int DensityDegree , class Data >
SparseNodeData< Data , IsotropicUIntPack< Dim , DataSig > > FEMTree< Dim , Real >::setSingleDepthDataField( const std::vector< PointSample >& samples , const std::vector< Data >& sampleData , const DensityEstimator< DensityDegree >* density )
{
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > densityKey;
	PointSupportKey< IsotropicUIntPack< Dim , FEMSignature< DataSig >::Degree > > dataKey;
	densityKey.set( _localToGlobal( maxDepth ) ) , dataKey.set( _localToGlobal( maxDepth ) );

	SparseNodeData< Data , IsotropicUIntPack< Dim , DataSig > > dataField;
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const ProjectiveData< Point< Real , Dim > , Real >& sample = samples[i].sample;
		const Data& data = sampleData[i];
		Point< Real , Dim > p = sample.weight==0 ? sample.data : sample.data / sample.weight;
		if( !_InBounds(p) )
		{
			WARN( "Point is out of bounds" );
			continue;
		}
		if( density ) _splatPointData< CreateNodes , DensityDegree , DataSig >( *density             , p , data * sample.weight , dataField , densityKey , dataKey , 0 , maxDepth , Dim );
		else          _splatPointData< CreateNodes ,                 DataSig >( leaf( p , maxDepth ) , p , data * sample.weight , dataField , dataKey );
	}
	MemoryUsage();
	return dataField;
}
template< unsigned int Dim , class Real >
template< unsigned int DataSig , bool CreateNodes , unsigned int DensityDegree , class Data >
SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > FEMTree< Dim , Real >::setDataField( const std::vector< PointSample >& samples , std::vector< Data >& sampleData , const DensityEstimator< DensityDegree >* density , bool nearest )
{
#ifdef NEW_CODE
	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
#endif // NEW_CODE
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > densityKey;
	PointSupportKey< IsotropicUIntPack< Dim , FEMSignature< DataSig >::Degree > > dataKey;
	densityKey.set( _localToGlobal( maxDepth ) ) , dataKey.set( _localToGlobal( maxDepth ) );

	SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > dataField;
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const ProjectiveData< Point< Real , Dim > , Real >& sample = samples[i].sample;
		const Data& data = sampleData[i];
		Point< Real , Dim > p = sample.weight==0 ? sample.data : sample.data / sample.weight;
		if( !_InBounds(p) )
		{
			WARN( "Point is out of bounds" );
			continue;
		}
		if( nearest ) _nearestMultiSplatPointData< DensityDegree >( density , (FEMTreeNode*)samples[i].node , p , ProjectiveData< Data , Real >( data , sample.weight ) , dataField , densityKey , 2 );
#ifdef NEW_CODE
		else          _multiSplatPointData< CreateNodes , DensityDegree >( nodeAllocator , density , (FEMTreeNode*)samples[i].node , p , ProjectiveData< Data , Real >( data , sample.weight ) , dataField , densityKey , dataKey , 2 );
#else // !NEW_CODE
		else          _multiSplatPointData< CreateNodes , DensityDegree >( density , (FEMTreeNode*)samples[i].node , p , ProjectiveData< Data , Real >( data , sample.weight ) , dataField , densityKey , dataKey , 2 );
#endif // NEW_CODE
	}
	MemoryUsage();
	return dataField;
}
template< unsigned int Dim , class Real >
template< unsigned int MaxDegree , class HasDataFunctor , class ... DenseOrSparseNodeData >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::finalizeForMultigrid( ThreadPool &tp , LocalDepth fullDepth , const HasDataFunctor F , DenseOrSparseNodeData* ... data )
#else // !NEW_THREADS
void FEMTree< Dim , Real >::finalizeForMultigrid( LocalDepth fullDepth , const HasDataFunctor F , DenseOrSparseNodeData* ... data )
#endif // NEW_THREADS
{
#ifdef NEW_CODE
	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
#endif // NEW_CODE
	_depthOffset = 1;
	while( _localInset( 0 ) + BSplineEvaluationData< FEMDegreeAndBType< MaxDegree >::Signature >::Begin( 0 )<0 || _localInset( 0 ) + BSplineEvaluationData< FEMDegreeAndBType< MaxDegree >::Signature >::End( 0 )>(1<<_depthOffset) )
	{
		//                       +-+-+-+-+-+-+-+-+
		//                       | | | | | | | | |
		//                       +-+-+-+-+-+-+-+-+
		//                       | | | | | | | | |
		//          +-+-+-+-+    +-+-+-+-+-+-+-+-+
		//          | | | | |    | | | | | | | | |
		// +-+-+    +-+-+-+-+    +-+-+-+-+-+-+-+-+
		// |*| |    | | | | |    | | | | | | | | |
		// +-o-+ -> +-+-o-+-+ -> +-+-+-+-o-+-+-+-+
		// | | |    | | |*| |    | | | | |*| | | |
		// +-+-+    +-+-+-+-+    +-+-+-+-+-+-+-+-+
		//          | | | | |    | | | | | | | | |
		//          +-+-+-+-+    +-+-+-+-+-+-+-+-+
		//                       | | | | | | | | |
		//                       +-+-+-+-+-+-+-+-+
		//                       | | | | | | | | |
		//                       +-+-+-+-+-+-+-+-+

		FEMTreeNode* newSpaceRootParent = FEMTreeNode::NewBrood( nodeAllocator , _NodeInitializer( *this ) );
		FEMTreeNode* oldSpaceRootParent = _spaceRoot->parent;
		int corner = _depthOffset<=1 ? (1<<Dim)-1 : 0;
		newSpaceRootParent[corner].children = _spaceRoot;
		oldSpaceRootParent->children = newSpaceRootParent;
		for( int c=0 ; c<(1<<Dim) ; c++ ) _spaceRoot[c].parent = newSpaceRootParent + corner , newSpaceRootParent[c].parent = oldSpaceRootParent;
		_depthOffset++;
	}
	int d=0 , off[Dim];
	for( int d=0 ; d<Dim ; d++ ) off[d] = 0;
	FEMTreeNode::ResetDepthAndOffset( _tree , d , off );
	_maxDepth = _spaceRoot->maxDepth();
	// Make the low-resolution part of the tree be complete
	fullDepth = std::max< LocalDepth >( 0 , std::min< LocalDepth >( _maxDepth , fullDepth ) );
#ifdef NEW_CODE
	_setFullDepth( IsotropicUIntPack< Dim , MaxDegree >() , nodeAllocator , fullDepth );
#else // !NEW_CODE
	_setFullDepth( IsotropicUIntPack< Dim , MaxDegree >() , fullDepth );
#endif // NEW_CODE
	// Clear all the flags and make everything that is not low-res a ghost node
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) node->nodeData.flags = 0 , SetGhostFlag< Dim >( node , _localDepth( node )>fullDepth );

	// Set the ghost nodes for the high-res part of the tree
#ifdef NEW_THREADS
	_clipTree( tp , F , fullDepth );
#else // !NEW_THREADS
	_clipTree( F , fullDepth );
#endif // NEW_THREADS

	const int OverlapRadius = -BSplineOverlapSizes< MaxDegree , MaxDegree >::OverlapStart;
	int maxDepth = _tree->maxDepth( );
	typedef typename FEMTreeNode::template NeighborKey< IsotropicUIntPack< Dim , OverlapRadius > , IsotropicUIntPack< Dim , OverlapRadius > > NeighborKey;

#ifdef NEW_THREADS
	std::vector< NeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< NeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( _maxDepth-1 ) );

	for( LocalDepth d=_maxDepth-1 ; d>=0 ; d-- )
	{
		std::vector< FEMTreeNode* > nodes;
		auto NodeTerminationLambda = [&]( const FEMTreeNode *node ){ return _localDepth( node )==d; };
		for( FEMTreeNode* node=_tree->nextNode( NodeTerminationLambda , NULL ) ; node ; node=_tree->nextNode( NodeTerminationLambda , node ) ) if( _localDepth( node )==d && IsActiveNode< Dim >( node->children ) ) nodes.push_back( node );
#ifdef NEW_THREADS
		tp.parallel_for( 0 , nodes.size() , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=0 ; i<(node_index_type)nodes.size() ; i++ )
#else // !NEW_CODE
		for( int i=0 ; i<nodes.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
#ifdef NEW_THREADS
			NeighborKey& neighborKey = neighborKeys[ thread ];
#else // !NEW_THREADS
			NeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			FEMTreeNode* node = nodes[i];
#ifdef NEW_CODE
#ifdef NEW_THREADS
			neighborKey.template getNeighbors< true >( node , nodeAllocators.size() ? nodeAllocators[ thread ] : NULL , _NodeInitializer( *this ) );
#else // !NEW_THREADS
			neighborKey.template getNeighbors< true >( node , nodeAllocators.size() ? nodeAllocators[ omp_get_thread_num() ] : NULL , _NodeInitializer( *this ) );
#endif // NEW_THREADS
#else // !NEW_CODE
			neighborKey.template getNeighbors< true >( node , nodeAllocator , _NodeInitializer( *this ) );
#endif // NEW_CODE
			Pointer( FEMTreeNode* ) nodes = neighborKey.neighbors[ _localToGlobal(d) ].neighbors().data;
			unsigned int size = neighborKey.neighbors[ _localToGlobal(d) ].neighbors.Size;
			for( unsigned int i=0 ; i<size ; i++ ) SetGhostFlag< Dim >( nodes[i] , false );
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
#ifdef NEW_CODE
	std::vector< node_index_type > map;
#else // !NEW_CODE
	std::vector< int > map;
#endif // NEW_CODE
	_sNodes.set( *_tree , &map );
#ifdef NEW_THREADS
	_setSpaceValidityFlags( tp );
#else // !NEW_THREADS
	_setSpaceValidityFlags();
#endif // NEW_THREADS
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( !IsActiveNode< Dim >( node ) ) node->nodeData.nodeIndex = -1;
	_reorderDenseOrSparseNodeData( &map[0] , _sNodes.size() , data ... );
	MemoryUsage();
}

template< unsigned int Dim , class Real >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_setSpaceValidityFlags( ThreadPool &tp ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_setSpaceValidityFlags( void ) const
#endif // NEW_THREADS
{
#ifdef NEW_THREADS
	tp.parallel_for( 0 , _sNodes.size() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)_sNodes.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<_sNodes.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		const unsigned char MASK = ~( FEMTreeNodeData::SPACE_FLAG );
		_sNodes.treeNodes[i]->nodeData.flags &= MASK;
		if( isValidSpaceNode( _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::SPACE_FLAG;
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs1 >
void FEMTree< Dim , Real >::_setFEM1ValidityFlags( UIntPack< FEMSigs1 ... > ) const
{
	bool needToReset;
	unsigned int femSigs1[] = { FEMSigs1 ... };
#ifdef NEW_THREADS
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock( m );
		needToReset = memcmp( femSigs1 , _femSigs1 , sizeof( _femSigs1 ) )!=0;
		if( needToReset ) memcpy( _femSigs1 , femSigs1 , sizeof( _femSigs1 ) );
	}
#else // !NEW_THREADS
#pragma omp critical (set_fem_1_validity_flags)
	{
		needToReset = memcmp( femSigs1 , _femSigs1 , sizeof( _femSigs1 ) )!=0;
		if( needToReset ) memcpy( _femSigs1 , femSigs1 , sizeof( _femSigs1 ) );
	}
#endif // NEW_THREADS
	if( needToReset )
#ifdef NEW_CODE
		for( node_index_type i=0 ; i<(node_index_type)_sNodes.size() ; i++ )
#else // !NEW_CODE
		for( int i=0 ; i<_sNodes.size() ; i++ )
#endif // NEW_CODE
		{
			const unsigned char MASK = ~( FEMTreeNodeData::FEM_FLAG_1 );
			_sNodes.treeNodes[i]->nodeData.flags &= MASK;
			if( isValidFEMNode( UIntPack< FEMSigs1 ... >() , _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::FEM_FLAG_1;
		}

}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs2 >
void FEMTree< Dim , Real >::_setFEM2ValidityFlags( UIntPack< FEMSigs2 ... > ) const
{
	bool needToReset;
	unsigned int femSigs2[] = { FEMSigs2 ... };
#ifdef NEW_THREADS
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock(m);
		needToReset = memcmp( femSigs2 , _femSigs2 , sizeof( _femSigs2 ) )!=0;
		if( needToReset ) memcpy( _femSigs2 , femSigs2 , sizeof( _femSigs2 ) );
	}
#else // !NEW_THREADS
#pragma omp critical (set_fem_2_validity_flags)
	{
		needToReset = memcmp( femSigs2 , _femSigs2 , sizeof( _femSigs2 ) )!=0;
		if( needToReset ) memcpy( _femSigs2 , femSigs2 , sizeof( _femSigs2 ) );
	}
#endif // NEW_THREADS
	if( needToReset )
#ifdef NEW_CODE
		for( node_index_type i=0 ; i<(node_index_type)_sNodes.size() ; i++ )
#else // !NEW_CODE
		for( int i=0 ; i<_sNodes.size() ; i++ )
#endif // NEW_CODE
		{
			const unsigned char MASK = ~( FEMTreeNodeData::FEM_FLAG_2 );
			_sNodes.treeNodes[i]->nodeData.flags &= MASK;
			if( isValidFEMNode( UIntPack< FEMSigs2 ... >() , _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::FEM_FLAG_2;
		}
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_setRefinabilityFlags( UIntPack< FEMSigs ... > , ThreadPool &tp ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_setRefinabilityFlags( UIntPack< FEMSigs ... > ) const
#endif // NEW_THREADS
{
	bool needToReset;
	unsigned int refinableSigs[] = { FEMSigs ... };
#ifdef NEW_THREADS
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock(m);
		needToReset = memcmp( refinableSigs , _refinableSigs , sizeof( _refinableSigs ) )!=0;
		if( needToReset ) memcpy( _refinableSigs , refinableSigs , sizeof( _refinableSigs ) );
	}
#else // !NEW_THREADS
#pragma omp critical (set_refinability_flags)
	{
		needToReset = memcmp( refinableSigs , _refinableSigs , sizeof( _refinableSigs ) )!=0;
		if( needToReset ) memcpy( _refinableSigs , refinableSigs , sizeof( _refinableSigs ) );
	}
#endif // NEW_THREADS
	if( needToReset )
	{
		typedef typename FEMTreeNode::template ConstNeighborKey< UIntPack< ( - BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleStart ) ... > , UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleEnd ... > > UpSampleKey;
		typedef typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleSize ... > > UpSampleNeighbors;
		static const int UpSampleStart[] = { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleStart ... };
#ifdef NEW_THREADS
		std::vector< UpSampleKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
		std::vector< UpSampleKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( _maxDepth ) );

		for( int d=0 ; d<_maxDepth ; d++ )
#ifdef NEW_THREADS
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#else // !NEW_CODE
			for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
#ifdef NEW_THREADS
				UpSampleKey& neighborKey = neighborKeys[ thread ];
#else // !NEW_THREADS
				UpSampleKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS

				// Clear the refinability flag
				const unsigned char MASK = ~( FEMTreeNodeData::REFINABLE_FLAG );
				_sNodes.treeNodes[i]->nodeData.flags &= MASK;

				LocalDepth d ; LocalOffset pOff;
				_localDepthAndOffset( _sNodes.treeNodes[i] , d , pOff );

				// Get the supporting child neighbors
				neighborKey.getNeighbors( _sNodes.treeNodes[i] );
				UpSampleNeighbors neighbors;
				neighborKey.getChildNeighbors( 0 , _localToGlobal( d ) , neighbors );

				// Check if the child neighbors exist (i.e. that the children nodes are not ghost-nodes if they correspond to valid coefficients)
				bool refinable = true;
				LocalOffset cOff;
				WindowLoop< Dim >::Run
				(
					IsotropicUIntPack< Dim , 0 >() , UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleSize ... >() ,
					[&]( int d , int i ){ cOff[d] = pOff[d]*2 + UpSampleStart[d] + i; } ,
					[&]( const FEMTreeNode* node ){ if( GetGhostFlag< Dim >( node ) && FEMIntegrator::IsValidFEMNode( UIntPack< FEMSigs ... >() , d+1 , cOff ) ) refinable = false; } ,
					neighbors.neighbors()
				);
				if( refinable ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::REFINABLE_FLAG;
			}
#ifdef NEW_THREADS
			);
#endif // NEW_THREADS
	}
}
template< unsigned int Dim , class Real >
template< class HasDataFunctor >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_clipTree( ThreadPool &tp , const HasDataFunctor& f , LocalDepth fullDepth )
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_clipTree( const HasDataFunctor& f , LocalDepth fullDepth )
#endif // NEW_THREADS
{
	std::vector< FEMTreeNode * > nodes;
	auto NodeTerminationLambda = [&]( const FEMTreeNode *node ){ return _localDepth( node )==fullDepth; };
	for( FEMTreeNode* temp=_tree->nextNode( NodeTerminationLambda , NULL ) ; temp ; temp=_tree->nextNode( NodeTerminationLambda , temp ) ) if( _localDepth( temp )==fullDepth ) nodes.push_back( temp );
#ifdef NEW_THREADS
	tp.parallel_for( 0 , nodes.size() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)nodes.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<nodes.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		for( FEMTreeNode* node=nodes[i]->nextNode() ; node ; node=nodes[i]->nextNode(node) ) if( node->children )
		{
			bool hasData = false;
			for( int c=0 ; c<(1<<Dim) && !hasData ; c++ ) hasData |= f( node->children + c );
			for( int c=0 ; c<(1<<Dim) ; c++ ) SetGhostFlag< Dim >( node->children+c , !hasData );
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< unsigned int Dim , class Real >
template< typename T , typename Data , unsigned int PointD , typename ConstraintDual , typename SystemDual >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_ExactPointAndDataInterpolationInfo< T , Data , PointD , ConstraintDual , SystemDual >::_init( ThreadPool &tp , const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , bool noRescale )
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_ExactPointAndDataInterpolationInfo< T , Data , PointD , ConstraintDual , SystemDual >::_init( const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , bool noRescale )
#endif // NEW_THREADS
{
	_sampleSpan.resize( tree.nodesSize() );
#ifdef NEW_CODE
#ifdef NEW_THREADS
	tp.parallel_for( 0 , tree.nodesSize() , [&]( unsigned int , size_t i ){ _sampleSpan[i] = std::pair< node_index_type , node_index_type >( 0 , 0 ); } );
#else // !NEW_THREADS
#pragma omp parallel for
	for( node_index_type i=0 ; i<(node_index_type)tree.nodesSize() ; i++ ) _sampleSpan[i] = std::pair< node_index_type , node_index_type >( 0 , 0 );
#endif // NEW_THREADS
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
#else // !NEW_CODE
#pragma omp parallel for
	for( int i=0 ; i<tree.nodesSize() ; i++ ) _sampleSpan[i] = std::pair< int , int >( 0 , 0 );
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) ) _sampleSpan[ leaf->nodeData.nodeIndex ].second++;
	}
	_iData.resize( samples.size() );

#ifdef NEW_CODE
	std::function< void ( FEMTreeNode* , node_index_type & ) > SetRange = [&] ( FEMTreeNode* node , node_index_type &start )
#else // !NEW_CODE
	std::function< void ( FEMTreeNode* , int& ) > SetRange = [&] ( FEMTreeNode* node , int &start )
#endif // NEW_CODE
	{
#ifdef NEW_CODE
		std::pair< node_index_type , node_index_type >& span = _sampleSpan[ node->nodeData.nodeIndex ];
#else // !NEW_CODE
		std::pair< int , int >& span = _sampleSpan[ node->nodeData.nodeIndex ];
#endif // NEW_CODE
		if( tree._isValidSpaceNode( node->children ) )
		{
			for( int c=0 ; c<(1<<Dim) ; c++ ) SetRange( node->children + c , start );
			span.first  = _sampleSpan[ node->children[0           ].nodeData.nodeIndex ].first;
			span.second = _sampleSpan[ node->children[ (1<<Dim)-1 ].nodeData.nodeIndex ].second;
		}
		else
		{
			span.second = start + span.second - span.first;
			span.first = start;
			start += span.second - span.first;
		}
	};

#ifdef NEW_CODE
	node_index_type start = 0;
#else // !NEW_CODE
	int start = 0;
#endif // NEW_CODE
	SetRange( tree._spaceRoot , start );
	for( FEMTreeNode* node=tree._spaceRoot->nextNode() ; node ; node=tree._spaceRoot->nextNode(node) )
		if( tree._isValidSpaceNode( node ) && !tree._isValidSpaceNode( node->children ) ) _sampleSpan[ node->nodeData.nodeIndex ].second = _sampleSpan[ node->nodeData.nodeIndex ].first;

#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) )
		{
			const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
			DualPointAndDataInfo< Dim , Real , Data , T , PointD >& _pData = _iData[ _sampleSpan[ leaf->nodeData.nodeIndex ].second++ ];
			_pData.pointInfo.position = pData.data;
			_pData.pointInfo.weight = pData.weight;
			_pData.pointInfo.dualValues = _constraintDual( pData.data/pData.weight , sampleData[i]/pData.weight ) * pData.weight;
			_pData.data = sampleData[i];
		}
	}

#ifdef NEW_THREADS
	tp.parallel_for( 0 , _iData.size() , [&]( unsigned int , size_t i  )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)_iData.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<(int)_iData.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		Real w = _iData[i].pointInfo.weight;
		_iData[i] /= w;
		if( noRescale ) _iData[i].pointInfo.weight = w;
		else            _iData[i].pointInfo.weight = w * ( 1<<tree._maxDepth );
		_iData[i].pointInfo.dualValues *= _iData[i].pointInfo.weight;
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< unsigned int Dim , class Real >
template< typename T , unsigned int PointD , typename ConstraintDual , typename SystemDual >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::ExactPointInterpolationInfo< T , PointD , ConstraintDual , SystemDual >::_init( ThreadPool &tp , const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , bool noRescale )
#else // !NEW_THREADS
void FEMTree< Dim , Real >::ExactPointInterpolationInfo< T , PointD , ConstraintDual , SystemDual >::_init( const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , bool noRescale )
#endif // NEW_THREADS
{
	_sampleSpan.resize( tree.nodesSize() );
#ifdef NEW_CODE
#ifdef NEW_THREADS
	tp.parallel_for( 0 , tree.nodesSize() , [&]( unsigned int , size_t i ){ _sampleSpan[i] = std::pair< node_index_type , node_index_type >( 0 , 0 ); } );
#else // !NEW_THREADS
#pragma omp parallel for
	for( node_index_type i=0 ; i<(node_index_type)tree.nodesSize() ; i++ ) _sampleSpan[i] = std::pair< node_index_type , node_index_type >( 0 , 0 );
#endif // NEW_THREADS
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
#else // !NEW_CODE
#pragma omp parallel for
	for( int i=0 ; i<tree.nodesSize() ; i++ ) _sampleSpan[i] = std::pair< int , int >( 0 , 0 );
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) ) _sampleSpan[ leaf->nodeData.nodeIndex ].second++;
	}
	_iData.resize( samples.size() );

#ifdef NEW_CODE
	std::function< void ( FEMTreeNode* , node_index_type & ) > SetRange = [&] ( FEMTreeNode* node , node_index_type &start )
#else // !NEW_CODE
	std::function< void ( FEMTreeNode* , int & ) > SetRange = [&] ( FEMTreeNode* node , int &start )
#endif // NEW_CODE
	{
#ifdef NEW_CODE
		std::pair< node_index_type , node_index_type >& span = _sampleSpan[ node->nodeData.nodeIndex ];
#else // !NEW_CODE
		std::pair< int , int >& span = _sampleSpan[ node->nodeData.nodeIndex ];
#endif // NEW_CODE
		if( tree._isValidSpaceNode( node->children ) )
		{
			for( int c=0 ; c<(1<<Dim) ; c++ ) SetRange( node->children + c , start );
			span.first  = _sampleSpan[ node->children[0           ].nodeData.nodeIndex ].first;
			span.second = _sampleSpan[ node->children[ (1<<Dim)-1 ].nodeData.nodeIndex ].second;
		}
		else
		{
			span.second = start + span.second - span.first;
			span.first = start;
			start += span.second - span.first;
		}
	};

#ifdef NEW_CODE
	node_index_type start=0;
#else // !NEW_CODE
	int start = 0;
#endif // NEW_CODE
	SetRange( tree._spaceRoot , start );
	for( FEMTreeNode* node=tree._spaceRoot->nextNode() ; node ; node=tree._spaceRoot->nextNode(node) )
		if( tree._isValidSpaceNode( node ) && !tree._isValidSpaceNode( node->children ) ) _sampleSpan[ node->nodeData.nodeIndex ].second = _sampleSpan[ node->nodeData.nodeIndex ].first;

#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) )
		{
			const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
			DualPointInfo< Dim , Real , T , PointD >& _pData = _iData[ _sampleSpan[ leaf->nodeData.nodeIndex ].second++ ];
			_pData.position = pData.data;
			_pData.dualValues = _constraintDual( pData.data/pData.weight ) * pData.weight;
			_pData.weight = pData.weight;
		}
	}

#ifdef NEW_THREADS
	tp.parallel_for( 0 , _iData.size() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)_iData.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<(int)_iData.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		Real w = _iData[i].weight;
		_iData[i] /= w;
		if( noRescale ) _iData[i].weight = w;
		else            _iData[i].weight = w * ( 1<<tree._maxDepth );
		_iData[i].dualValues *= _iData[i].weight;
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}
template< unsigned int Dim , class Real >
template< unsigned int PointD , typename ConstraintDual , typename SystemDual >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::ExactPointInterpolationInfo< double , PointD , ConstraintDual , SystemDual >::_init( ThreadPool &tp , const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , bool noRescale )
#else // !NEW_THREADS
void FEMTree< Dim , Real >::ExactPointInterpolationInfo< double , PointD , ConstraintDual , SystemDual >::_init( const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , bool noRescale )
#endif // NEW_THREADS
{
	_sampleSpan.resize( tree.nodesSize() );
#ifdef NEW_CODE
#ifdef NEW_THREADS
	tp.parallel_for( 0 , tree.nodesSize() , [&]( unsigned int , size_t i ){ _sampleSpan[i] = std::pair< node_index_type , node_index_type >( 0 , 0 ); } );
#else // !NEW_THREADS
#pragma omp parallel for
	for( node_index_type i=0 ; i<tree.nodesSize() ; i++ ) _sampleSpan[i] = std::pair< node_index_type , node_index_type >( 0 , 0 );
#endif // NEW_THREADS
	for( node_index_type i=0 ; i<samples.size() ; i++ )
#else // !NEW_CODE
#pragma omp parallel for
	for( int i=0 ; i<tree.nodesSize() ; i++ ) _sampleSpan[i] = std::pair< int , int >( 0 , 0 );
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) ) _sampleSpan[ leaf->nodeData.nodeIndex ].second++;
	}
	_iData.resize( samples.size() );

#ifdef NEW_CODE
	std::function< void ( FEMTreeNode* , node_index_type & ) > SetRange = [&] ( FEMTreeNode *node , node_index_type &start )
#else // !NEW_CODE
	std::function< void ( FEMTreeNode* , int & ) > SetRange = [&] ( FEMTreeNode* node , int& start )
#endif // NEW_CODE
	{
#ifdef NEW_CODE
		std::pair< node_index_type , node_index_type >& span = _sampleSpan[ node->nodeData.nodeIndex ];
#else // !NEW_CODE
		std::pair< int , int >& span = _sampleSpan[ node->nodeData.nodeIndex ];
#endif // NEW_CODE
		if( tree._isValidSpaceNode( node->children ) )
		{
			for( int c=0 ; c<(1<<Dim) ; c++ ) SetRange( node->children + c , start );
			span.first  = _sampleSpan[ node->children[0           ].nodeData.nodeIndex ].first;
			span.second = _sampleSpan[ node->children[ (1<<Dim)-1 ].nodeData.nodeIndex ].second;
		}
		else
		{
			span.second = start + span.second - span.first;
			span.first = start;
			start += span.second - span.first;
		}
	};

#ifdef NEW_CODE
	node_index_type start = 0;
#else // !NEW_CODE
	int start = 0;
#endif // NEW_CODE
	SetRange( tree._spaceRoot , start );
	for( FEMTreeNode* node=tree._spaceRoot->nextNode() ; node ; node=tree._spaceRoot->nextNode(node) )
		if( tree._isValidSpaceNode( node ) && !tree._isValidSpaceNode( node->children ) ) _sampleSpan[ node->nodeData.nodeIndex ].second = _sampleSpan[ node->nodeData.nodeIndex ].first;

#ifdef NEW_CODE
	for( node_index_type i=0 ; i<samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) )
		{
			const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
			DualPointInfo< Dim , Real , T , PointD >& _pData = _iData[ _sampleSpan[ leaf->nodeData.nodeIndex ].second++ ];
			_pData.position = pData.data;
			_pData.dualValues = _constraintDual( pData.data/pData.weight ) * pData.weight;
			_pData.weight = pData.weight;
		}
	}

#ifdef NEW_THREADS
	tp.parallel_for( 0 , _iData.size() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<_iData.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<(int)_iData.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		Real w = _iData[i].weight;
		_iData[i] /= w;
		if( noRescale ) _iData[i].weight = w;
		else            _iData[i].weight = w * ( 1<<tree._maxDepth );
		_iData[i].dualValues *= _iData[i].weight;
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}
template< unsigned int Dim , class Real >
template< typename T >
bool FEMTree< Dim , Real >::_setInterpolationInfoFromChildren( FEMTreeNode* node , SparseNodeData< T , IsotropicUIntPack< Dim , FEMTrivialSignature > >& interpolationInfo ) const
{
	if( IsActiveNode< Dim >( node->children ) )
	{
		bool hasChildData = false;
		T t = {};
		for( int c=0 ; c<(1<<Dim) ; c++ )
			if( _setInterpolationInfoFromChildren( node->children + c , interpolationInfo ) )
			{
				t += interpolationInfo[ node->children + c ];
				hasChildData = true;
			}
		if( hasChildData && IsActiveNode< Dim >( node ) ) interpolationInfo[ node ] += t;
		return hasChildData;
	}
	else return interpolationInfo( node )!=NULL;
}
template< unsigned int Dim , class Real >
template< typename T , unsigned int PointD , typename ConstraintDual >
#ifdef NEW_THREADS
SparseNodeData< DualPointInfo< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyInterpolationInfoAndSetDualConstraints( ThreadPool &tp , const std::vector< PointSample >& samples , ConstraintDual constraintDual , int adaptiveExponent ) const
#else // !NEW_THREADS
SparseNodeData< DualPointInfo< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyInterpolationInfoAndSetDualConstraints( const std::vector< PointSample >& samples , ConstraintDual constraintDual , int adaptiveExponent ) const
#endif // NEW_THREADS
{
	SparseNodeData< DualPointInfo< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > iInfo;
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* node = samples[i].node;
		const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
		while( !IsActiveNode< Dim >( node ) ) node = node->parent;
		if( pData.weight )
		{
			DualPointInfo< Dim , Real , T , PointD >& _pData = iInfo[node];
			_pData.position += pData.data;
			_pData.weight += pData.weight;
			_pData.dualValues += constraintDual( pData.data/pData.weight ) * pData.weight;
		}
	}

	// Set the interior values
	_setInterpolationInfoFromChildren( _spaceRoot , iInfo );

#ifdef NEW_THREADS
	tp.parallel_for( 0 , iInfo.size() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)iInfo.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<(int)iInfo.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THRADS
	{
		Real w = iInfo[i].weight;
		iInfo[i] /= w ; iInfo[i].weight = w;
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	LocalDepth maxDepth = _spaceRoot->maxDepth();

	// Set the average position and scale the weights
	for( const FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode(node) ) if( IsActiveNode< Dim >( node ) )
	{
		DualPointInfo< Dim , Real , T , PointD >* pData = iInfo( node );
		if( pData )
		{
			int e = _localDepth( node ) * adaptiveExponent - ( maxDepth ) * (adaptiveExponent-1);
			if( e<0 ) pData->weight /= Real( 1<<(-e) );
			else      pData->weight *= Real( 1<<  e  );
			pData->dualValues *= pData->weight;
		}
	}
	return iInfo;
}
template< unsigned int Dim , class Real >
template< typename T , typename Data , unsigned int PointD , typename ConstraintDual >
#ifdef NEW_THREADS
SparseNodeData< DualPointAndDataInfo< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyInterpolationInfoAndSetDualConstraints( ThreadPool &tp , const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , ConstraintDual constraintDual , int adaptiveExponent ) const
#else // !NEW_THREADS
SparseNodeData< DualPointAndDataInfo< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyInterpolationInfoAndSetDualConstraints( const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , ConstraintDual constraintDual , int adaptiveExponent ) const
#endif // NEW_THREADS
{
	SparseNodeData< DualPointAndDataInfo< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > iInfo;
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* node = samples[i].node;
		const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
		while( !IsActiveNode< Dim >( node ) ) node = node->parent;
		if( pData.weight )
		{
			DualPointAndDataInfo< Dim , Real , Data , T , PointD >& _pData = iInfo[node];
			_pData.pointInfo.position += pData.data;
			_pData.pointInfo.dualValues += constraintDual( pData.data/pData.weight , sampleData[i]/pData.weight ) * pData.weight;
			_pData.pointInfo.weight += pData.weight;
			_pData.data += sampleData[i];
		}
	}

	// Set the interior values
	_setInterpolationInfoFromChildren( _spaceRoot , iInfo );

#ifdef NEW_THREADS
	tp.parallel_for( 0 , iInfo.size() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)iInfo.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<(int)iInfo.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		Real w = iInfo[i].pointInfo.weight;
		iInfo[i] /= w ; iInfo[i].pointInfo.weight = w;
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	LocalDepth maxDepth = _spaceRoot->maxDepth();

	// Set the average position and scale the weights
	for( const FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode(node) ) if( IsActiveNode< Dim >( node ) )
	{
		DualPointAndDataInfo< Dim , Real , Data , T , PointD >* pData = iInfo( node );
		if( pData )
		{
			int e = _localDepth( node ) * adaptiveExponent - ( maxDepth ) * (adaptiveExponent-1);
			if( e<0 ) pData->pointInfo.weight /= Real( 1<<(-e) );
			else      pData->pointInfo.weight *= Real( 1<<  e  );
			pData->pointInfo.dualValues *= pData->pointInfo.weight;
		}
	}
	return iInfo;
}
template< unsigned int Dim , class Real >
template< typename T , unsigned int PointD , typename ConstraintDual >
#ifdef NEW_THREADS
SparseNodeData< DualPointInfoBrood< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyChildInterpolationInfoAndSetDualConstraints( ThreadPool &tp , const std::vector< PointSample >& samples , ConstraintDual constraintDual , bool noRescale ) const
#else // !NEW_THREADS
SparseNodeData< DualPointInfoBrood< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyChildInterpolationInfoAndSetDualConstraints( const std::vector< PointSample >& samples , ConstraintDual constraintDual , bool noRescale ) const
#endif // NEW_THREADS
{
	SparseNodeData< DualPointInfoBrood< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > iInfo;
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* node = samples[i].node;
		const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
		while( !IsActiveNode< Dim >( node ) ) node = node->parent;
		if( pData.weight )
		{
			DualPointInfoBrood< Dim , Real , T , PointD >& _pData = iInfo[node];
			Point< Real , Dim > p = pData.data/pData.weight;
			int cIdx = _childIndex( node , p );
			_pData[cIdx].position += pData.data;
			_pData[cIdx].weight += pData.weight;
			_pData[cIdx].dualValues += constraintDual( p ) * pData.weight;
		}
	}

	// Set the interior values
	_setInterpolationInfoFromChildren( _spaceRoot , iInfo );

#ifdef NEW_THREADS
	tp.parallel_for( 0 , iInfo.size() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<iInfo.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<(int)iInfo.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		iInfo[i].finalize();
#ifdef NEW_CODE
		for( size_t c=0 ; c<iInfo[i].size() ; c++ )
#else // !NEW_CODE
		for( int c=0 ; c<(int)iInfo[i].size() ; c++ )
#endif // NEW_CODE
		{
			iInfo[i][c].position /= iInfo[i][c].weight;
			if( !noRescale )
			{
				iInfo[i][c].weight     *= ( 1<<_maxDepth );
				iInfo[i][c].dualValues *= ( 1<<_maxDepth );
			}
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	return iInfo;
}
template< unsigned int Dim , class Real >
template< typename T , typename Data , unsigned int PointD , typename ConstraintDual >
#ifdef NEW_THREADS
SparseNodeData< DualPointAndDataInfoBrood< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyChildInterpolationInfoAndSetDualConstraints( ThreadPool &tp , const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , ConstraintDual constraintDual , bool noRescale ) const
#else // !NEW_THREADS
SparseNodeData< DualPointAndDataInfoBrood< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyChildInterpolationInfoAndSetDualConstraints( const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , ConstraintDual constraintDual , bool noRescale ) const
#endif // NEW_THREADS
{
	SparseNodeData< DualPointAndDataInfoBrood< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > iInfo;
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<samples.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<samples.size() ; i++ )
#endif // NEW_CODE
	{
		const FEMTreeNode* node = samples[i].node;
		const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
		while( !IsActiveNode< Dim >( node ) ) node = node->parent;
		if( pData.weight )
		{
			DualPointAndDataInfoBrood< Dim , Real , Data , T , PointD >& _pData = iInfo[node];
			Point< Real , Dim > p = pData.data/pData.weight;
			int cIdx = _childIndex( node , p );
			_pData[cIdx].pointInfo.position += pData.data;
			_pData[cIdx].pointInfo.dualValues += constraintDual( p , sampleData[i]/pData.weight ) * pData.weight;
			_pData[cIdx].pointInfo.weight += pData.weight;
			_pData[cIdx].data += sampleData[i];
		}
	}

	// Set the interior values
	_setInterpolationInfoFromChildren( _spaceRoot , iInfo );

#ifdef NEW_THREADS
	tp.parallel_for( 0 , iInfo.size() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<iInfo.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<(int)iInfo.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		iInfo[i].finalize();
#ifdef NEW_CODE
		for( size_t c=0 ; c<iInfo[i].size() ; c++ )
#else // !NEW_CODE
		for( int c=0 ; c<(int)iInfo[i].size() ; c++ )
#endif // NEW_CODE
		{
			iInfo[i][c].pointInfo.position /= iInfo[i][c].pointInfo.weight;
			iInfo[i][c].data /= iInfo[i][c].pointInfo.weight;
			if( !noRescale )
			{
				iInfo[i][c].pointInfo.weight     *= ( 1<<_maxDepth );
				iInfo[i][c].pointInfo.dualValues *= ( 1<<_maxDepth );
				iInfo[i][c].data                 *= ( 1<<_maxDepth );
			}
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	return iInfo;
}



template< unsigned int Dim , class Real >
#ifdef NEW_CODE
std::vector< node_index_type > FEMTree< Dim , Real >::merge( FEMTree* tree )
#else // !NEW_CODE
std::vector< int > FEMTree< Dim , Real >::merge( FEMTree* tree )
#endif // NEW_CODE
{
#ifdef NEW_CODE
	std::vector< node_index_type > map;
#else // !NEW_CODE
	std::vector< int > map;
#endif // NEW_CODE
	if( _depthOffset!=tree->_depthOffset ) ERROR_OUT( "depthOffsets don't match: %d != %d" , _depthOffset , tree->_depthOffset );

	// Compute the next available index
#ifdef NEW_CODE
	node_index_type nextIndex = 0;
	for( const FEMTreeNode* node=_tree->nextNode() ; node!=NULL ; node=_tree->nextNode( node ) ) nextIndex = std::max< node_index_type >( nextIndex , node->nodeData.nodeIndex+1 );
#else // !NEW_CODE
	int nextIndex = 0;
	for( const FEMTreeNode* node=_tree->nextNode() ; node!=NULL ; node=_tree->nextNode( node ) ) nextIndex = std::max< int >( nextIndex , node->nodeData.nodeIndex+1 );
#endif // NEW_CODE

	// Set the size of the map
	{
#ifdef NEW_CODE
		node_index_type mapSize = 0;
		for( const FEMTreeNode* node=tree->_tree->nextNode() ; node!=NULL ; node=tree->_tree->nextNode( node ) ) mapSize = std::max< node_index_type >( mapSize , node->nodeData.nodeIndex+1 );
#else // !NEW_CODE
		int mapSize = 0;
		for( const FEMTreeNode* node=tree->_tree->nextNode() ; node!=NULL ; node=tree->_tree->nextNode( node ) ) mapSize = std::max< int >( mapSize , node->nodeData.nodeIndex+1 );
#endif // NEW_CODE
		map.resize( mapSize );
	}

#ifdef NEW_CODE
	std::function< void ( FEMTreeNode* , FEMTreeNode* , std::vector< node_index_type > & , node_index_type & ) > MergeNodes = [&]( FEMTreeNode* node1 , FEMTreeNode* node2 , std::vector< node_index_type > &map , node_index_type &nextIndex )
#else // !NEW_CODE
	std::function< void ( FEMTreeNode* , FEMTreeNode* , std::vector< int > & , int & ) > MergeNodes = [&]( FEMTreeNode* node1 , FEMTreeNode* node2 , std::vector< int > &map , int &nextIndex )
#endif // NEW_CODE
	{
		if( node1 && node2 )
		{
			if( node2->nodeData.nodeIndex>=0 )
			{
				if( node1->nodeData.nodeIndex<0 ) node1->nodeData.nodeIndex = nextIndex++;
				map[ node2->nodeData.nodeIndex ] = node1->nodeData.nodeIndex;
			}
			if( node1->children && node2->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) MergeNodes( node1->children+c , node2->children+c , map , nextIndex );
			else if( node2->children )
			{
				for( int c=0 ; c<(1<<Dim) ; c++ ) MergeNodes( NULL , node2->children+c , map , nextIndex );
				node1->children = node2->children;
				node2->children = NULL;
				for( int c=0 ; c<(1<<Dim) ; c++ ) node1->children[c].parent = node1;
			}
		}
		else if( node2 )
		{
			if( node2->nodeData.nodeIndex>=0 ){ map[ node2->nodeData.nodeIndex ] = nextIndex ; node2->nodeData.nodeIndex = nextIndex++; }
			if( node2->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) MergeNodes( NULL , node2->children+c , map , nextIndex );
		}
	};

	MergeNodes( _tree , tree->_tree , map , nextIndex );
	return map;
}

