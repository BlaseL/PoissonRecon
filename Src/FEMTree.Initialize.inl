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

////////////////////////
// FEMTreeInitializer //
////////////////////////
template< unsigned int Dim , class Real >
#ifdef NEW_CODE
size_t FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& node , int maxDepth , std::function< bool ( int , int[] ) > Refine , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#else // !NEW_CODE
int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& node , int maxDepth , std::function< bool ( int , int[] ) > Refine , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#endif // NEW_CODE
{
#ifdef NEW_CODE
	size_t count = 0;
#else // !NEW_CODE
	int count = 0;
#endif // NEW_CODE
	int d , off[3];
	node.depthAndOffset( d , off );
	if( node.depth()<maxDepth && Refine( d , off ) )
	{
#ifdef THREAD_SAFE_CHILD_INIT
		node.initChildren< false >( nodeAllocator , NodeInitializer ) , count += 1<<Dim;
#else // !THREAD_SAFE_CHILD_INIT
		node.initChildren( nodeAllocator , NodeInitializer ) , count += 1<<Dim;
#endif // THREAD_SAFE_CHILD_INIT
		for( int c=0 ; c<(1<<Dim) ; c++ ) count += Initialize( node.children[c] , maxDepth , Refine , nodeAllocator , NodeInitializer );
	}
	return count;
}

template< unsigned int Dim , class Real >
#ifdef NEW_CODE
size_t FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , InputPointStream< Real , Dim >& pointStream , int maxDepth , std::vector< PointSample >& samplePoints , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#else // !NEW_CODE
int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , InputPointStream< Real , Dim >& pointStream , int maxDepth , std::vector< PointSample >& samplePoints , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#endif // NEW_CODE
{
	auto Leaf = [&]( FEMTreeNode& root , Point< Real , Dim > p , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return (FEMTreeNode*)NULL;
		Point< Real , Dim > center;
		for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
		Real width = Real(1.0);
		FEMTreeNode* node = &root;
		int d = 0;
		while( d<maxDepth )
		{
#ifdef THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->template initChildren< false >( nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
			int cIndex = FEMTreeNode::ChildIndex( center , p );
			node = node->children + cIndex;
			d++;
			width /= 2;
			for( int dd=0 ; dd<Dim ; dd++ )
				if( (cIndex>>dd) & 1 ) center[dd] += width/2;
				else                   center[dd] -= width/2;
		}
		return node;
	};

	// Add the point data
#ifdef NEW_CODE
	size_t outOfBoundPoints = 0 , pointCount = 0;
#else // !NEW_CODE
	int outOfBoundPoints = 0 , pointCount = 0;
#endif // NEW_CODE
	{
#ifdef NEW_CODE
		std::vector< node_index_type > nodeToIndexMap;
#else // !NEW_CODE
		std::vector< int > nodeToIndexMap;
#endif // NEW_CODE
		Point< Real , Dim > p;
		while( pointStream.nextPoint( p ) )
		{
			Real weight = (Real)1.;
			FEMTreeNode* temp = Leaf( root , p , maxDepth );
			if( !temp ){ outOfBoundPoints++ ; continue; }
#ifdef NEW_CODE
			node_index_type nodeIndex = temp->nodeData.nodeIndex;
			if( nodeIndex>=(node_index_type)nodeToIndexMap.size() ) nodeToIndexMap.resize( nodeIndex+1 , -1 );
			node_index_type idx = nodeToIndexMap[ nodeIndex ];
#else // !NEW_CODE
			int nodeIndex = temp->nodeData.nodeIndex;
			if( nodeIndex>=nodeToIndexMap.size() ) nodeToIndexMap.resize( nodeIndex+1 , -1 );
			int idx = nodeToIndexMap[ nodeIndex ];
#endif // NEW_CODE
			if( idx==-1 )
			{
#ifdef NEW_CODE
				idx = (node_index_type)samplePoints.size();
#else // !NEW_CODE
				idx = (int)samplePoints.size();
#endif // NEW_CODE
				nodeToIndexMap[ nodeIndex ] = idx;
				samplePoints.resize( idx+1 ) , samplePoints[idx].node = temp;
			}
			samplePoints[idx].sample += ProjectiveData< Point< Real , Dim > , Real >( p*weight , weight );
			pointCount++;
		}
		pointStream.reset();
	}
	if( outOfBoundPoints  ) WARN( "Found out-of-bound points: " , outOfBoundPoints );
	FEMTree< Dim , Real >::MemoryUsage();
	return pointCount;
}

template< unsigned int Dim , class Real >
template< class Data >
#ifdef NEW_CODE
size_t FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , InputPointStreamWithData< Real , Dim , Data >& pointStream , int maxDepth , std::vector< PointSample >& samplePoints , std::vector< Data >& sampleData , bool mergeNodeSamples , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer , std::function< Real ( const Point< Real , Dim >& , Data& ) > ProcessData )
#else // !NEW_CODE
int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , InputPointStreamWithData< Real , Dim , Data >& pointStream , int maxDepth , std::vector< PointSample >& samplePoints , std::vector< Data >& sampleData , bool mergeNodeSamples , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer , std::function< Real ( const Point< Real , Dim >& , Data& ) > ProcessData )
#endif // NEW_CODE
{
	auto Leaf = [&]( FEMTreeNode& root , Point< Real , Dim > p , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return (FEMTreeNode*)NULL;
		Point< Real , Dim > center;
		for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
		Real width = Real(1.0);
		FEMTreeNode* node = &root;
		int d = 0;
		while( d<maxDepth )
		{
#ifdef THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->template initChildren< false >( nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
			int cIndex = FEMTreeNode::ChildIndex( center , p );
			node = node->children + cIndex;
			d++;
			width /= 2;
			for( int dd=0 ; dd<Dim ; dd++ )
				if( (cIndex>>dd) & 1 ) center[dd] += width/2;
				else                   center[dd] -= width/2;
		}
		return node;
	};

	// Add the point data
#ifdef NEW_CODE
	size_t outOfBoundPoints = 0 , badData = 0 , pointCount = 0;
#else //!NEW_CODE
	int outOfBoundPoints = 0 , badData = 0 , pointCount = 0;
#endif // NEW_CODE
	{
#ifdef NEW_CODE
		std::vector< node_index_type > nodeToIndexMap;
#else // !NEW_CODE
		std::vector< int > nodeToIndexMap;
#endif // NEW_CODE
		Point< Real , Dim > p;
		Data d;

		while( pointStream.nextPoint( p , d ) )
		{
			Real weight = ProcessData( p , d );
			if( weight<=0 ){ badData++ ; continue; }
			FEMTreeNode* temp = Leaf( root , p , maxDepth );
			if( !temp ){ outOfBoundPoints++ ; continue; }
#ifdef NEW_CODE
			node_index_type nodeIndex = temp->nodeData.nodeIndex;
#else // !NEW_CODE
			int nodeIndex = temp->nodeData.nodeIndex;
#endif // NEW_CODE
			if( mergeNodeSamples )
			{
#ifdef NEW_CODE
				if( nodeIndex>=(node_index_type)nodeToIndexMap.size() ) nodeToIndexMap.resize( nodeIndex+1 , -1 );
				node_index_type idx = nodeToIndexMap[ nodeIndex ];
#else // !NEW_CODE
				if( nodeIndex>=nodeToIndexMap.size() ) nodeToIndexMap.resize( nodeIndex+1 , -1 );
				int idx = nodeToIndexMap[ nodeIndex ];
#endif // NEW_CODE
				if( idx==-1 )
				{
#ifdef NEW_CODE
					idx = (node_index_type)samplePoints.size();
#else // !NEW_CODE
					idx = (int)samplePoints.size();
#endif // NEW_CODE
					nodeToIndexMap[ nodeIndex ] = idx;
					samplePoints.resize( idx+1 ) , samplePoints[idx].node = temp;
					sampleData.resize( idx+1 );
				}
				samplePoints[idx].sample += ProjectiveData< Point< Real , Dim > , Real >( p*weight , weight );
				sampleData[ idx ] += d*weight;
			}
			else
			{
#ifdef NEW_CODE
				node_index_type idx = (node_index_type)samplePoints.size();
#else // !NEW_CODE
				int idx = (int)samplePoints.size();
#endif // NEW_CODE
				samplePoints.resize( idx+1 ) , sampleData.resize( idx+1 );
				samplePoints[idx].node = temp;
				samplePoints[idx].sample = ProjectiveData< Point< Real , Dim > , Real >( p*weight , weight );
				sampleData[ idx ] = d*weight;
			}
			pointCount++;
		}
		pointStream.reset();
	}
	if( outOfBoundPoints  ) WARN( "Found out-of-bound points: " , outOfBoundPoints );
	if( badData           ) WARN( "Found bad data: " , badData );
	FEMTree< Dim , Real >::MemoryUsage();
	return pointCount;
}
template< unsigned int Dim , class Real >
#ifdef NEW_CODE
#ifdef NEW_THREADS
void FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , const std::vector< Point< Real , Dim > >& vertices , const std::vector< SimplexIndex< Dim-1 , node_index_type > >& simplices , int maxDepth , std::vector< PointSample >& samples , bool mergeNodeSamples , std::vector< Allocator< FEMTreeNode > * > &nodeAllocators , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#else // !NEW_THREADS
void FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , const std::vector< Point< Real , Dim > >& vertices , const std::vector< SimplexIndex< Dim-1 , node_index_type > >& simplices , int maxDepth , std::vector< PointSample >& samples , bool mergeNodeSamples , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#endif // NEW_THREADS
#else // !NEW_CODE
void FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , const std::vector< Point< Real , Dim > >& vertices , const std::vector< SimplexIndex< Dim-1 > >& simplices , int maxDepth , std::vector< PointSample >& samples , bool mergeNodeSamples , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#endif // NEW_CODE
{
#ifdef NEW_CODE
	std::vector< node_index_type > nodeToIndexMap;
#else // !NEW_CODE
	std::vector< int > nodeToIndexMap;
#endif // NEW_CODE
#ifdef NEW_THREADS
	ThreadPool::Parallel_for( 0 , simplices.size() , [&]( unsigned int t , size_t  i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( long long i=0 ; i<(long long)simplices.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<simplices.size() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		Simplex< Real , Dim , Dim-1 > s;
		for( int k=0 ; k<Dim ; k++ ) s[k] = vertices[ simplices[i][k] ];
#ifdef NEW_CODE
#ifdef NEW_THREADS
#ifdef THREAD_SAFE_CHILD_INIT
		if( mergeNodeSamples ) _AddSimplex< true >( root , s , maxDepth , samples , &nodeToIndexMap , nodeAllocators.size() ? nodeAllocators[t] : NULL , NodeInitializer );
		else                   _AddSimplex< true >( root , s , maxDepth , samples , NULL ,            nodeAllocators.size() ? nodeAllocators[t] : NULL , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		if( mergeNodeSamples ) _AddSimplex( root , s , maxDepth , samples , &nodeToIndexMap , nodeAllocators.size() ? nodeAllocators[t] : NULL , NodeInitializer );
		else                   _AddSimplex( root , s , maxDepth , samples , NULL ,            nodeAllocators.size() ? nodeAllocators[t] : NULL , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
#else // !NEW_THREADS
#ifdef THREAD_SAFE_CHILD_INIT
		if( mergeNodeSamples ) _AddSimplex< true >( root , s , maxDepth , samples , &nodeToIndexMap , nodeAllocator , NodeInitializer );
		else                   _AddSimplex< true >( root , s , maxDepth , samples , NULL ,            nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		if( mergeNodeSamples ) _AddSimplex( root , s , maxDepth , samples , &nodeToIndexMap , nodeAllocator , NodeInitializer );
		else                   _AddSimplex( root , s , maxDepth , samples , NULL ,            nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
#endif // NEW_THREADS
#else // !NEW_CODE
		int sCount;
#ifdef THREAD_SAFE_CHILD_INIT
		if( mergeNodeSamples ) sCount = _AddSimplex< true >( root , s , maxDepth , samples , &nodeToIndexMap , nodeAllocator , NodeInitializer );
		else                   sCount = _AddSimplex< true >( root , s , maxDepth , samples , NULL ,            nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		if( mergeNodeSamples ) sCount = _AddSimplex( root , s , maxDepth , samples , &nodeToIndexMap , nodeAllocator , NodeInitializer );
		else                   sCount = _AddSimplex( root , s , maxDepth , samples , NULL ,            nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
#endif // NEW_CODE
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	FEMTree< Dim , Real >::MemoryUsage();
}

template< unsigned int Dim , class Real >
#ifdef NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
template< bool ThreadSafe >
#endif // THREAD_SAFE_CHILD_INIT
size_t FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode& root , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< PointSample >& samples , std::vector< node_index_type >* nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#else // !NEW_CODE
int FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode& root , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< PointSample >& samples , std::vector< int >* nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#endif // NEW_CODE
{
	std::vector< Simplex< Real , Dim , Dim-1 > > subSimplices;
	subSimplices.push_back( s );

	// Clip the simplex to the unit cube
	{
		for( int d=0 ; d<Dim ; d++ )
		{
			Point< Real , Dim > n;
			n[d] = 1;
			{
				std::vector< Simplex< Real , Dim , Dim-1 > > back , front;
				for( int i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 0 , back , front );
				subSimplices = front;
			}
			{
				std::vector< Simplex< Real , Dim , Dim-1 > > back , front;
				for( int i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 1 , back , front );
				subSimplices = back;
			}
		}
	}

	struct RegularGridIndex
	{
		int idx[Dim];
		bool operator != ( const RegularGridIndex& i ) const
		{
			for( int d=0 ; d<Dim ; d++ ) if( idx[d]!=i.idx[d] ) return true;
			return false;
		}
	};

	auto Leaf = [&]( Point< Real , Dim > p , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return (FEMTreeNode*)NULL;
		Point< Real , Dim > center;
		for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
		Real width = Real(1.0);
		FEMTreeNode* node = &root;
		int d=0;
		while( d<maxDepth )
		{
#ifdef THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->template initChildren< ThreadSafe >( nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
#ifdef NEW_THREADS
			{
				static std::mutex m;
				std::lock_guard< std::mutex > lock( m );
				if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
			}
#else // !NEW_THREADS
#pragma omp critical
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
#endif // NEW_THREADS
#endif // THREAD_SAFE_CHILD_INIT
			int cIndex = FEMTreeNode::ChildIndex( center , p );
			node = node->children + cIndex;
			d++;
			width /= 2;
			for( int d=0 ; d<Dim ; d++ )
				if( (cIndex>>d) & 1 ) center[d] += width/2;
				else                  center[d] -= width/2;
		}
		return node;
	};


#ifdef NEW_CODE
	size_t sCount = 0;
#else // !NEW_CODE
	int sCount = 0;
#endif // NEW_CODE
	for( int i=0 ; i<subSimplices.size() ; i++ )
	{
		// Find the finest depth at which the simplex is entirely within a node
		int tDepth;
		RegularGridIndex idx0 , idx;
		for( tDepth=0 ; tDepth<maxDepth ; tDepth++ )
		{
			// Get the grid index of the first vertex of the simplex
			for( int d=0 ; d<Dim ; d++ ) idx0.idx[d] = idx.idx[d] = (int)( subSimplices[i][0][d] * (1<<(tDepth+1)) );
			bool done = false;
			for( int k=0 ; k<Dim && !done ; k++ )
			{
				for( int d=0 ; d<Dim ; d++ ) idx.idx[d] = (int)( subSimplices[i][k][d] * (1<<(tDepth+1)) );
				if( idx!=idx0 ) done = true;
			}
			if( done ) break;
		}

		// Generate a point in the middle of the simplex
#ifdef THREAD_SAFE_CHILD_INIT
		for( int i=0 ; i<subSimplices.size() ; i++ ) sCount += _AddSimplex< ThreadSafe >( Leaf( subSimplices[i].center() , tDepth ) , subSimplices[i] , maxDepth , samples , nodeToIndexMap , nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		for( int i=0 ; i<subSimplices.size() ; i++ ) sCount += _AddSimplex( Leaf( subSimplices[i].center() , tDepth ) , subSimplices[i] , maxDepth , samples , nodeToIndexMap , nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
	}
	return sCount;
}
template< unsigned int Dim , class Real >
#ifdef NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
template< bool ThreadSafe >
#endif // THREAD_SAFE_CHILD_INIT
size_t FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode* node , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< PointSample >& samples , std::vector< node_index_type >* nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#else // !NEW_CODE
int FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode* node , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< PointSample >& samples , std::vector< int >* nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#endif // NEW_CODE
{
	int d = node->depth();
	if( d==maxDepth )
	{
		Real weight = s.measure();
		Point< Real , Dim > position = s.center() , normal;
		{
			Point< Real , Dim > v[Dim-1];
			for( int k=0 ; k<Dim-1 ; k++ ) v[k] = s[k+1]-s[0];
			normal = Point< Real , Dim >::CrossProduct( v );
		}
		if( weight && weight==weight )
		{
			if( nodeToIndexMap )
			{
#ifdef NEW_CODE
				node_index_type nodeIndex = node->nodeData.nodeIndex;
#else // !NEW_CODE
				int nodeIndex = node->nodeData.nodeIndex;
#endif // NEW_CODE
#ifdef NEW_THREADS
				{
					static std::mutex m;
					std::lock_guard< std::mutex > lock(m);
#else // !NEW_THREADS
#pragma omp critical
				{
#endif // NEW_THREADS
					if( nodeIndex>=(node_index_type)nodeToIndexMap->size() ) nodeToIndexMap->resize( nodeIndex+1 , -1 );
#ifdef NEW_CODE
					node_index_type idx = (*nodeToIndexMap)[ nodeIndex ];
#else // !NEW_CODE
					int idx = (*nodeToIndexMap)[ nodeIndex ];
#endif // NEW_CODE
					if( idx==-1 )
					{
#ifdef NEW_CODE
						idx = (node_index_type)samples.size();
#else // !NEW_CODE
						idx = (int)samples.size();
#endif // NEW_CODE
						(*nodeToIndexMap)[ nodeIndex ] = idx;
						samples.resize( idx+1 );
						samples[idx].node = node;
					}
					samples[idx].sample += ProjectiveData< Point< Real , Dim > , Real >( position*weight , weight );
				}
			}
			else
			{
#ifdef NEW_THREADS
				{
					static std::mutex m;
					std::lock_guard< std::mutex > lock(m);
#else // !NEW_THREADS
#pragma omp critical
				{
#endif // NEW_THREADS
#ifdef NEW_CODE
					node_index_type idx = (node_index_type)samples.size();
#else // !NEW_CODE
					int idx = (int)samples.size();
#endif // NEW_CODE
					samples.resize( idx+1 );
					samples[idx].node = node;
					samples[idx].sample = ProjectiveData< Point< Real , Dim > , Real >( position*weight , weight );
				}
			}
		}
		return 1;
	}
	else
	{
#ifdef NEW_CODE
		size_t sCount = 0;
#else // !NEW_CODE
		int sCount = 0;
#endif // NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
		if( !node->children ) node->template initChildren< ThreadSafe >( nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
#ifdef NEW_THREADS
		{
			static std::mutex m;
			std::lock_guard< std::mutex > lock(m);
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
		}
#else // !NEW_THREADS
#pragma omp critical
		if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
#endif // NEW_THREADS
#endif // THREAD_SAFE_CHILD_INIT

		// Split up the simplex and pass the parts on to the children
		Point< Real , Dim > center;
		Real width;
		node->centerAndWidth( center , width );

		std::vector< std::vector< Simplex< Real , Dim , Dim-1 > > > childSimplices( 1 );
		childSimplices[0].push_back( s );
		for( int d=0 ; d<Dim ; d++ )
		{
			Point< Real , Dim > n ; n[Dim-d-1] = 1;
			std::vector< std::vector< Simplex< Real , Dim , Dim-1 > > > temp( (int)( 1<<(d+1) ) );
#ifdef NEW_CODE
			for( int c=0 ; c<(1<<d) ; c++ ) for( size_t i=0 ; i<childSimplices[c].size() ; i++ ) childSimplices[c][i].split( n , center[Dim-d-1] , temp[2*c] , temp[2*c+1] );
#else // !NEW_CODE
			for( int c=0 ; c<(1<<d) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) childSimplices[c][i].split( n , center[Dim-d-1] , temp[2*c] , temp[2*c+1] );
#endif // NEW_CODE
			childSimplices = temp;
		}
#ifdef NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
		for( int c=0 ; c<(1<<Dim) ; c++ ) for( size_t i=0 ; i<childSimplices[c].size() ; i++ ) if( childSimplices[c][i].measure() ) sCount += _AddSimplex< ThreadSafe >( node->children+c , childSimplices[c][i] , maxDepth , samples , nodeToIndexMap , nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		for( int c=0 ; c<(1<<Dim) ; c++ ) for( size_t i=0 ; i<childSimplices[c].size() ; i++ ) if( childSimplices[c][i].measure() ) sCount += _AddSimplex( node->children+c , childSimplices[c][i] , maxDepth , samples , nodeToIndexMap , nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
#else // !NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
		for( int c=0 ; c<(1<<Dim) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) if( childSimplices[c][i].measure() ) sCount += _AddSimplex< ThreadSafe >( node->children+c , childSimplices[c][i] , maxDepth , samples , nodeToIndexMap , nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		for( int c=0 ; c<(1<<Dim) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) if( childSimplices[c][i].measure() ) sCount += _AddSimplex( node->children+c , childSimplices[c][i] , maxDepth , samples , nodeToIndexMap , nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
#endif // NEW_CODE
		return sCount;
	}
}

template< unsigned int Dim , class Real >
#ifdef NEW_CODE
void FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , const std::vector< Point< Real , Dim > >& vertices , const std::vector< SimplexIndex< Dim-1 , node_index_type > >& simplices , int maxDepth , std::vector< NodeSimplices< Dim , Real > >& nodeSimplices , Allocator< FEMTreeNode > *nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#else // !NEW_CODE
void FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , const std::vector< Point< Real , Dim > >& vertices , const std::vector< SimplexIndex< Dim-1 > >& simplices , int maxDepth , std::vector< NodeSimplices< Dim , Real > >& nodeSimplices , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#endif // NEW_CODE
{
#ifdef NEW_CODE
	std::vector< size_t > nodeToIndexMap;
	for( size_t i=0 ; i<simplices.size() ; i++ )
#else // !NEW_CODE
	std::vector< int > nodeToIndexMap;
	for( int i=0 ; i<simplices.size() ; i++ )
#endif // NEW_CODE
	{
		Simplex< Real , Dim , Dim-1 > s;
		for( int k=0 ; k<Dim ; k++ ) s[k] = vertices[ simplices[i][k] ];
#ifdef THREAD_SAFE_CHILD_INIT
		_AddSimplex< false >( root , s , maxDepth , nodeSimplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		_AddSimplex( root , s , maxDepth , nodeSimplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
	}
	FEMTree< Dim , Real >::MemoryUsage();
}

template< unsigned int Dim , class Real >
#ifdef NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
template< bool ThreadSafe >
#endif // THREAD_SAFE_CHILD_INIT
size_t FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode& root , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< NodeSimplices< Dim , Real > >& simplices , std::vector< node_index_type >& nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#else // !NEW_CODE
int FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode& root , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< NodeSimplices< Dim , Real > >& simplices , std::vector< int >& nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#endif // NEW_CODE
{
	std::vector< Simplex< Real , Dim , Dim-1 > > subSimplices;
	subSimplices.push_back( s );

	// Clip the simplex to the unit cube
	{
		for( int d=0 ; d<Dim ; d++ )
		{
			Point< Real , Dim > n;
			n[d] = 1;
			{
				std::vector< Simplex< Real , Dim , Dim-1 > > back , front;
#ifdef NEW_CODE
				for( size_t i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 0 , back , front );
#else // !NEW_CODE
				for( int i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 0 , back , front );
#endif // NEW_CODE
				subSimplices = front;
			}
			{
				std::vector< Simplex< Real , Dim , Dim-1 > > back , front;
#ifdef NEW_CODE
				for( size_t i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 1 , back , front );
#else // !NEW_CODE
				for( int i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 1 , back , front );
#endif // NEW_CODE
				subSimplices = back;
			}
		}
	}

	struct RegularGridIndex
	{
		int idx[Dim];
		bool operator != ( const RegularGridIndex& i ) const
		{
			for( int d=0 ; d<Dim ; d++ ) if( idx[d]!=i.idx[d] ) return true;
			return false;
		}
	};

	auto Leaf = [&]( Point< Real , Dim > p , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return (FEMTreeNode*)NULL;
		Point< Real , Dim > center;
		for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
		Real width = Real(1.0);
		FEMTreeNode* node = &root;
		int d=0;
		while( d<maxDepth )
		{
#ifdef THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->template initChildren< ThreadSafe >( nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
			int cIndex = FEMTreeNode::ChildIndex( center , p );
			node = node->children + cIndex;
			d++;
			width /= 2;
			for( int d=0 ; d<Dim ; d++ )
				if( (cIndex>>d) & 1 ) center[d] += width/2;
				else                  center[d] -= width/2;
		}
		return node;
	};

#ifdef NEW_CODE
	size_t sCount = 0;
#else // !NEW_CODE
	int sCount = 0;
#endif // NEW_CODE

#ifdef NEW_CODE
	for( size_t i=0 ; i<subSimplices.size() ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<subSimplices.size() ; i++ )
#endif // NEW_CODE
	{
		// Find the finest depth at which the simplex is entirely within a node
		int tDepth;
		RegularGridIndex idx0 , idx;
		for( tDepth=0 ; tDepth<maxDepth ; tDepth++ )
		{
			// Get the grid index of the first vertex of the simplex
			for( int d=0 ; d<Dim ; d++ ) idx0.idx[d] = (int)( subSimplices[i][0][d] * (1<<(tDepth+1)) );
			bool done = false;
			for( int k=0 ; k<Dim && !done ; k++ )
			{
				for( int d=0 ; d<Dim ; d++ ) idx.idx[d] = (int)( subSimplices[i][k][d] * (1<<(tDepth+1)) );
				if( idx!=idx0 ) done = true;
			}
			if( done ) break;
		}

		// Add the simplex to the node
		FEMTreeNode* subSimplexNode = Leaf( subSimplices[i].center() , tDepth );
#ifdef NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
		for( size_t i=0 ; i<subSimplices.size() ; i++ ) sCount += _AddSimplex< ThreadSafe >( subSimplexNode , subSimplices[i] , maxDepth , simplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		for( size_t i=0 ; i<subSimplices.size() ; i++ ) sCount += _AddSimplex( subSimplexNode , subSimplices[i] , maxDepth , simplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
#else // !NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
		for( int i=0 ; i<subSimplices.size() ; i++ ) sCount += _AddSimplex< ThreadSafe >( subSimplexNode , subSimplices[i] , maxDepth , simplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		for( int i=0 ; i<subSimplices.size() ; i++ ) sCount += _AddSimplex( subSimplexNode , subSimplices[i] , maxDepth , simplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
#endif // NEW_CODE
	}
	return sCount;
}
template< unsigned int Dim , class Real >
#ifdef NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
template< bool ThreadSafe >
#endif // THREAD_SAFE_CHILD_INIT
size_t FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode* node , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< NodeSimplices< Dim , Real > >& simplices , std::vector< node_index_type >& nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#else // !NEW_CODE
int FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode* node , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< NodeSimplices< Dim , Real > >& simplices , std::vector< int >& nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
#endif // NEW_CODE
{
	int d = node->depth();
	if( d==maxDepth )
	{
		// If the simplex has non-zero size, add it to the list
		Real weight = s.measure();
		if( weight && weight==weight )
		{
#ifdef NEW_CODE
			node_index_type nodeIndex = node->nodeData.nodeIndex;
#else // !NEW_CODE
			int nodeIndex = node->nodeData.nodeIndex;
#endif // NEW_CODE
			if( nodeIndex>=nodeToIndexMap.size() ) nodeToIndexMap.resize( nodeIndex+1 , -1 );
#ifdef NEW_CODE
			node_index_type idx = nodeToIndexMap[ nodeIndex ];
#else // !NEW_CODE
			int idx = nodeToIndexMap[ nodeIndex ];
#endif // NEW_CODE
			if( idx==-1 )
			{
#ifdef NEW_CODE
				idx = (node_index_type)simplices.size();
#else // !NEW_CODE
				idx = (int)simplices.size();
#endif // NEW_CODE
				nodeToIndexMap[ nodeIndex ] = idx;
				simplices.resize( idx+1 );
				simplices[idx].node = node;
			}
			simplices[idx].data.push_back( s );
		}
		return 1;
	}
	else
	{
#ifdef NEW_CODE
		size_t sCount = 0;
#else // !NEW_CODE
		int sCount = 0;
#endif // NEW_CODE
#ifdef THREAD_SAFE_CHILD_INIT
		if( !node->children ) node->template initChildren< ThreadSafe >( nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT

		// Split up the simplex and pass the parts on to the children
		Point< Real , Dim > center;
		Real width;
		node->centerAndWidth( center , width );

		std::vector< std::vector< Simplex< Real , Dim , Dim-1 > > > childSimplices( 1 );
		childSimplices[0].push_back( s );
		for( int d=0 ; d<Dim ; d++ )
		{
			Point< Real , Dim > n ; n[Dim-d-1] = 1;
			std::vector< std::vector< Simplex< Real , Dim , Dim-1 > > > temp( (int)( 1<<(d+1) ) );
			for( int c=0 ; c<(1<<d) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) childSimplices[c][i].split( n , center[Dim-d-1] , temp[2*c] , temp[2*c+1] );
			childSimplices = temp;
		}
#ifdef THREAD_SAFE_CHILD_INIT
		for( int c=0 ; c<(1<<Dim) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) sCount += _AddSimplex< ThreadSafe >( node->children+c , childSimplices[c][i] , maxDepth , simplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
		for( int c=0 ; c<(1<<Dim) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) sCount += _AddSimplex( node->children+c , childSimplices[c][i] , maxDepth , simplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
		return sCount;
	}
}

template< unsigned int Dim , class Real >
template< class Data , class _Data , bool Dual >
#ifdef NEW_CODE
size_t FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , ConstPointer( Data ) values , ConstPointer( int ) labels , int resolution[Dim] , std::vector< NodeSample< Dim , _Data > > derivatives[Dim] , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer , std::function< _Data ( const Data& ) > DataConverter )
#else // !NEW_CODE
int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , ConstPointer( Data ) values , ConstPointer( int ) labels , int resolution[Dim] , std::vector< NodeSample< Dim , _Data > > derivatives[Dim] , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer , std::function< _Data ( const Data& ) > DataConverter )
#endif // NEW_CODE
{
	auto Leaf = [&]( FEMTreeNode& root , const int idx[Dim] , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( idx[d]<0 || idx[d]>=(1<<maxDepth) ) return (FEMTreeNode*)NULL;
		FEMTreeNode* node = &root;
		for( int d=0 ; d<maxDepth ; d++ )
		{
#ifdef THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->template initChildren< false >( nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
			int cIndex = 0;
			for( int dd=0 ; dd<Dim ; dd++ ) if( idx[dd]&(1<<(maxDepth-d-1)) ) cIndex |= 1<<dd;
			node = node->children + cIndex;
		}
		return node;
	};
	auto FactorIndex = []( size_t i , const int resolution[Dim] , int idx[Dim] )
	{
		size_t ii = i;
		for( int d=0 ; d<Dim ; d++ ) idx[d] = ii % resolution[d] , ii /= resolution[d];
	};
	auto MakeIndex = [] ( const int idx[Dim] , const int resolution[Dim] )
	{
		size_t i = 0;
		for( int d=0 ; d<Dim ; d++ ) i = i * resolution[Dim-1-d] + idx[Dim-1-d];
		return i;
	};


	int maxResolution = resolution[0];
	for( int d=1 ; d<Dim ; d++ ) maxResolution = std::max< int >( maxResolution , resolution[d] );
	int maxDepth = 0;
	while( ( (1<<maxDepth) + ( Dual ? 0 : 1 ) )<maxResolution ) maxDepth++;

	size_t totalRes = 1;
	for( int d=0 ; d<Dim ; d++ ) totalRes *= resolution[d];

	// Iterate over each direction
	for( int d=0 ; d<Dim ; d++ ) for( size_t i=0 ; i<totalRes ; i++ )
	{
		// Factor the index into directional components and get the index of the next cell
		int idx[Dim] ; FactorIndex( i , resolution , idx ) ; idx[d]++;

		if( idx[d]<resolution[d] )
		{
			// Get the index of the next cell
			size_t ii = MakeIndex( idx , resolution );

			// [NOTE] There are no derivatives across negative labels
			if( labels[i]!=labels[ii] && labels[i]>=0 && labels[ii]>=0 )
			{
				if( !Dual ) idx[d]--;
				NodeSample< Dim , _Data > nodeSample;
				nodeSample.node = Leaf( root , idx , maxDepth );
				nodeSample.data = DataConverter( values[ii] ) - DataConverter( values[i] );
				if( nodeSample.node ) derivatives[d].push_back( nodeSample );
			}
		}
	}
	return maxDepth;
}

template< unsigned int Dim , class Real >
template< bool Dual , class Data >
unsigned int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , DerivativeStream< Data >& dStream , std::vector< NodeSample< Dim , Data > > derivatives[Dim] , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	// Note:
	// --   Dual: The difference between [i] and [i+1] is stored at cell [i+1]
	// -- Primal: The difference between [i] and [i+1] is stored at cell [i]

	// Find the leaf containing the specified cell index
	auto Leaf = [&]( FEMTreeNode& root , const unsigned int idx[Dim] , unsigned int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( idx[d]<0 || idx[d]>=(unsigned int)(1<<maxDepth) ) return (FEMTreeNode*)NULL;
		FEMTreeNode* node = &root;
		for( unsigned int d=0 ; d<maxDepth ; d++ )
		{
#ifdef THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->template initChildren< false >( nodeAllocator , NodeInitializer );
#else // !THREAD_SAFE_CHILD_INIT
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
#endif // THREAD_SAFE_CHILD_INIT
			int cIndex = 0;
			for( int dd=0 ; dd<Dim ; dd++ ) if( idx[dd]&(1<<(maxDepth-d-1)) ) cIndex |= 1<<dd;
			node = node->children + cIndex;
		}
		return node;
	};

	unsigned int resolution[Dim];
	dStream.resolution( resolution );
	unsigned int maxResolution = resolution[0];
	for( int d=1 ; d<Dim ; d++ ) maxResolution = std::max< unsigned int >( maxResolution , resolution[d] );
	unsigned int maxDepth = 0;

	// If we are using a dual formulation, we need at least maxResolution cells.
	// Otherwise, we need at least maxResolution-1 cells.
	while( (unsigned int)( (1<<maxDepth) + ( Dual ? 0 : 1 ) )<maxResolution ) maxDepth++;

	unsigned int idx[Dim] , dir;
	Data dValue;
	while( dStream.nextDerivative( idx , dir , dValue ) )
	{
		if( Dual ) idx[dir]++;
		NodeSample< Dim , Data > nodeSample;
		nodeSample.node = Leaf( root , idx , maxDepth );
		nodeSample.data = dValue;
		if( nodeSample.node ) derivatives[dir].push_back( nodeSample );
	}
	return maxDepth;
}
