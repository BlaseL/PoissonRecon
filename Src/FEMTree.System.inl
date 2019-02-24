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

///////////////////////////////////
// BaseFEMIntegrator::Constraint //
///////////////////////////////////
template< unsigned int ... TDegrees , unsigned int ... CDegrees , unsigned int CDim >
template< bool IterateFirst >
void BaseFEMIntegrator::Constraint< UIntPack< TDegrees ... > , UIntPack< CDegrees ... > , CDim >::setStencil( CCStencil & stencil ) const
{
	static const int Dim = sizeof ... ( TDegrees );
	int center = ( 1<<_highDepth )>>1;
	int femOffset[Dim] , cOffset[Dim];
	static const int overlapStart[] = { ( IterateFirst ? BSplineOverlapSizes< CDegrees , TDegrees >::OverlapStart : BSplineOverlapSizes< TDegrees , CDegrees >::OverlapStart ) ... };
	if( IterateFirst )
	{
		for( int d=0 ; d<Dim ; d++ ) cOffset[d] = center;
		WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , UIntPack< BSplineOverlapSizes< TDegrees , CDegrees >::OverlapSize ... >() , [&]( int d , int i ){ femOffset[d] = i + center + overlapStart[d]; } , [&]( Point< double , CDim >& p ){ p = ccIntegrate( femOffset , cOffset ); } , stencil() );
	}
	else
	{
		for( int d=0 ; d<Dim ; d++ ) femOffset[d] = center;
		WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , UIntPack< BSplineOverlapSizes< TDegrees , CDegrees >::OverlapSize ... >() , [&]( int d , int i ){   cOffset[d] = i + center + overlapStart[d]; } , [&]( Point< double , CDim >& p ){ p = ccIntegrate( femOffset , cOffset );} , stencil() );
	}
}
template< unsigned int ... TDegrees , unsigned int ... CDegrees , unsigned int CDim >
template< bool IterateFirst >
void BaseFEMIntegrator::Constraint< UIntPack< TDegrees ... > , UIntPack< CDegrees ... > , CDim >::setStencils( PCStencils& stencils ) const
{
	static const int Dim = sizeof ... ( TDegrees );
	typedef UIntPack< BSplineOverlapSizes< TDegrees, CDegrees >::OverlapSize ... > OverlapSizes;
	// [NOTE] We want the center to be at the first node of the brood, which is not the case when childDepth is 1.
	int center = ( 1<<_highDepth )>>1 ; center = ( center>>1 )<<1;	
	int fineCenter[Dim] , femOffset[Dim] , cOffset[Dim];
	static const int overlapStart[] = { ( IterateFirst ? BSplineOverlapSizes< CDegrees , TDegrees >::OverlapStart : BSplineOverlapSizes< TDegrees , CDegrees >::OverlapStart ) ... };
	std::function< void ( int , int )               > outerUpdateState = [&]( int d , int i ){ fineCenter[Dim-d-1] = i+center; };
	std::function< void ( Point< double , CDim >& ) > innerFunction    = [&]( Point< double , CDim >& p ){ p = pcIntegrate( femOffset , cOffset ); };
	std::function< void ( int , int )               > innerUpdateState = [&]( int d , int i ){ femOffset[d] = IterateFirst ? (i+center/2+overlapStart[d]) : center/2 , cOffset[d] = IterateFirst ? fineCenter[d] : (i+fineCenter[d]+overlapStart[d]); };
	std::function< void ( CCStencil& )              > outerFunction    = [&]( CCStencil& s )
	{
		WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , OverlapSizes() , innerUpdateState , innerFunction , s() );
	};
	WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , IsotropicUIntPack< Dim , 2 >() , outerUpdateState , outerFunction , stencils() );
}
template< unsigned int ... TDegrees , unsigned int ... CDegrees , unsigned int CDim >
template< bool IterateFirst >
void BaseFEMIntegrator::Constraint< UIntPack< TDegrees ... > , UIntPack< CDegrees ... > , CDim >::setStencils( CPStencils& stencils ) const
{
	static const int Dim = sizeof ... ( TDegrees );
	typedef UIntPack< BSplineOverlapSizes< TDegrees , CDegrees >::OverlapSize ... > OverlapSizes;
	// [NOTE] We want the center to be at the first node of the brood, which is not the case when childDepth is 1.
	int center = ( 1<<_highDepth )>>1 ; center = ( center>>1 )<<1;
	static const int overlapStart[] = { ( IterateFirst ? BSplineOverlapSizes< CDegrees , TDegrees >::OverlapStart : BSplineOverlapSizes< TDegrees , CDegrees >::OverlapStart ) ... };
	int fineCenter[Dim] , femOffset[Dim] , cOffset[Dim];
	std::function< void ( int , int )               > outerUpdateState = [&]( int d , int i ){ fineCenter[Dim-d-1] = i+center; };
	std::function< void ( Point< double , CDim >& ) > innerFunction    = [&]( Point< double , CDim >& p ){ p = cpIntegrate( femOffset , cOffset ); };
	std::function< void ( int , int )               > innerUpdateState = [&]( int d , int i ){ femOffset[d] = IterateFirst ? (i+fineCenter[d]+overlapStart[d]) : fineCenter[d] , cOffset[d] = IterateFirst ? center/2 : (i+center/2+overlapStart[d]); };
	std::function< void ( CCStencil& )              > outerFunction    = [&]( CCStencil& s )
	{
		WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , OverlapSizes() , innerUpdateState , innerFunction , s() );
	};
	WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , IsotropicUIntPack< Dim , 2 >() , outerUpdateState , outerFunction , stencils() );
}

///////////////////////////////
// BaseFEMIntegrator::System //
///////////////////////////////
template< unsigned int ... TDegrees >
template< bool IterateFirst >
void BaseFEMIntegrator::System< UIntPack< TDegrees ... > >::setStencil( CCStencil & stencil ) const
{
	static const int Dim = sizeof ... ( TDegrees );
	int center = ( 1<<_highDepth )>>1;
	int offset1[Dim] , offset2[Dim];
	static const int overlapStart[] = { BSplineOverlapSizes< TDegrees , TDegrees >::OverlapStart ... };
	if( IterateFirst )
	{
		for( int d=0 ; d<Dim ; d++ ) offset2[d] = center;
		WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , UIntPack< BSplineOverlapSizes< TDegrees , TDegrees >::OverlapSize ... >() , [&]( int d , int i ){ offset1[d] = i + center + overlapStart[d]; } , [&]( double& v ){ v = ccIntegrate( offset1 , offset2 ); } , stencil() );
	}
	else
	{
		for( int d=0 ; d<Dim ; d++ ) offset1[d] = center;
		WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , UIntPack< BSplineOverlapSizes< TDegrees , TDegrees >::OverlapSize ... >() , [&]( int d , int i ){ offset2[d] = i + center + overlapStart[d]; } , [&]( double& v ){ v = ccIntegrate( offset1 , offset2 ); } , stencil() );
	}
}
template< unsigned int ... TDegrees >
template< bool IterateFirst >
void BaseFEMIntegrator::System< UIntPack< TDegrees ... > >::setStencils( PCStencils& stencils ) const
{
	static const int Dim = sizeof ... ( TDegrees );
	typedef UIntPack< BSplineOverlapSizes< TDegrees , TDegrees >::OverlapSize ... > OverlapSizes;
	// [NOTE] We want the center to be at the first node of the brood
	// Which is not the case when childDepth is 1.
	int center = ( 1<<_highDepth )>>1 ; center = ( center>>1 )<<1;
	static const int overlapStart[] = { BSplineOverlapSizes< TDegrees , TDegrees >::OverlapStart ... };
	int fineCenter[Dim] , offset1[Dim] , offset2[Dim];
	std::function< void ( int , int )  > outerUpdateState = [&]( int d , int i ){ fineCenter[Dim-d-1] = i+center; };
	std::function< void ( double& )    > innerFunction    = [&]( double& v ){ v = pcIntegrate( offset1 , offset2 ); };
	std::function< void ( int , int )  > innerUpdateState = [&]( int d , int i ){ offset1[d] = IterateFirst ? (i+center/2+overlapStart[d]) : center/2 , offset2[d] = IterateFirst ? fineCenter[d] : (i+fineCenter[d]+overlapStart[d]); };
	std::function< void ( CCStencil& ) > outerFunction    = [&]( CCStencil& s )
	{
		WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , OverlapSizes() , innerUpdateState , innerFunction , s() );
	};
	WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , IsotropicUIntPack< Dim , 2 >() , outerUpdateState , outerFunction , stencils() );
}
/////////////////////////////////
// BaseFEMIntegrator::UpSample //
/////////////////////////////////
template< unsigned int ... TDegrees >
void BaseFEMIntegrator::RestrictionProlongation< UIntPack< TDegrees ... > >::setStencil( UpSampleStencil & stencil ) const
{
	static const int Dim = sizeof ... ( TDegrees );
	int highCenter = ( 1<<_highDepth )>>1;
	int pOff[Dim] , cOff[Dim];
	static const int upSampleStart[] = { BSplineSupportSizes< TDegrees >::UpSampleStart ... };
	for( int d=0 ; d<Dim ; d++ ) pOff[d] = highCenter/2;
	WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , UIntPack< BSplineSupportSizes< TDegrees >::UpSampleSize ... >() , [&]( int d , int i ){ cOff[d] = i + highCenter + upSampleStart[d]; } , [&]( double& v ){ v = upSampleCoefficient( pOff , cOff ); } , stencil() );
}
template< unsigned int ... TDegrees >
void BaseFEMIntegrator::RestrictionProlongation< UIntPack< TDegrees ... > >::setStencils( DownSampleStencils& stencils ) const
{
	static const int Dim = sizeof ... ( TDegrees );
	// [NOTE] We want the center to be at the first node of the brood, which is not the case when childDepth is 1.
	int highCenter = ( 1<<_highDepth )>>1 ; highCenter = ( highCenter>>1 )<<1;	
	int pOff[Dim] , cOff[Dim];
	static const int offsets[] = { BSplineSupportSizes< TDegrees >::DownSample0Start ... };
	std::function< void ( double& )            > innerFunction    = [&]( double& v ){ v = upSampleCoefficient( pOff , cOff ); };
	std::function< void ( int , int )          > innerUpdateState = [&]( int d , int i ){ pOff[d] = cOff[d]/2 + i + offsets[d]; };
	std::function< void ( int , int )          > outerUpdateState = [&]( int d , int i ){ cOff[Dim-d-1] = i+highCenter; };
	std::function< void ( DownSampleStencil& ) > outerFunction    = [&]( DownSampleStencil& s )
	{
		WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , UIntPack< ( - BSplineSupportSizes< TDegrees >::DownSample0Start + BSplineSupportSizes< TDegrees >::DownSample1End + 1 ) ... >() , innerUpdateState , innerFunction , s() );
	};
	WindowLoop< Dim >::Run( IsotropicUIntPack< Dim , 0 >() , IsotropicUIntPack< Dim , 2 >() , outerUpdateState , outerFunction , stencils() );
}

///////////////////////////////
// FEMIntegrator::Constraint //
///////////////////////////////

template< unsigned int ... TSignatures , unsigned int ... TDerivatives , unsigned int ... CSignatures , unsigned int ... CDerivatives , unsigned int CDim >
Point< double , CDim > FEMIntegrator::Constraint< UIntPack< TSignatures ... > , UIntPack< TDerivatives ... > , UIntPack< CSignatures ... > , UIntPack< CDerivatives ... > , CDim >::_integrate( IntegrationType iType , const int off1[] , const int off2[] ) const
{
	Point< double , CDim > integral;
	for( unsigned int i=0 ; i<_weightedIndices.size() ; i++ )
	{
		const _WeightedIndices& w = _weightedIndices[i];
		unsigned int _d1[Dim] , _d2[Dim];
		TFactorDerivatives( w.d1 , _d1 );
		CFactorDerivatives( w.d2 , _d2 );
		double __integral = _integral( iType , off1 , off2 , _d1 , _d2 );
		for( unsigned int j=0 ; j<w.indices.size() ; j++ ) integral[ w.indices[j].first ] += w.indices[j].second * __integral;
	}
	return integral;
}

#ifndef MOD
#define MOD( a , b ) ( (a)>0 ? (a) % (b) : ( (b) - ( -(a) % (b) ) ) % (b) )
#endif // MOD

/////////////
// FEMTree //
/////////////
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
#ifdef NEW_CODE
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::setMultiColorIndices( UIntPack< FEMSigs ... > , ThreadPool &tp , int depth , std::vector< std::vector< size_t > >& indices ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::setMultiColorIndices( UIntPack< FEMSigs ... > , int depth , std::vector< std::vector< size_t > >& indices ) const
#endif // NEW_THREADS
#else // !NEW_CODE
void FEMTree< Dim , Real >::setMultiColorIndices( UIntPack< FEMSigs ... > , int depth , std::vector< std::vector< int > >& indices ) const
#endif // NEW_CODE
{
#ifdef NEW_THREADS
	_setMultiColorIndices( UIntPack< FEMSigs ... >() , tp , _sNodesBegin(depth) , _sNodesEnd(depth) , indices );
#else // !NEW_THREADS
	_setMultiColorIndices( UIntPack< FEMSigs ... >() , _sNodesBegin(depth) , _sNodesEnd(depth) , indices );
#endif // NEW_THREADS
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
#ifdef NEW_CODE
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_setMultiColorIndices( UIntPack< FEMSigs ... > , ThreadPool &tp , node_index_type start , node_index_type end , std::vector< std::vector< size_t > >& indices ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_setMultiColorIndices( UIntPack< FEMSigs ... > , node_index_type start , node_index_type end , std::vector< std::vector< size_t > >& indices ) const
#endif // NEW_THREADS
#else // !NEW_CODE
void FEMTree< Dim , Real >::_setMultiColorIndices( UIntPack< FEMSigs ... > , int start , int end , std::vector< std::vector< int > >& indices ) const
#endif // NEW_CODE
{
	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	typedef UIntPack< ( 1 - BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree , FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > Moduli;
	static const unsigned int Colors = WindowSize< Moduli >::Size;
	indices.resize( Colors );
#ifdef NEW_CODE
#ifdef NEW_THREADS
	struct ColorCount
	{
		size_t count[ Colors ];
		ColorCount( void ){ memset( count , 0 , sizeof(count) ); }
	};
	std::vector< ColorCount > counts( tp.threadNum() );
#endif // NEW_THREADS
	size_t count[ Colors ];
	memset( count , 0 , sizeof(count) );
#else // !NEW_CODE
	int count[ Colors ];
	memset( count , 0 , sizeof(count) );
#endif // NEW_CODE
	auto MCIndex = [&] ( const FEMTreeNode* node )
	{
		LocalDepth d ; LocalOffset off ; _localDepthAndOffset( node , d , off );
		int index = 0;
		for( int dd=0 ; dd<Dim ; dd++ ) index = index * Moduli::Values[Dim-dd-1] + MOD( off[Dim-dd-1] , Moduli::Values[Dim-dd-1] );
		return index;
	};
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( start , end , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( start , end , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=start ; i<end ; i++ )
#else // !NEW_CODE
	for( int i=start ; i<end ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
		{
			int idx = MCIndex( _sNodes.treeNodes[i] );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			counts[thread()].count[idx]++;
#else // !NEW_THREAD_NUM
			counts[thread].count[idx]++;
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp atomic
			count[idx]++;
#endif // NEW_THREADS
		}
	}
#ifdef NEW_THREADS
	);
	for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) for( int i=0 ; i<Colors ; i++) count[i] += counts[t].count[i];
#endif // NEW_THREADS

	for( int i=0 ; i<Colors ; i++ ) indices[i].reserve( count[i] ) , count[i]=0;

#ifdef NEW_CODE
	for( node_index_type i=start ; i<end ; i++ ) if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
#else // !NEW_CODE
	for( int i=start ; i<end ; i++ ) if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
#endif // NEW_CODE
	{
		int idx = MCIndex( _sNodes.treeNodes[i] );
		indices[idx].push_back( i - start );
	}
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , typename TDotT , typename SORWeights , unsigned int ... PointDs >
#ifdef NEW_THREADS
int FEMTree< Dim , Real >::_solveFullSystemGS( UIntPack< FEMSigs ... > , ThreadPool &tp , const typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , Pointer( T ) solution , ConstPointer( T ) prolongedSolution , ConstPointer( T ) constraints , TDotT Dot , int iters , bool coarseToFine , SORWeights sorWeights , _SolverStats& stats , bool computeNorms , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
int FEMTree< Dim , Real >::_solveFullSystemGS( UIntPack< FEMSigs ... > , const typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , Pointer( T ) solution , ConstPointer( T ) prolongedSolution , ConstPointer( T ) constraints , TDotT Dot , int iters , bool coarseToFine , SORWeights sorWeights , _SolverStats& stats , bool computeNorms , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_THREADS
{
	double& systemTime = stats.systemTime;
	double&  solveTime = stats. solveTime;
	systemTime = solveTime = 0.;

	CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > > ccStencil;
	PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > > pcStencils;
	F.template setStencil< false >( ccStencil );
	F.template setStencils< true >( pcStencils );
	double bNorm=0 , inRNorm=0 , outRNorm=0;
	if( depth>=0 )
	{
#ifdef NEW_CODE
		SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > M;
#else // !NEW_CODE
		SparseMatrix< Real , int , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > M;
#endif // NEW_CODE
		double t = Time();
		Pointer( Real ) D = AllocPointer< Real >( _sNodesEnd( depth ) - _sNodesBegin( depth ) );
		Pointer( T ) _constraints = AllocPointer< T >( _sNodesSize( depth ) );
#ifdef NEW_THREADS
		_getSliceMatrixAndProlongationConstraints( UIntPack< FEMSigs ... >() , tp , F , M , D , bsData , depth , _sNodesBegin( depth ) , _sNodesEnd( depth ) , prolongedSolution , _constraints , ccStencil , pcStencils , interpolationInfo... );
#else // !NEW_THREADS
		_getSliceMatrixAndProlongationConstraints( UIntPack< FEMSigs ... >() , F , M , D , bsData , depth , _sNodesBegin( depth ) , _sNodesEnd( depth ) , prolongedSolution , _constraints , ccStencil , pcStencils , interpolationInfo... );
#endif // NEW_THREADS
#ifdef NEW_CODE
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
		tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd(depth) , [&]( const ThreadPool::ThreadNum & , size_t i ){ _constraints[ i - _sNodesBegin(depth) ] = constraints[ _sNodes.treeNodes[i]->nodeData.nodeIndex ] - _constraints[ i - _sNodesBegin(depth) ]; } );
#else // !NEW_THREAD_NUM
		tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd(depth) , [&]( unsigned int , size_t i ){ _constraints[ i - _sNodesBegin(depth) ] = constraints[ _sNodes.treeNodes[i]->nodeData.nodeIndex ] - _constraints[ i - _sNodesBegin(depth) ]; } );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
		for( node_index_type i=_sNodesBegin(depth) ; i<_sNodesEnd(depth) ; i++ ) _constraints[ i - _sNodesBegin(depth) ] = constraints[ _sNodes.treeNodes[i]->nodeData.nodeIndex ] - _constraints[ i - _sNodesBegin(depth) ];
#endif // NEW_THREADS
		{
			node_index_type begin = _sNodesBegin( depth ) , end = _sNodesEnd( depth );
			for( node_index_type i=begin ; i<end ; i++ ) if( M.rowSize( i-begin ) ) D[i-begin] *= sorWeights[i];
		}
#else // !NEW_CODE
#pragma omp parallel for
		for( int i=_sNodesBegin(depth) ; i<_sNodesEnd(depth) ; i++ ) _constraints[ i - _sNodesBegin(depth) ] = constraints[ _sNodes.treeNodes[i]->nodeData.nodeIndex ] - _constraints[ i - _sNodesBegin(depth) ];
		{
			int begin = _sNodesBegin( depth ) , end = _sNodesEnd( depth );
			for( int i=begin ; i<end ; i++ ) if( M.rowSize( i-begin ) ) D[i-begin] *= sorWeights[i];
		}
#endif // NEW_CODE

		systemTime += Time()-t;
		// The list of multi-colored indices  for each in-memory slice
#ifdef NEW_CODE
		std::vector< std::vector< size_t > > mcIndices;
#else // !NEW_CODE
		std::vector< std::vector< int > > mcIndices;
#endif // NEW_CODE
#ifdef NEW_THREADS
		_setMultiColorIndices( UIntPack< FEMSigs ... >() , tp , _sNodesBegin( depth ) , _sNodesEnd( depth ) , mcIndices );
#else // !NEW_THREADS
		_setMultiColorIndices( UIntPack< FEMSigs ... >() , _sNodesBegin( depth ) , _sNodesEnd( depth ) , mcIndices );
#endif // NEW_THREADS

		ConstPointer( T ) B = _constraints;
		Pointer( T ) X = GetPointer( &solution[0] + _sNodesBegin( depth ) , _sNodesSize( depth ) );
		if( computeNorms )
		{
#ifdef NEW_THREADS
			std::vector< double > bNorms( tp.threadNum() , 0 ) , inRNorms( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
			tp.parallel_for( 0 , M.rows() , [&]( const ThreadPool::ThreadNum &thread , size_t j )
#else // !NEW_THREAD_NUM
			tp.parallel_for( 0 , M.rows() , [&]( unsigned int thread , size_t j )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : bNorm , inRNorm )
#ifdef NEW_CODE
			for( matrix_index_type j=0 ; j<(matrix_index_type)M.rows() ; j++ )
#else // !NEW_CODE
			for( int j=0 ; j<M.rows() ; j++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
				T temp = {};
#ifdef NEW_CODE_SPARSE_MATRIX
				ConstPointer( MatrixEntry< Real , matrix_index_type > ) start = M[j];
				ConstPointer( MatrixEntry< Real , matrix_index_type > ) end = start + M.rowSize(j);
				ConstPointer( MatrixEntry< Real , matrix_index_type > ) e;
#else // !NEW_CODE_SPARSE_MATRIX
				ConstPointer( MatrixEntry< Real > ) start = M[j];
				ConstPointer( MatrixEntry< Real > ) end = start + (unsigned long long)M.rowSize(j);
				ConstPointer( MatrixEntry< Real > ) e;
#endif // NEW_CODE_SPARSE_MATRIX
				for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
				bNorms[thread()] += Dot( B[j] , B[j] );
				inRNorms[thread()] += Dot( temp - B[j] , temp - B[j] );
#else // !NEW_THREAD_NUM
				bNorms[thread] += Dot( B[j] , B[j] );
				inRNorms[thread] += Dot( temp - B[j] , temp - B[j] );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
				bNorm += Dot( B[j] , B[j] );
				inRNorm += Dot( temp - B[j] , temp - B[j] );
#endif // NEW_THREADS
			}
#ifdef NEW_THREADS
			);
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) bNorm += bNorms[t] , inRNorm += inRNorms[t];
#endif // NEW_THREADS
		}

		t = Time();
		MemoryUsage();
#ifdef NEW_THREADS
		for( int i=0 ; i<iters ; i++ ) M.gsIteration( tp , mcIndices , ( ConstPointer( Real ) )D , B , X , coarseToFine , true );
#else // !NEW_THREADS
		for( int i=0 ; i<iters ; i++ ) M.gsIteration( mcIndices , ( ConstPointer( Real ) )D , B , X , coarseToFine , true );
#endif // NEW_THREADS
		FreePointer( D );
		solveTime += Time() - t;

		if( computeNorms )
		{
#ifdef NEW_THREADS
			std::vector< double > outRNorms( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
			tp.parallel_for( 0 , M.rows() , [&]( const ThreadPool::ThreadNum &thread , size_t j )
#else // !NEW_THREAD_NUM
			tp.parallel_for( 0 , M.rows() , [&]( unsigned int thread , size_t j )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : outRNorm )
#ifdef NEW_CODE
			for( matrix_index_type j=0 ; j<(matrix_index_type)M.rows() ; j++ )
#else // !NEW_CODE
			for( int j=0 ; j<M.rows() ; j++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
				T temp = {};
#ifdef NEW_CODE_SPARSE_MATRIX
				ConstPointer( MatrixEntry< Real , matrix_index_type > ) start = M[j];
				ConstPointer( MatrixEntry< Real , matrix_index_type > ) end = start + M.rowSize(j);
				ConstPointer( MatrixEntry< Real , matrix_index_type > ) e;
#else // !NEW_CODE_SPARSE_MATRIX
				ConstPointer( MatrixEntry< Real > ) start = M[j];
				ConstPointer( MatrixEntry< Real > ) end = start + (unsigned long long)M.rowSize(j);
				ConstPointer( MatrixEntry< Real > ) e;
#endif // NEW_CODE_SPARSE_MATRIX
				for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
				outRNorms[thread()] += Dot( temp-B[j] , temp-B[j] );
#else // !NEW_THREAD_NUM
				outRNorms[thread] += Dot( temp-B[j] , temp-B[j] );
#endif //NEW_THREAD_NUM
#else // !NEW_THREADS
				outRNorm += Dot( temp-B[j] , temp-B[j] );
#endif // NEW_THREADS
			}
#ifdef NEW_THREADS
			);
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) outRNorm += outRNorms[t];
#endif // NEW_THREADS
		}
		FreePointer( _constraints );
	}
	if( computeNorms ) stats.bNorm2 = bNorm , stats.inRNorm2 = inRNorm , stats.outRNorm2 = outRNorm;
	MemoryUsage();

	return iters;
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , typename TDotT , typename SORWeights , unsigned int ... PointDs >
#ifdef NEW_THREADS
int FEMTree< Dim , Real >::_solveSlicedSystemGS( UIntPack< FEMSigs ... > , ThreadPool &tp , const typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , Pointer( T ) solution , ConstPointer( T ) prolongedSolution , ConstPointer( T ) constraints , TDotT Dot , int iters , bool coarseToFine , unsigned int sliceBlockSize , SORWeights sorWeights , _SolverStats& stats , bool computeNorms , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
int FEMTree< Dim , Real >::_solveSlicedSystemGS( UIntPack< FEMSigs ... > , const typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , Pointer( T ) solution , ConstPointer( T ) prolongedSolution , ConstPointer( T ) constraints , TDotT Dot , int iters , bool coarseToFine , unsigned int sliceBlockSize , SORWeights sorWeights , _SolverStats& stats , bool computeNorms , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_THREADS
{
#ifdef NEW_THREADS
	if( sliceBlockSize<=0 ) return _solveFullSystemGS( UIntPack< FEMSigs ... >() , tp , F , bsData , depth , solution , prolongedSolution , constraints , Dot , iters , coarseToFine , sorWeights , stats , computeNorms , interpolationInfo ... );
#else // !NEW_THREADS
	if( sliceBlockSize<=0 ) return _solveFullSystemGS( UIntPack< FEMSigs ... >() , F , bsData , depth , solution , prolongedSolution , constraints , Dot , iters , coarseToFine , sorWeights , stats , computeNorms , interpolationInfo ... );
#endif // NEW_THREADS
	CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > > ccStencil;
	PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > > pcStencils;
	F.template setStencil< false >( ccStencil );
	F.template setStencils< true >( pcStencils );

	{
		// Assuming Degree=2 and we are solving forward using two iterations, the pattern of relaxations should look like:
		//      +--+--+--+--+--+
		//      *  |  |  |  |  |
		//     o|  |  |  |  |  |
		//    o |  |  |  |  |  |
		//   o  |  |  |  |  |  |
		//  o   |  |  |  |  |  |
		// o    |  |  |  |  |  |
		//      |  *  |  |  |  |
		//      | *|  |  |  |  |
		//      |* |  |  |  |  |
		//      *  |  |  |  |  |
		//     o|  |  |  |  |  |
		//    o |  |  |  |  |  |
		//      |  |  *  |  |  |
		//      |  | *|  |  |  |
		//      |  |* |  |  |  |
		//      |  *  |  |  |  |
		//      | *|  |  |  |  |
		//      |* |  |  |  |  |
		//      |  |  |  *  |  |
		//      |  |  | *|  |  |
		//      |  |  |* |  |  |
		//      |  |  *  |  |  |
		//      |  | *|  |  |  |
		//      |  |* |  |  |  |
		//      |  |  |  |  *  |
		//      |  |  |  | *|  |
		//      |  |  |  |* |  |
		//      |  |  |  *  |  |
		//      |  |  | *|  |  |
		//      |  |  |* |  |  |
		//      |  |  *  |  |  |
		//      |  |  |  |  |  *
		//      |  |  |  |  | *|
		//      |  |  |  |  |* |
		//      |  |  |  |  *  |
		//      |  |  |  | *|  |
		//      |  |  |  |* |  |
		//      |  |  |  |  |  |  o
		//      |  |  |  |  |  | o
		//      |  |  |  |  |  |o
		//      |  |  |  |  |  *
		//      |  |  |  |  | *|
		//      |  |  |  |  |* |

		const int SliceBlockSize = (int)sliceBlockSize;
		// OverlapRadius = Degree
		const int OverlapRadii[] = { ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... };
		const int OverlapBlockRadius = ( OverlapRadii[Dim-1] + SliceBlockSize - 1 ) / SliceBlockSize;
		static const int LastFEMSig = UIntPack< FEMSigs ... >::template Get< Dim-1 >();
		int _sliceBegin = _BSplineBegin< LastFEMSig >( depth ) , _sliceEnd = _BSplineEnd< LastFEMSig >( depth );

		int blockBegin = ( _sliceBegin - ( SliceBlockSize - 1 ) ) / SliceBlockSize , blockEnd = ( _sliceEnd + ( SliceBlockSize - 1 ) ) / SliceBlockSize;
		std::function< int ( int ) > BlockFirst = [&]( int b ){ return std::max< int >( b * SliceBlockSize , _sliceBegin ); };
		std::function< int ( int ) > BlockLast  = [&]( int b ){ return std::min< int >( b * SliceBlockSize + SliceBlockSize - 1 , _sliceEnd - 1 ); };

		auto BBlock = [&]( int d , int b , ConstPointer( T ) B )
		{
			return GetPointer( &B[0] + _sNodesBegin( d , BlockFirst( b ) ) , _sNodesEnd( d , BlockLast( b ) ) - _sNodesBegin( d , BlockFirst( b ) ) );
		};
		auto XBlocks = [&]( int d , int b , Pointer( T ) X )
		{
			return GetPointer( &X[0] + _sNodesBegin( d , BlockFirst( b ) ) , _sNodesBegin( d , BlockFirst( b - OverlapBlockRadius ) ) - _sNodesBegin( d , BlockFirst( b ) ) , _sNodesEnd( d , BlockLast( b + OverlapBlockRadius ) ) - _sNodesBegin( d , BlockFirst( b ) ) );
		};

		double& systemTime = stats.systemTime;
		double&  solveTime = stats. solveTime;
		systemTime = solveTime = 0.;

		struct BlockWindow
		{
		protected:
			int _begin , _end;
		public:
			BlockWindow( int begin , int end )
			{
				if( begin<=end ) _begin = begin , _end = end;
				else             _begin = end+1 , _end = begin+1;
			}
			int size( void ) const { return _end-_begin; }
			BlockWindow& operator += ( int off ){ _begin += off , _end += off ; return *this; }
			BlockWindow& operator -= ( int off ){ _begin -= off , _end -= off ; return *this; }
			BlockWindow& operator++ ( void ){ _begin++ , _end++ ; return *this; }
			BlockWindow& operator-- ( void ){ _begin-- , _end-- ; return *this; }
			int begin( bool forward ) const { return forward ? _begin : _end-1; }
			int end  ( bool forward ) const { return forward ? _end : _begin-1; }
			bool inBlock( int b ) const { return b>=_begin && b<_end; }
		};
		double bNorm=0 , inRNorm=0 , outRNorm=0;
		bool forward = !coarseToFine;
		int residualOffset = computeNorms ? OverlapBlockRadius : 0;
		// Set the number of in-memory blocks required for a temporally blocked solver
		const int ColorModulus = OverlapBlockRadius;
		// The number of in-core blocks over which we relax
		// [WARNING] If the block size is larger than one, we may be able to use fewer blocks
		int solveBlocks = std::max< int >( 0 , std::min< int >( ColorModulus*iters - ( ColorModulus-1 ) , blockEnd-blockBegin ) );
		// The number of in-core blocks over which we either solve or compute residuals
		int matrixBlocks = std::max< int >( 1 , std::min< int >( solveBlocks+2*residualOffset , blockEnd-blockBegin ) );
		// The list of matrices for each in-memory block
#ifdef NEW_CODE
		Pointer( SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > ) _M = NewPointer< SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > >( matrixBlocks );
#else // !NEW_CODE
		Pointer( SparseMatrix< Real , int , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > ) _M = NewPointer< SparseMatrix< Real , int , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > >( matrixBlocks );
#endif // NEW_CODE
		Pointer( Pointer( Real ) ) _D = AllocPointer< Pointer( Real ) >( matrixBlocks );
		std::vector< Pointer( T ) > _constraints( matrixBlocks );
		for( int i=0 ; i<matrixBlocks ; i++ ) _D[i] = NullPointer( Real ) , _constraints[i] = NullPointer( T );
		// The list of multi-colored indices  for each in-memory block
#ifdef NEW_CODE
		Pointer( std::vector< std::vector< size_t > > ) mcIndices = NewPointer< std::vector< std::vector< size_t > > >( solveBlocks );
#else // !NEW_CODE
		Pointer( std::vector< std::vector< int > > ) mcIndices = NewPointer< std::vector< std::vector< int > > >( solveBlocks );
#endif // NEW_CODE
		int dir = forward ? 1 : -1 , start = forward ? blockBegin : blockEnd-1 , end = forward ? blockEnd : blockBegin-1;
		const BlockWindow FullWindow( blockBegin , blockEnd );
		BlockWindow residualWindow( FullWindow.begin(forward) , FullWindow.begin(forward) - ( ColorModulus*iters - ( ColorModulus-1 ) ) * dir - 2*residualOffset*dir );
		BlockWindow solveWindow( FullWindow.begin(forward) - residualOffset*dir , FullWindow.begin(forward) - residualOffset*dir - ( ColorModulus*iters - ( ColorModulus-1 ) ) * dir );
		// If we are solving forward we start in a block S with S mod ColorModulus = ColorModulus-1
		// and end in a block E with E mod ColorModulus = 0
		while( MOD( solveWindow.begin(!forward) , ColorModulus )!=( forward ? ColorModulus-1 : 0 ) ) solveWindow -= dir , residualWindow -= dir;
#ifdef NEW_CODE
		size_t maxBlockSize = 0;
#else // !NEW_CODE
		int maxBlockSize = 0;
#endif // NEW_CODE
		BlockWindow _residualWindow = residualWindow;
		for( ; _residualWindow.end(!forward)*dir<FullWindow.end(forward)*dir ; _residualWindow += dir )
		{
			int b = _residualWindow.begin(!forward);
#ifdef NEW_CODE
			if( FullWindow.inBlock( b ) ) maxBlockSize = std::max< size_t >( maxBlockSize , _sNodesEnd( depth , BlockLast( b ) ) - _sNodesBegin( depth , BlockFirst( b ) ) );
#else // !NEW_CODE
			if( FullWindow.inBlock( b ) ) maxBlockSize = std::max< int >( maxBlockSize , _sNodesEnd( depth , BlockLast( b ) ) - _sNodesBegin( depth , BlockFirst( b ) ) );
#endif // NEW_CODE
		}
#ifdef NEW_CODE
		if( maxBlockSize>std::numeric_limits< matrix_index_type >::max() ) ERROR_OUT( "more entries in a block than can be indexed in " , sizeof(matrix_index_type) , " bytes" );
#endif // NEW_CODE
		for( int i=0 ; i<matrixBlocks ; i++ ) _constraints[i] = AllocPointer< T >( maxBlockSize ) , _D[i] = AllocPointer< Real >( maxBlockSize );
		for( ; residualWindow.end(!forward)*dir<FullWindow.end(forward)*dir ; residualWindow += dir , solveWindow += dir )
		{
			double t;
			{
				int frontSolveBlock =    solveWindow.begin(!forward);
				int residualBlock   = residualWindow.begin(!forward);
				// Get the leading matrix and compute the constraint norm / initial residual
				// [WARNNG] This is likely wrong. We probably have to pull this into its own for "for( int _c=0 ; _c<ColorModulus ; _c++ )" loop
				//          to ensure that adjacent read-only blocks have not been updated yet.
				if( FullWindow.inBlock( residualBlock ) )
				{
					int b = residualBlock , _b = MOD( b , matrixBlocks );

					t = Time();
#ifdef NEW_THREADS
					_getSliceMatrixAndProlongationConstraints( UIntPack< FEMSigs ... >() , tp , F , _M[_b] , _D[_b] , bsData , depth , _sNodesBegin( depth , BlockFirst( b ) ) , _sNodesEnd( depth , BlockLast( b ) ) , prolongedSolution , _constraints[_b] , ccStencil , pcStencils , interpolationInfo... );
#else // !NEW_THREADS
					_getSliceMatrixAndProlongationConstraints( UIntPack< FEMSigs ... >() , F , _M[_b] , _D[_b] , bsData , depth , _sNodesBegin( depth , BlockFirst( b ) ) , _sNodesEnd( depth , BlockLast( b ) ) , prolongedSolution , _constraints[_b] , ccStencil , pcStencils , interpolationInfo... );
#endif // NEW_THREADS
#ifdef NEW_CODE
#ifdef NEW_THREADS
					size_t begin = _sNodesBegin( depth , BlockFirst( b ) ) , end = _sNodesEnd( depth , BlockLast( b ) );
#ifdef NEW_THREAD_NUM
					tp.parallel_for( begin , end , [&]( const ThreadPool::ThreadNum & , size_t i ){  _constraints[_b][ i-begin ] = constraints[i] - _constraints[_b][ i-begin ]; } );
#else // !NEW_THREAD_NUM
					tp.parallel_for( begin , end , [&]( unsigned int , size_t i ){  _constraints[_b][ i-begin ] = constraints[i] - _constraints[_b][ i-begin ]; } );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
					for( node_index_type i=_sNodesBegin( depth , BlockFirst( b ) ) ; i<_sNodesEnd( depth , BlockLast( b ) ) ; i++ ) _constraints[_b][ i - _sNodesBegin( depth , BlockFirst( b ) ) ] = constraints[i] - _constraints[_b][ i - _sNodesBegin( depth , BlockFirst( b ) ) ];
#endif // NEW_THREADS
					{
						node_index_type begin = _sNodesBegin( depth , BlockFirst( b ) ) , end = _sNodesEnd( depth , BlockLast( b ) );
						for( node_index_type i=begin ; i<end ; i++ ) if( _M[_b].rowSize( i-begin ) ) _D[_b][i-begin] *= sorWeights[i];
					}
#else // !NEW_CODE
#pragma omp parallel for
					for( int i=_sNodesBegin( depth , BlockFirst( b ) ) ; i<_sNodesEnd( depth , BlockLast( b ) ) ; i++ ) _constraints[_b][ i - _sNodesBegin( depth , BlockFirst( b ) ) ] = constraints[i] - _constraints[_b][ i - _sNodesBegin( depth , BlockFirst( b ) ) ];
					{
						int begin = _sNodesBegin( depth , BlockFirst( b ) ) , end = _sNodesEnd( depth , BlockLast( b ) );
						for( int i=begin ; i<end ; i++ ) if( _M[_b].rowSize( i-begin ) ) _D[_b][i-begin] *= sorWeights[i];
					}
#endif // NEW_CODE
					systemTime += Time()-t;
					if( computeNorms )
					{
						ConstPointer( T ) B = _constraints[_b];
						ConstPointer( T ) X = XBlocks( depth , b , solution );
#ifdef NEW_THREADS
						std::vector< double > bNorms( tp.threadNum() , 0 ) , inRNorms( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
						tp.parallel_for( 0 , _M[_b].rows() , [&]( const ThreadPool::ThreadNum &thread , size_t j )
#else // !NEW_THREAD_NUM
						tp.parallel_for( 0 , _M[_b].rows() , [&]( unsigned int thread , size_t j )
#endif //NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : bNorm , inRNorm )
#ifdef NEW_CODE
						for( matrix_index_type j=0 ; j<(matrix_index_type)_M[_b].rows() ; j++ )
#else // !NEW_CODE
						for( int j=0 ; j<_M[_b].rows() ; j++ )
#endif // NEW_CODE
#endif // NEW_THREADS
						{
							T temp = {};
#ifdef NEW_CODE_SPARSE_MATRIX
							ConstPointer( MatrixEntry< Real , matrix_index_type > ) start = _M[_b][j];
							ConstPointer( MatrixEntry< Real , matrix_index_type  > ) end = start + _M[_b].rowSize(j);
							ConstPointer( MatrixEntry< Real , matrix_index_type  > ) e;
#else // !NEW_CODE_SPARSE_MATRIX
							ConstPointer( MatrixEntry< Real > ) start = _M[_b][j];
							ConstPointer( MatrixEntry< Real > ) end = start + (unsigned long long)_M[_b].rowSize(j);
							ConstPointer( MatrixEntry< Real > ) e;
#endif // NEW_CODE_SPARSE_MATRIX
							for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
							bNorms[thread()] += Dot( B[j] , B[j] );
							inRNorms[thread()] += Dot( temp - B[j] , temp - B[j] );
#else // !NEW_THREAD_NUM
							bNorms[thread] += Dot( B[j] , B[j] );
							inRNorms[thread] += Dot( temp - B[j] , temp - B[j] );
#endif //NEW_THREAD_NUM
#else // !NEW_THREADS
							bNorm += Dot( B[j] , B[j] );
							inRNorm += Dot( temp - B[j] , temp - B[j] );
#endif // NEW_THREADS
						}
#ifdef NEW_THREADS
						);
						for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) bNorm += bNorms[t] , inRNorm += inRNorms[t];
#endif // NEW_THREADS
					}
				}
				t = Time();
				// Get the leading multi-color indices
				if( iters && FullWindow.inBlock( frontSolveBlock ) )
				{
					int b = frontSolveBlock , _b = MOD( b , matrixBlocks ) , __b = MOD( b , solveBlocks );
					for( int i=0 ; i<int( mcIndices[__b].size() ) ; i++ ) mcIndices[__b][i].clear();
#ifdef NEW_THREADS
					_setMultiColorIndices( UIntPack< FEMSigs ... >() , tp , _sNodesBegin( depth , BlockFirst( b ) ) , _sNodesEnd( depth , BlockLast( b ) ) , mcIndices[__b] );
#else // !NEW_THREADS
					_setMultiColorIndices( UIntPack< FEMSigs ... >() , _sNodesBegin( depth , BlockFirst( b ) ) , _sNodesEnd( depth , BlockLast( b ) ) , mcIndices[__b] );
#endif // NEW_THREADS
				}
			}

			// Relax the system
			for( int block=solveWindow.begin(!forward) ; solveWindow.inBlock(block) ; block-=dir*ColorModulus ) if( FullWindow.inBlock( block ) )
			{
				int b = block , _b = MOD( b , matrixBlocks ) , __b = MOD( b , solveBlocks );
				ConstPointer( T ) B = _constraints[_b];
				Pointer( T ) X = XBlocks( depth , b , solution );
#ifdef NEW_THREADS
				_M[_b].gsIteration( tp , mcIndices[__b] , ( ConstPointer( Real ) )_D[_b] , B , X , coarseToFine , true );
#else // !NEW_THREADS
				_M[_b].gsIteration( mcIndices[__b] , ( ConstPointer( Real ) )_D[_b] , B , X , coarseToFine , true );
#endif // NEW_THREADS
			}
			solveTime += Time() - t;

			// Compute the final residual
			{
				int residualBlock = residualWindow.begin(forward);
				if( computeNorms && FullWindow.inBlock( residualBlock ) )
				{
					int b = residualBlock , _b = MOD( b , matrixBlocks );
					ConstPointer( T ) B = _constraints[_b];
					ConstPointer( T ) X = XBlocks( depth , b , solution );
#ifdef NEW_THREADS
					std::vector< double > outRNorms( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
					tp.parallel_for( 0 , _M[_b].rows() , [&]( const ThreadPool::ThreadNum &thread , size_t j )
#else // !NEW_THREAD_NUM
					tp.parallel_for( 0 , _M[_b].rows() , [&]( unsigned int thread , size_t j )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : outRNorm )
#ifdef NEW_CODE
					for( matrix_index_type j=0 ; j<(matrix_index_type)_M[_b].rows() ; j++ )
#else // !NEW_CODE
					for( int j=0 ; j<_M[_b].rows() ; j++ )
#endif // NEW_CODE
#endif // NEW_THREADS
					{
						T temp = {};
#ifdef NEW_CODE_SPARSE_MATRIX
						ConstPointer( MatrixEntry< Real , matrix_index_type > ) start = _M[_b][j];
						ConstPointer( MatrixEntry< Real , matrix_index_type > ) end = start + _M[_b].rowSize(j);
						ConstPointer( MatrixEntry< Real , matrix_index_type > ) e;
#else // !NEW_CODE_SPARSE_MATRIX
						ConstPointer( MatrixEntry< Real > ) start = _M[_b][j];
						ConstPointer( MatrixEntry< Real > ) end = start + (unsigned long long)_M[_b].rowSize(j);
						ConstPointer( MatrixEntry< Real > ) e;
#endif // NEW_CODE_SPARSE_MATRIX
						for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
						outRNorms[thread()] += Dot( temp-B[j] , temp-B[j] );
#else // !NEW_THREAD_NUM
						outRNorms[thread] += Dot( temp-B[j] , temp-B[j] );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
						outRNorm += Dot( temp-B[j] , temp-B[j] );
#endif // NEW_THREADS
					}
#ifdef NEW_THREADS
					);
					for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) outRNorm += outRNorms[t];
#endif // NEW_THREADS
				}
			}
		}
		for( int i=0 ; i<matrixBlocks ; i++ ) FreePointer( _D[i] );
		for( int i=0 ; i<matrixBlocks ; i++ ) FreePointer( _constraints[i] );

		if( computeNorms ) stats.bNorm2 = bNorm , stats.inRNorm2 = inRNorm , stats.outRNorm2 = outRNorm;
		DeletePointer( _M );
		DeletePointer( mcIndices );
		FreePointer( _D );
	}
	MemoryUsage();
	return iters;
}
#undef MOD

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , typename TDotT , unsigned int ... PointDs >
#ifdef NEW_THREADS
int FEMTree< Dim , Real >::_solveSystemCG( UIntPack< FEMSigs ... > , ThreadPool &tp , const typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , Pointer( T ) solution , ConstPointer( T ) prolongedSolution , ConstPointer( T ) constraints , TDotT Dot , int iters , bool coarseToFine , _SolverStats& stats , bool computeNorms , double accuracy , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
int FEMTree< Dim , Real >::_solveSystemCG( UIntPack< FEMSigs ... > , const typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , Pointer( T ) solution , ConstPointer( T ) prolongedSolution , ConstPointer( T ) constraints , TDotT Dot , int iters , bool coarseToFine , _SolverStats& stats , bool computeNorms , double accuracy , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_THREADS
{
	int iter = 0;
#ifdef NEW_CODE
	Pointer( T ) X = GetPointer( &solution[0] + _sNodesBegin(depth) , _sNodesSize(depth) );
	ConstPointer( T ) B = GetPointer( &constraints[0] + _sNodesBegin(depth) , _sNodesSize(depth) );
	SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > M;
#else // !NEW_CODE
	Pointer( T ) X = GetPointer( &solution[0] + _sNodesBegin(depth) , _sNodesSize(depth) );
	ConstPointer( T ) B = GetPointer( &constraints[0] + _sNodesBegin(depth) , _sNodesSize(depth) );
	SparseMatrix< Real , int , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > M;
#endif // NEW_CODE

	double& systemTime = stats.systemTime;
	double&  solveTime = stats. solveTime;
	systemTime = solveTime = 0.;
	// Get the system matrix (and adjust the right-hand-side based on the coarser solution if prolonging)
	systemTime = Time();
	Pointer( T ) _constraints = AllocPointer< T >( _sNodesSize( depth ) );
	B = _constraints;
	CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > > ccStencil;
	PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > > pcStencils;
	F.template setStencil< false >( ccStencil );
	F.template setStencils< true >( pcStencils );
#ifdef NEW_THREADS
	_getSliceMatrixAndProlongationConstraints( UIntPack< FEMSigs ... >() , tp , F , M , NullPointer( Real ) , bsData , depth , _sNodesBegin( depth ) , _sNodesEnd( depth ) , prolongedSolution , _constraints , ccStencil , pcStencils , interpolationInfo... );
#else // !NEW_THREADS
	_getSliceMatrixAndProlongationConstraints( UIntPack< FEMSigs ... >() , F , M , NullPointer( Real ) , bsData , depth , _sNodesBegin( depth ) , _sNodesEnd( depth ) , prolongedSolution , _constraints , ccStencil , pcStencils , interpolationInfo... );
#endif // NEW_THREADS
#ifdef NEW_CODE
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd(depth) , [&]( const ThreadPool::ThreadNum & , size_t i ){ _constraints[ i - _sNodesBegin(depth) ] = constraints[i] - _constraints[ i - _sNodesBegin(depth) ]; } );
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd(depth) , [&]( unsigned int , size_t i ){ _constraints[ i - _sNodesBegin(depth) ] = constraints[i] - _constraints[ i - _sNodesBegin(depth) ]; } );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
	for( node_index_type i=_sNodesBegin(depth) ; i<_sNodesEnd(depth) ; i++ ) _constraints[ i - _sNodesBegin(depth) ] = constraints[i] - _constraints[ i - _sNodesBegin(depth) ];
#endif // NEW_THREADS
#else // !NEW_CODE
#pragma omp parallel for
	for( int i=_sNodesBegin(depth) ; i<_sNodesEnd(depth) ; i++ ) _constraints[ i - _sNodesBegin(depth) ] = constraints[i] - _constraints[ i - _sNodesBegin(depth) ];
#endif // NEW_CODE
	systemTime = Time()-systemTime;
	solveTime = Time();
	// Solve the linear system
	accuracy = Real( accuracy / 100000 ) * M.rows();
	int dims[] = { ( _BSplineEnd< FEMSigs >( depth ) - _BSplineBegin< FEMSigs >( depth ) ) ... };
#ifdef NEW_CODE
	size_t nonZeroRows = 0;
	for( matrix_index_type i=0 ; i<(matrix_index_type)M.rows() ; i++ ) if( M.rowSize(i) ) nonZeroRows++;
	size_t totalDim = 1;
#else // !NEW_CODE
	int nonZeroRows = 0;
	for( int i=0 ; i<M.rows() ; i++ ) if( M.rowSize(i) ) nonZeroRows++;
	int totalDim = 1;
#endif // NEW_CODE
	for( int d=0 ; d<Dim ; d++ ) totalDim *= dims[d];
	BoundaryType bTypes[] = { FEMSignature< FEMSigs >::BType ... };
	bool hasPartitionOfUnity = true;
	for( int d=0 ; d<Dim ; d++ ) hasPartitionOfUnity &= HasPartitionOfUnity( bTypes[d] );
	bool addDCTerm = ( nonZeroRows==totalDim && !ConstrainsDCTerm( interpolationInfo... ) && hasPartitionOfUnity && F.vanishesOnConstants() );
	double bNorm = 0 , inRNorm = 0 , outRNorm = 0;
	if( computeNorms )
	{
#ifdef NEW_THREADS
		std::vector< double > bNorms( tp.threadNum() , 0 ) , inRNorms( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
		tp.parallel_for( 0 , M.rows() , [&]( const ThreadPool::ThreadNum &thread , size_t j )
#else // !NEW_THREAD_NUM
		tp.parallel_for( 0 , M.rows() , [&]( unsigned int thread , size_t j )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : bNorm , inRNorm )
#ifdef NEW_CODE
		for( matrix_index_type j=0 ; j<(matrix_index_type)M.rows() ; j++ )
#else // !NEW_CODE
		for( int j=0 ; j<M.rows() ; j++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			T temp = {};
#ifdef NEW_CODE_SPARSE_MATRIX
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) start = M[j];
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) end = start + M.rowSize(j);
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) e;
#else // !NEW_CODE_SPARSE_MATRIX
			ConstPointer( MatrixEntry< Real > ) start = M[j];
			ConstPointer( MatrixEntry< Real > ) end = start + (unsigned long long)M.rowSize(j);
			ConstPointer( MatrixEntry< Real > ) e;
#endif // NEW_CODE_SPARSE_MATRIX
			for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			bNorms[thread()] += Dot( B[j] , B[j] );
			inRNorms[thread()] += Dot( temp-B[j] , temp-B[j] );
#else // !NEW_THREAD_NUM
			bNorms[thread] += Dot( B[j] , B[j] );
			inRNorms[thread] += Dot( temp-B[j] , temp-B[j] );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			bNorm += Dot( B[j] , B[j] );
			inRNorm += Dot( temp-B[j] , temp-B[j] );
#endif // NEW_THREADS
		}
#ifdef NEW_THREADS
		);
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) bNorm += bNorms[t] , inRNorm += inRNorms[t];
#endif // NEW_THREADS
	}

#ifdef NEW_CODE
	iters = (int)std::min< size_t >( nonZeroRows , iters );
#else // !NEW_CODE
	iters = std::min< int >( nonZeroRows , iters );
#endif // NEW_CODE
	struct SPDFunctor
	{
	protected:
#ifdef NEW_CODE
		const SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size >& _M;
#else // !NEW_CODE
		const SparseMatrix< Real , int , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size >& _M;
#endif // NEW_CODE
#ifdef NEW_THREADS
		ThreadPool &_tp;
#endif // NEW_THREADS
		bool _addDCTerm;
	public:
#ifdef NEW_CODE
#ifdef NEW_THREADS
		SPDFunctor( ThreadPool &tp , const SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size >& M , bool addDCTerm ) : _M(M) , _addDCTerm(addDCTerm) , _tp(tp){ }
#else // !NEW_THREADS
		SPDFunctor( const SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size >& M , bool addDCTerm ) : _M(M) , _addDCTerm(addDCTerm){ }
#endif // NEW_THREADS
#else // !NEW_CODE
		SPDFunctor( const SparseMatrix< Real , int , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size >& M , bool addDCTerm ) : _M(M) , _addDCTerm(addDCTerm){ }
#endif // NEW_CODE
		void operator()( ConstPointer( T ) in , Pointer( T ) out ) const
		{
#ifdef NEW_THREADS
			_M.multiply( _tp , in , out );
#else // !NEW_THREADS
			_M.multiply( in , out );
#endif // NEW_THREADS
			if( _addDCTerm )
			{
				T average = {};
#ifdef NEW_CODE
				for( matrix_index_type i=0 ; i<(matrix_index_type)_M.rows() ; i++ ) average += in[i];
				average /= _M.rows();
				for( matrix_index_type i=0 ; i<(matrix_index_type)_M.rows() ; i++ ) out[i] += average;
#else // !NEW_CODE
				for( int i=0 ; i<_M.rows() ; i++ ) average += in[i];
				average /= _M.rows();
				for( int i=0 ; i<_M.rows() ; i++ ) out[i] += average;
#endif // NEW_CODE
			}
		}
	};
#ifdef NEW_CODE
#ifdef NEW_THREADS
	if( iters ) iter = (int)SolveCG< SPDFunctor , T , Real >( tp , SPDFunctor( tp , M , addDCTerm ) , M.rows() , ( ConstPointer( T ) )B , iters , X , Real( accuracy ) , Dot );
#else // !NEW_THREADS
	if( iters ) iter = (int)SolveCG< SPDFunctor , T , Real >( SPDFunctor( M , addDCTerm ) , M.rows() , ( ConstPointer( T ) )B , iters , X , Real( accuracy ) , Dot );
#endif // NEW_THREADS
#else // !NEW_CODE
	if( iters ) iter = SolveCG< SPDFunctor , T , Real >( SPDFunctor( M , addDCTerm ) , (int)M.rows() , ( ConstPointer( T ) )B , iters , X , Real( accuracy ) , Dot );
#endif // NEW_CODE

	solveTime = Time()-solveTime;
	if( computeNorms )
	{
#ifdef NEW_THREADS
		std::vector< double > outRNorms( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
		tp.parallel_for( 0 , M.rows() , [&]( const ThreadPool::ThreadNum &thread , size_t j )
#else // !NEW_THREAD_NUM
		tp.parallel_for( 0 , M.rows() , [&]( unsigned int thread , size_t j )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : outRNorm )
#ifdef NEW_CODE
		for( matrix_index_type j=0 ; j<(matrix_index_type)M.rows() ; j++ )
#else // !NEW_CODE
		for( int j=0 ; j<M.rows() ; j++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			T temp = {};
#ifdef NEW_CODE_SPARSE_MATRIX
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) start = M[j];
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) end = start + M.rowSize(j);
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) e;
#else // !NEW_CODE_SPARSE_MATRIX
			ConstPointer( MatrixEntry< Real > ) start = M[j];
			ConstPointer( MatrixEntry< Real > ) end = start + (unsigned long long)M.rowSize(j);
			ConstPointer( MatrixEntry< Real > ) e;
#endif // NEW_CODE_SPARSE_MATRIX
			for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			outRNorms[ thread() ] += Dot( temp-B[j] , temp-B[j] );
#else // !NEW_THREAD_NUM
			outRNorms[ thread ] += Dot( temp-B[j] , temp-B[j] );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			outRNorm += Dot( temp-B[j] , temp-B[j] );
#endif // NEW_THREADS
		}
#ifdef NEW_THREADS
		);
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) outRNorm += outRNorms[t];
#endif // NEW_THREADS
		stats.bNorm2 = bNorm , stats.inRNorm2 = inRNorm , stats.outRNorm2 = outRNorm;
	}
	FreePointer( _constraints );

	MemoryUsage();
	return iter;
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , typename TDotT , unsigned int ... PointDs >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_solveRegularMG( UIntPack< FEMSigs ... > , ThreadPool &tp , typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , Pointer( T ) solution , ConstPointer( T ) constraints , TDotT Dot , int vCycles , int iters , _SolverStats& stats , bool computeNorms , double cgAccuracy , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_solveRegularMG( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , Pointer( T ) solution , ConstPointer( T ) constraints , TDotT Dot , int vCycles , int iters , _SolverStats& stats , bool computeNorms , double cgAccuracy , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_THREADS
{
	double& systemTime = stats.systemTime;
	double&  solveTime = stats. solveTime;

#ifdef NEW_CODE
	std::vector< SparseMatrix< Real , matrix_index_type > > P( depth ) , R( depth ) , M( depth+1 );
#else // !NEW_CODE
	std::vector< SparseMatrix< Real , int > > P( depth ) , R( depth ) , M( depth+1 );
#endif // NEW_CODE
	std::vector< Pointer( Real ) > D( depth+1 );
	std::vector< Pointer( T ) > B( depth+1 ) , X( depth+1 ) , MX( depth+1 );
#ifdef NEW_CODE
	std::vector< std::vector< std::vector< size_t > > > multiColorIndices( depth+1 );
#else // !NEW_CODE
	std::vector< std::vector< std::vector< int > > > multiColorIndices( depth+1 );
#endif // NEW_CODE

	systemTime = Time();
#ifdef NEW_THREADS
	M.back() = systemMatrix< Real >( UIntPack< FEMSigs ... >() , tp , F , depth , interpolationInfo ... );
#else // !NEW_THREADS
	M.back() = systemMatrix< Real >( UIntPack< FEMSigs ... >() , F , depth , interpolationInfo ... );
#endif // NEW_THREADS
	for( int d=depth ; d>0 ; d-- )
	{
#ifdef NEW_THREADS
		R[d-1] = downSampleMatrix( UIntPack< FEMSigs ... >() , tp , d );
#else // !NEW_THREADS
		R[d-1] = downSampleMatrix( UIntPack< FEMSigs ... >() , d );
#endif // NEW_THREADS
		P[d-1] = R[d-1].transpose();
		M[d-1] = R[d-1] * M[d] * P[d-1];
	}
	for( int d=0 ; d<=depth ; d++ )
	{
		size_t dim = M[d].rows();
		D[d]  = AllocPointer< Real >( dim );
		MX[d] = AllocPointer< T >( dim );
#ifdef NEW_THREADS
		M[d].setDiagonalR( tp , D[d] );
		setMultiColorIndices( UIntPack< FEMSigs ... >() , tp , d , multiColorIndices[d] );
#else // !NEW_THREADS
		M[d].setDiagonalR( D[d] );
		setMultiColorIndices( UIntPack< FEMSigs ... >() , d , multiColorIndices[d] );
#endif // NEW_THREADS
		if( d<depth )
		{
			X[d]  = AllocPointer< T >( dim );
			B[d]  = AllocPointer< T >( dim );
		}
	}
	X.back() = solution + nodesBegin( depth );
	ConstPointer( T ) _B = constraints + nodesBegin( depth );
	systemTime = Time() - systemTime;

	solveTime = Time();

	double bNorm = 0 , inRNorm = 0 , outRNorm = 0;
	if( computeNorms )
	{
#ifdef NEW_CODE
		const SparseMatrix< Real , matrix_index_type >& _M = M.back();
#else // !NEW_CODE
		const SparseMatrix< Real , int >& _M = M.back();
#endif // NEW_CODE
		ConstPointer( T ) _X = X.back();
#ifdef NEW_THREADS
		std::vector< double > bNorms( tp.threadNum() , 0 ) , inRNorms( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( const ThreadPool::ThreadNum &thread , size_t j )
#else // !NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( unsigned int thread , size_t j )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : bNorm , inRNorm )
#ifdef NEW_CODE
		for( matrix_index_type j=0 ; j<_M.rows() ; j++ )
#else // !NEW_CODE
		for( int j=0 ; j<_M.rows() ; j++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			T temp = {};
#ifdef NEW_CODE_SPARSE_MATRIX
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) start = _M[j];
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) end = start + _M.rowSize(j);
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) e;
#else // !NEW_CODE_SPARSE_MATRIX
			ConstPointer( MatrixEntry< Real > ) start = _M[j];
			ConstPointer( MatrixEntry< Real > ) end = start + (unsigned long long)_M.rowSize(j);
			ConstPointer( MatrixEntry< Real > ) e;
#endif // NEW_CODE_SPARSE_MATRIX
			for( e=start ; e!=end ; e++ ) temp += _X[ e->N ] * e->Value;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			bNorms[thread()] += Dot( _B[j] , _B[j] );
			inRNorms[thread()] += Dot( temp-_B[j] , temp-_B[j] );
#else // !NEW_THREAD_NUM
			bNorms[thread] += Dot( _B[j] , _B[j] );
			inRNorms[thread] += Dot( temp-_B[j] , temp-_B[j] );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			bNorm += Dot( _B[j] , _B[j] );
			inRNorm += Dot( temp-_B[j] , temp-_B[j] );
#endif // NEW_THREADS
		}
#ifdef NEW_THREADS
		);
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) bNorm += bNorms[t] , inRNorm += inRNorms[t];
#endif // NEW_THREADS
	}

	for( int v=0 ; v<vCycles ; v++ )
	{
		// Restriction
		for( int d=depth ; d>0 ; d-- )
		{
			ConstPointer( T ) __B = d==depth ? _B : B[d];
#ifdef NEW_THREADS
			for( int i=0 ; i<iters ; i++ ) M[d].gsIteration( tp , multiColorIndices[d] , D[d] , __B , X[d] , true , true );
			M[d].multiply( tp , X[d] , MX[d] );
#else // !NEW_THREADS
			for( int i=0 ; i<iters ; i++ ) M[d].gsIteration( multiColorIndices[d] , D[d] , __B , X[d] , true , true );
			M[d].multiply( X[d] , MX[d] );
#endif // NEW_THREADS
#ifdef NEW_CODE
			for( matrix_index_type i=0 ; i<(matrix_index_type)M[d].rows() ; i++ ) MX[d][i] = __B[i] - MX[d][i];
#else // !NEW_CODE
			for( int i=0 ; i<M[d].rows() ; i++ ) MX[d][i] = __B[i] - MX[d][i];
#endif // NEW_CODE
#ifdef NEW_THREADS
			R[d-1].multiply( tp , MX[d] , B[d-1] );
#else // !NEW_THREADS
			R[d-1].multiply( MX[d] , B[d-1] );
#endif // NEW_THREADS
			memset( X[d-1] , 0 , sizeof( T )*M[d-1].rows() );
		}

		// Base
		{
			int d = 0;
			ConstPointer( T ) __B = d==depth ? _B : B[d];
			struct SPDFunctor
			{
			protected:
#ifdef NEW_CODE
				const SparseMatrix< Real , matrix_index_type >& _M;
#else // !NEW_CODE
				const SparseMatrix< Real , int >& _M;
#endif // NEW_CODE
#ifdef NEW_THREADS
				ThreadPool &_tp;
#endif // NEW_THREADS
				bool _addDCTerm;
			public:
#ifdef NEW_CODE
#ifdef NEW_THREADS
				SPDFunctor( ThreadPool &tp , const SparseMatrix< Real , matrix_index_type  >& M , bool addDCTerm ) : _M(M) , _addDCTerm(addDCTerm) , _tp(tp){ }
#else // !NEW_THREADS
				SPDFunctor( const SparseMatrix< Real , matrix_index_type  >& M , bool addDCTerm ) : _M(M) , _addDCTerm(addDCTerm){ }
#endif // NEW_THREADS
#else // !NEW_CODE
				SPDFunctor( const SparseMatrix< Real , int  >& M , bool addDCTerm ) : _M(M) , _addDCTerm(addDCTerm){ }
#endif // NEW_CODE
				void operator()( ConstPointer( T ) in , Pointer( T ) out ) const
				{
#ifdef NEW_THREADS
					_M.multiply( _tp , in , out );
#else // !NEW_THREADS
					_M.multiply( in , out );
#endif // NEW_THREADS
					if( _addDCTerm )
					{
						T average = {};
#ifdef NEW_CODE
						for( matrix_index_type i=0 ; i<(matrix_index_type)_M.rows() ; i++ ) average += in[i];
						average /= _M.rows();
						for( matrix_index_type i=0 ; i<(matrix_index_type)_M.rows() ; i++ ) out[i] += average;
#else // !NEW_CODE
						for( int i=0 ; i<_M.rows() ; i++ ) average += in[i];
						average /= _M.rows();
						for( int i=0 ; i<_M.rows() ; i++ ) out[i] += average;
#endif // NEW_CODE
					}
				}
			};
#ifdef NEW_CODE
			size_t nonZeroRows = 0;
			for( matrix_index_type i=0 ; i<(matrix_index_type)M[d].rows() ; i++ ) if( M[d].rowSize(i) ) nonZeroRows++;
			size_t totalDim = 1;
#else // !NEW_CODE
			int nonZeroRows = 0;
			for( int i=0 ; i<M[d].rows() ; i++ ) if( M[d].rowSize(i) ) nonZeroRows++;
			int totalDim = 1;
#endif // NEW_CODE
			int dims[] = { ( _BSplineEnd< FEMSigs >( depth ) - _BSplineBegin< FEMSigs >( depth ) ) ... };
			for( int dd=0 ; dd<Dim ; dd++ ) totalDim *= dims[dd];
			BoundaryType bTypes[] = { FEMSignature< FEMSigs >::BType ... };
			bool hasPartitionOfUnity = true;
			for( int dd=0 ; dd<Dim ; dd++ ) hasPartitionOfUnity &= HasPartitionOfUnity( bTypes[dd] );
			bool addDCTerm = ( nonZeroRows==totalDim && !ConstrainsDCTerm( interpolationInfo... ) && hasPartitionOfUnity && F.vanishesOnConstants() );

#ifdef NEW_CODE
#ifdef NEW_THREADS
			SolveCG< SPDFunctor , T , Real >( tp , SPDFunctor( tp,  M[d] , addDCTerm ) , M[d].rows() , ( ConstPointer( T ) )__B , nonZeroRows , X[d] , Real( cgAccuracy ) , Dot );
#else // !NEW_THREADS
			SolveCG< SPDFunctor , T , Real >( SPDFunctor( M[d] , addDCTerm ) , M[d].rows() , ( ConstPointer( T ) )__B , nonZeroRows , X[d] , Real( cgAccuracy ) , Dot );
#endif // NEW_THREADS
#else // !NEW_CODE
			SolveCG< SPDFunctor , T , Real >( SPDFunctor( M[d] , addDCTerm ) , (int)M[d].rows() , ( ConstPointer( T ) )__B , nonZeroRows , X[d] , Real( cgAccuracy ) , Dot );
#endif // NEW_CODE
		}

		// Prolongation
		for( int d=1 ; d<=depth ; d++ )
		{
			ConstPointer( T ) __B = d==depth ? _B : B[d];
#ifdef NEW_THREADS
			P[d-1].multiply( tp , X[d-1] , X[d] , MULTIPLY_ADD );
			for( int i=0 ; i<iters ; i++ ) M[d].gsIteration( tp , multiColorIndices[d] , D[d] , __B , X[d] , false , true );
#else // !NEW_THREADS
			P[d-1].multiply( X[d-1] , X[d] , MULTIPLY_ADD );
			for( int i=0 ; i<iters ; i++ ) M[d].gsIteration( multiColorIndices[d] , D[d] , __B , X[d] , false , true );
#endif // NEW_THREADS
		}
	}
	if( computeNorms )
	{
#ifdef NEW_CODE
		const SparseMatrix< Real , matrix_index_type >& _M = M.back();
#else // !NEW_CODE
		const SparseMatrix< Real , int >& _M = M.back();
#endif // NEW_CODE
		ConstPointer( T ) _X = X.back();
#ifdef NEW_THREADS
		std::vector< double > outRNorms( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( const ThreadPool::ThreadNum &thread , size_t j )
#else // !NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( unsigned int thread , size_t j )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : outRNorm )
#ifdef NEW_CODE
		for( matrix_index_type j=0 ; j<(matrix_index_type)_M.rows() ; j++ )
#else // !NEW_CODE
		for( int j=0 ; j<_M.rows() ; j++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			T temp = {};
#ifdef NEW_CODE_SPARSE_MATRIX
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) start = _M[j];
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) end = start + _M.rowSize(j);
			ConstPointer( MatrixEntry< Real , matrix_index_type > ) e;
#else // !NEW_CODE_SPARSE_MATRIX
			ConstPointer( MatrixEntry< Real > ) start = _M[j];
			ConstPointer( MatrixEntry< Real > ) end = start + (unsigned long long)_M.rowSize(j);
			ConstPointer( MatrixEntry< Real > ) e;
#endif // NEW_CODE_SPARSE_MATRIX
			for( e=start ; e!=end ; e++ ) temp += _X[ e->N ] * e->Value;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			outRNorms[thread()] += Dot( temp-_B[j] , temp-_B[j] );
#else // !NEW_THREAD_NUM
			outRNorms[thread] += Dot( temp-_B[j] , temp-_B[j] );
#endif //NEW_THREAD_NUM
#else // !NEW_THREADS
			outRNorm += Dot( temp-_B[j] , temp-_B[j] );
#endif // NEW_THREADS
		}
#ifdef NEW_THREADS
		);
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) outRNorm += outRNorms[t];
#endif // NEW_THREADS
		stats.bNorm2 = bNorm , stats.inRNorm2 = inRNorm , stats.outRNorm2 = outRNorm;
	}
	solveTime = Time() - solveTime;

	MemoryUsage();

	for( int d=0 ; d<=depth ; d++ )
	{
		FreePointer( D[d] );
		FreePointer( MX[d] );
		if( d<depth )
		{
			FreePointer( X[d] );
			FreePointer( B[d] );
		}
	}
}

#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
int FEMTree< Dim , Real >::_getMatrixRowSize( UIntPack< FEMSigs ... > , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors ) const
#else // !__GNUC__ || __GNUC__ >=5
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
int FEMTree< Dim , Real >::_getMatrixRowSize( const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors ) const
#endif // __GNUC__ || __GNUC__ < 4
{
	typedef UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > OverlapSizes;

	int count = 0;
	const FEMTreeNode* const * _nodes = neighbors.neighbors.data;
	for( int i=0 ; i<WindowSize< OverlapSizes >::Size ; i++ ) if( _isValidFEM1Node( _nodes[i] ) ) count++;
	return count;
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
int FEMTree< Dim , Real >::_getProlongedMatrixRowSize( const FEMTreeNode* node , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors ) const
{
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > FEMDegrees;
	typedef UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > OverlapSizes;
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] This change needs to be validated" )
#endif // SHOW_WARNINGS
	int count = 0;
	static const WindowLoopData< OverlapSizes > loopData( []( int c , int* start , int*end ){ _SetParentOverlapBounds( FEMDegrees() , FEMDegrees() , c , start , end );} );
	if( node->parent )
	{
		int c =  (int)( node - node->parent->children );
		const unsigned int size = loopData.size[c];
		const unsigned int* indices = loopData.indices[c];
		ConstPointer( FEMTreeNode * const ) nodes = pNeighbors.neighbors().data;
		for( unsigned int i=0 ; i<size ; i++ ) if( _isValidFEM1Node( nodes[ indices[i] ] ) ) count++;
	}
	return count;
}


// Given a node:
// -- For each of its neighbors:
// ---- Compute the weighted sum of the product of the evaluations of the associated basis functions over the points
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , unsigned int PointD >
void FEMTree< Dim , Real >::_addPointValues( UIntPack< FEMSigs ... > , StaticWindow< Real , UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pointValues , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , const InterpolationInfo< T , PointD >* interpolationInfo ) const
{
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > FEMDegrees;
	typedef UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > OverlapSizes;
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	typedef UIntPack< ( -BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportStart ) ... > LeftSupportRadii;
	typedef UIntPack<    BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportEnd     ... > RightSupportRadii;
	typedef UIntPack<    BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportEnd     ... > LeftPointSupportRadii;
	typedef UIntPack< ( -BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportStart ) ... > RightPointSupportRadii;
	typedef UIntPack<    BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize    ... > SupportSizes;

	if( !( FEMDegrees() >= IsotropicUIntPack< Dim , PointD >() ) ) ERROR_OUT( "Insufficient derivatives" );
	if( !interpolationInfo ) return;
	const InterpolationInfo< T , PointD >& iInfo = *interpolationInfo;


	const FEMTreeNode* node = neighbors.neighbors.data[ WindowIndex< OverlapSizes , OverlapRadii >::Index ];
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );

	PointEvaluatorState< UIntPack< FEMSigs ... > , IsotropicUIntPack< Dim , PointD > > peState;

	int idx[Dim];	// The coordinates of the node containing the point _relative_ to the center node
	int _idx[Dim==1 ? 1 : Dim-1];
	CumulativeDerivativeValues< double , Dim , PointD > dualValues;

	auto outerFunction = [&]( const FEMTreeNode* _node  )
	{
		if( _isValidSpaceNode( _node ) )
		{
			LocalOffset pOff;	// The coordinates of the node containing the point
			for( int d=0 ; d<Dim ; d++ ) pOff[d] = off[d] + idx[d];
			size_t begin , end;
			iInfo.range( _node , begin , end );
			for( size_t pIndex=begin ; pIndex<end ; pIndex++ )
			{
				const DualPointInfo< Dim , Real , T , PointD >& pData = iInfo[ pIndex ];
				CumulativeDerivativeValues< double , Dim , PointD > values;
				{
					Real weight = pData.weight;
					Point< Real , Dim > p = pData.position;
					// Compute the partial evaluation of all B-splines (and derivatives) that are supported on the point
					bsData.initEvaluationState( p , d , pOff , peState );

					// The value (and derivatives) of the function of the center node at this point
					values = peState.template dValues< Real , CumulativeDerivatives< Dim , PointD > >( off );
				}
				dualValues = iInfo( pIndex , values ) * pData.weight;
				int start[Dim==1 ? 1 : Dim-1] , end[Dim==1 ? 1 : Dim-1];
				// Compute the bounds of nodes which can be supported on the point
				for( int d=0 ; d<Dim-1 ; d++ ) start[d] = idx[d] + (int)OverlapRadii::Values[d] - (int)LeftPointSupportRadii::Values[d] , end[d] = idx[d] + (int)OverlapRadii::Values[d] + (int)RightPointSupportRadii::Values[d] + 1;
				WindowLoop< Dim , Dim-1 >::Run
				(
					start , end , 
					[&]( int d , int i ){ _idx[d] = i - (int)OverlapRadii::Values[d] + off[d]; } ,
					[&]( const WindowSlice< Real , UIntPack< OverlapSizes::template Get< Dim-1 >() > > pointValues , ConstWindowSlice< const FEMTreeNode* , UIntPack< OverlapSizes::template Get< Dim-1 >() > > neighbors )
					{
						Point< double , PointD+1 > partialDot = peState.template partialDotDValues< Real , CumulativeDerivatives< Dim , PointD > >( dualValues , _idx );
						Pointer( Real ) _pointValues = pointValues.data + idx[Dim-1] + OverlapRadii::Values[Dim-1];

						int _i = idx[Dim-1] + (int)OverlapRadii::Values[Dim-1] - (int)LeftPointSupportRadii::Values[Dim-1];
						const double (*splineValues)[PointD+1] = peState.template values< Dim-1 >();
						for( unsigned int i=0 ; i<SupportSizes::Values[Dim-1] ; i++ ) if( _isValidFEM1Node( neighbors[ _i + i ] ) )
							for( int d=0 ; d<=PointD ; d++ ) _pointValues[(int)i-(int)LeftPointSupportRadii::Values[Dim-1]] += (Real)( splineValues[i][d] * partialDot[d] );
					} ,
					pointValues() , neighbors.neighbors()
				);
			}
		}
	};
	// Loop over all nodes which are supported on the center
	WindowLoop< Dim >::Run
	(
		OverlapRadii() - LeftSupportRadii() , OverlapRadii() + RightSupportRadii() + IsotropicUIntPack< Dim , 1 >() ,
		[&]( int d , int i ){ idx[d] = i - (int)OverlapRadii::Values[d]; } ,
		outerFunction ,
		neighbors.neighbors()
	);
}

template< unsigned int Dim , class Real >
template< typename T , unsigned int ... PointDs , unsigned int ... FEMSigs >
#ifdef NEW_CODE
T FEMTree< Dim , Real >::_setMatrixRowAndGetConstraintFromProlongation( UIntPack< FEMSigs ... > , const BaseSystem< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , size_t idx , SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > &M , node_index_type offset , const PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& pcStencils , const CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& ccStencil , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , ConstPointer( T ) prolongedSolution , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_CODE
T FEMTree< Dim , Real >::_setMatrixRowAndGetConstraintFromProlongation( UIntPack< FEMSigs ... > , const BaseSystem< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , size_t idx , SparseMatrix< Real , int , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size > &M , int offset , const PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& pcStencils , const CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& ccStencil , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , ConstPointer( T ) prolongedSolution , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_CODE
{
	T constraint ={};
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	typedef UIntPack<    BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize    ... > OverlapSizes;

	int count = 0;
	const FEMTreeNode* node = neighbors.neighbors.data[ WindowIndex< OverlapSizes , OverlapRadii >::Index ];
#ifdef NEW_CODE_SPARSE_MATRIX
	Pointer( MatrixEntry< Real , matrix_index_type > ) row = M[idx];
#else // !NEW_CODE_SPARSE_MATRIX
	Pointer( MatrixEntry< Real > ) row = M[idx];
#endif // NEW_CODE_SPARSE_MATRIX

	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	if( d>0 && prolongedSolution )
	{
		int cIdx = (int)( node - node->parent->children );
		constraint = _getConstraintFromProlongedSolution( UIntPack< FEMSigs ... >() , F , neighbors , pNeighbors , node , prolongedSolution , pcStencils.data[cIdx] , bsData , interpolationInfo... );
	}

	bool isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< FEMSignature< FEMSigs >::Degree ... >() , UIntPack< FEMSignature< FEMSigs >::Degree ... >() , d , off );

	StaticWindow< Real , UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > > pointValues;
	memset( pointValues.data , 0 , sizeof(Real)*WindowSize< OverlapSizes >::Size );
	_addPointValues( UIntPack< FEMSigs ... >() , pointValues , neighbors , bsData , interpolationInfo ... );
#ifdef NEW_CODE
	node_index_type nodeIndex = node->nodeData.nodeIndex;
#else // !NEW_CODE
	int nodeIndex = node->nodeData.nodeIndex;
#endif // NEW_CODE
	if( isInterior ) // General case, so try to make fast
	{
		const FEMTreeNode* const * _nodes = neighbors.neighbors.data;
		ConstPointer( double ) _stencil = ccStencil.data;
		Real* _values = pointValues.data;
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
		row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( nodeIndex-offset ) , (Real)( _values[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] + _stencil[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] ) );
#else // !NEW_CODE_SPARSE_MATRIX
		row[count++] = MatrixEntry< Real >( (matrix_index_type)( nodeIndex-offset ) , (Real)( _values[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] + _stencil[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] ) );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
		row[count++] = MatrixEntry< Real >( nodeIndex-offset , (Real)( _values[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] + _stencil[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] ) );
#endif // NEW_CODE
		for( int i=0 ; i<WindowSize< OverlapSizes >::Size ; i++ ) if( _isValidFEM1Node( _nodes[i] ) )
		{
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
			if( i!=WindowIndex< OverlapSizes , OverlapRadii >::Index ) row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( _nodes[i]->nodeData.nodeIndex-offset ) , (Real)( _values[i] + _stencil[i] ) );
#else // !NEW_CODE_SPARSE_MATRIX
			if( i!=WindowIndex< OverlapSizes , OverlapRadii >::Index ) row[count++] = MatrixEntry< Real >( (matrix_index_type)( _nodes[i]->nodeData.nodeIndex-offset ) , (Real)( _values[i] + _stencil[i] ) );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
			if( i!=WindowIndex< OverlapSizes , OverlapRadii >::Index ) row[count++] = MatrixEntry< Real >( _nodes[i]->nodeData.nodeIndex-offset , (Real)( _values[i] + _stencil[i] ) );
#endif // NEW_CODE
		}
	}
	else
	{
		LocalDepth d ; LocalOffset off;
		_localDepthAndOffset( node , d , off );
		Real temp = (Real)F.ccIntegrate( off , off ) + pointValues.data[ WindowIndex< OverlapSizes , OverlapRadii >::Index ];

#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
		row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( nodeIndex-offset ) , temp );
#else // !NEW_CODE_SPARSE_MATRIX
		row[count++] = MatrixEntry< Real >( (matrix_index_type)( nodeIndex-offset ) , temp );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
		row[count++] = MatrixEntry< Real >( nodeIndex-offset , temp );
#endif // NEW_CODE
		LocalOffset _off;
		WindowLoop< Dim >::Run
		(
			ZeroUIntPack< Dim >() , OverlapSizes() ,
			[&]( int d , int i ){ _off[d] = off[d] - (int)OverlapRadii::Values[d] + i; } ,
			[&]( const FEMTreeNode* _node , Real pointValue )
		{
			if( node!=_node && FEMIntegrator::IsValidFEMNode( UIntPack< FEMSigs ... >() , d , _off ) )
			{
				Real temp = (Real)F.ccIntegrate( _off , off ) + pointValue;
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
				if( _isValidFEM1Node( _node ) ) row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( _node->nodeData.nodeIndex-offset ) , temp );
#else // !NEW_CODE_SPARSE_MATRIX
				if( _isValidFEM1Node( _node ) ) row[count++] = MatrixEntry< Real >( (matrix_index_type)( _node->nodeData.nodeIndex-offset ) , temp );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
				if( _isValidFEM1Node( _node ) ) row[count++] = MatrixEntry< Real >( _node->nodeData.nodeIndex-offset , temp );
#endif // NEW_CODE
			}
		} ,
			neighbors.neighbors() , pointValues()
			);
	}
	M.setRowSize( idx , count );
	return constraint;
}

template< unsigned int Dim , class Real >
template< typename T , unsigned int ... PointDs , unsigned int ... FEMSigs >
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
T FEMTree< Dim , Real >::_setMatrixRowAndGetConstraintFromProlongation( UIntPack< FEMSigs ... > , const BaseSystem< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , Pointer( MatrixEntry< Real , matrix_index_type > ) row , node_index_type offset , const PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& pcStencils , const CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& ccStencil , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , ConstPointer( T ) prolongedSolution , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_CODE_SPARSE_MATRIX
T FEMTree< Dim , Real >::_setMatrixRowAndGetConstraintFromProlongation( UIntPack< FEMSigs ... > , const BaseSystem< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , Pointer( MatrixEntry< Real > ) row , node_index_type offset , const PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& pcStencils , const CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& ccStencil , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , ConstPointer( T ) prolongedSolution , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
T FEMTree< Dim , Real >::_setMatrixRowAndGetConstraintFromProlongation( UIntPack< FEMSigs ... > , const BaseSystem< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , Pointer( MatrixEntry< Real > ) row , int offset , const PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& pcStencils , const CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& ccStencil , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , ConstPointer( T ) prolongedSolution , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_CODE
{
	T constraint ={};
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	typedef UIntPack<    BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize    ... > OverlapSizes;

	int count = 0;
	const FEMTreeNode* node = neighbors.neighbors.data[ WindowIndex< OverlapSizes , OverlapRadii >::Index ];
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	if( d>0 && prolongedSolution )
	{
		int cIdx = (int)( node - node->parent->children );
		constraint = _getConstraintFromProlongedSolution( UIntPack< FEMSigs ... >() , F , neighbors , pNeighbors , node , prolongedSolution , pcStencils.data[cIdx] , bsData , interpolationInfo... );
	}

	bool isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< FEMSignature< FEMSigs >::Degree ... >() , UIntPack< FEMSignature< FEMSigs >::Degree ... >() , d , off );

	StaticWindow< Real , UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > > pointValues;
	memset( pointValues.data , 0 , sizeof(Real)*WindowSize< OverlapSizes >::Size );
	_addPointValues( UIntPack< FEMSigs ... >() , pointValues , neighbors , bsData , interpolationInfo ... );
#ifdef NEW_CODE
	node_index_type nodeIndex = node->nodeData.nodeIndex;
#else // !NEW_CODE
	int nodeIndex = node->nodeData.nodeIndex;
#endif // NEW_CODE
	if( isInterior ) // General case, so try to make fast
	{
		const FEMTreeNode* const * _nodes = neighbors.neighbors.data;
		ConstPointer( double ) _stencil = ccStencil.data;
		Real* _values = pointValues.data;
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
		row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( nodeIndex-offset ) , (Real)( _values[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] + _stencil[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] ) );
#else // !NEW_CODE_SPARSE_MATRIX
		row[count++] = MatrixEntry< Real >( (matrix_index_type)( nodeIndex-offset ) , (Real)( _values[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] + _stencil[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] ) );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
		row[count++] = MatrixEntry< Real >( nodeIndex-offset , (Real)( _values[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] + _stencil[ WindowIndex< OverlapSizes , OverlapRadii >::Index ] ) );
#endif // NEW_CODE
		for( int i=0 ; i<WindowSize< OverlapSizes >::Size ; i++ ) if( _isValidFEM1Node( _nodes[i] ) )
		{
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
			if( i!=WindowIndex< OverlapSizes , OverlapRadii >::Index ) row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( _nodes[i]->nodeData.nodeIndex-offset ) , (Real)( _values[i] + _stencil[i] ) );
#else // !NEW_CODE_SPARSE_MATRIX
			if( i!=WindowIndex< OverlapSizes , OverlapRadii >::Index ) row[count++] = MatrixEntry< Real >( (matrix_index_type)( _nodes[i]->nodeData.nodeIndex-offset ) , (Real)( _values[i] + _stencil[i] ) );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
			if( i!=WindowIndex< OverlapSizes , OverlapRadii >::Index ) row[count++] = MatrixEntry< Real >( _nodes[i]->nodeData.nodeIndex-offset , (Real)( _values[i] + _stencil[i] ) );
#endif // NEW_CODE
		}
	}
	else
	{
		LocalDepth d ; LocalOffset off;
		_localDepthAndOffset( node , d , off );
		Real temp = (Real)F.ccIntegrate( off , off ) + pointValues.data[ WindowIndex< OverlapSizes , OverlapRadii >::Index ];
	
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
		row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( nodeIndex-offset ) , temp );
#else // !NEW_CODE_SPARSE_MATRIX
		row[count++] = MatrixEntry< Real >( (matrix_index_type)( nodeIndex-offset ) , temp );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
		row[count++] = MatrixEntry< Real >( nodeIndex-offset , temp );
#endif // NEW_CODE
		LocalOffset _off;
		WindowLoop< Dim >::Run
		(
			ZeroUIntPack< Dim >() , OverlapSizes() ,
			[&]( int d , int i ){ _off[d] = off[d] - (int)OverlapRadii::Values[d] + i; } ,
			[&]( const FEMTreeNode* _node , Real pointValue )
			{
				if( node!=_node && FEMIntegrator::IsValidFEMNode( UIntPack< FEMSigs ... >() , d , _off ) )
				{
					Real temp = (Real)F.ccIntegrate( _off , off ) + pointValue;
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
					if( _isValidFEM1Node( _node ) ) row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( _node->nodeData.nodeIndex-offset ) , temp );
#else // !NEW_CODE_SPARSE_MATRIX
					if( _isValidFEM1Node( _node ) ) row[count++] = MatrixEntry< Real >( (matrix_index_type)( _node->nodeData.nodeIndex-offset ) , temp );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
					if( _isValidFEM1Node( _node ) ) row[count++] = MatrixEntry< Real >( _node->nodeData.nodeIndex-offset , temp );
#endif // NEW_CODE
				}
			} ,
			neighbors.neighbors() , pointValues()
		);
	}
	return constraint;
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , unsigned int PointD >
void FEMTree< Dim , Real >::_addProlongedPointValues( UIntPack< FEMSigs ... > , WindowSlice< Real , UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > > pointValues , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , const InterpolationInfo< T , PointD >* interpolationInfo ) const
{
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] This code is broken" )
#endif // SHOW_WARNINGS
#if 1
	ERROR_OUT( "Broken code" );
#else
	if( !interpolationInfo ) return;
	const InterpolationInfo< T , PointD >& iInfo = *interpolationInfo;
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	typedef UIntPack<    BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize    ... > OverlapSizes;

	const FEMTreeNode* node = neighbors.neighbors.data[ WindowIndex< OverlapSizes , OverlapRadii >::Index ];

	LocalDepth d , parentD ; LocalOffset off , parentOff;
	_localDepthAndOffset( node , d , off );
	_localDepthAndOffset( node->parent , parentD , parentOff );
	int fStart , fEnd;
	BSplineData< FEMSig >::FunctionSpan( d , fStart , fEnd );

	int fIdx[Dim];
	functionIndex( IsotropicUIntPack< Dim , FEMSig >() , node , fIdx );
	double       splineValues[ Dim ]               [ PointD+1 ];
	double parentSplineValues[ Dim ][ SupportSize ][ PointD+1 ];
	int s[Dim];
	CumulativeDerivativeValues< Real , Dim , PointD > dualValues;
	std::function< void ( const FEMTreeNode* , Real& ) > innerFunction = [&]( const FEMTreeNode* pNode , Real& pointValue )
	{
		if( _isValidFEM1Node( pNode ) )
		{
			CumulativeDerivativeValues< Real , Dim , PointD > values = Evaluate< SupportSize , Dim , Real , PointD >( s , parentSplineValues );
			pointValue += CumulativeDerivativeValues< Real , Dim , PointD >::Dot( dualValues , values );
		};
	};
	std::function< void ( const FEMTreeNode* ) > outerFunction = [&]( const FEMTreeNode* _node )
	{
		if( _isValidSpaceNode( _node ) ) for( const PointData< Dim , Real , T , PointD >* _pData=iInfo.begin( _node ) ; _pData!=iInfo.end( _node ) ; _pData++ )
		{
			// Evaluate the node's basis function at the sample
			const PointData< Dim , Real , T , PointD >& pData = *_pData;
			_setDValues< FEMSig , PointD , FEMDegree >( pData.position , _node , node , bsData , splineValues );
			_setDValues< FEMSig , PointD , FEMDegree >( pData.position , _node->parent , bsData , parentSplineValues );
			dualValues = iInfo.weights * Evaluate< Dim , Real , PointD >( splineValues ) * pData.weight;

			// Get the indices of the parent
			LocalDepth _parentD ; LocalOffset _parentOff;
			_localDepthAndOffset( _node->parent , _parentD , _parentOff );

			int _off[Dim];
			for( int dd=0 ; dd<Dim ; dd++ ) _off[dd] = _parentOff[dd] - parentOff[dd];

			int _start[Dim] , _end[Dim];
			for( int dd=0 ; dd<Dim ; dd++ ) _start[dd] = OverlapRadius + _off[dd] - LeftPointSupportRadius , _end[dd] = _start[dd] + SupportSize;
			WindowLoop< Dim >::Run
			(
				_start , _end ,
				[&]( int d , int i ){ s[d] = i + LeftPointSupportRadius - _off[d] - OverlapRadius; } ,
				innerFunction ,
				pNeighbors.neighbors() , pointValues
			);
		}
	};
	int start[Dim] , end[Dim];
	for( int dd=0 ; dd<Dim ; dd++ ) start[dd] = OverlapRadius - LeftSupportRadius , end[dd] = start[dd] + SupportSize;
	WindowLoop< Dim >::Run
	(
		start , end ,
		[&]( int , int ){;} ,
		outerFunction ,
		neighbors.neighbors()
	);
#endif
}

template< unsigned int Dim , class Real >
template< typename T , unsigned int ... PointDs , unsigned int ... FEMSigs >
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
int FEMTree< Dim , Real >::_setProlongedMatrixRow( const typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , Pointer( MatrixEntry< Real , matrix_index_type > ) row , node_index_type offset , const DynamicWindow< double , UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& stencil , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_CODE_SPARSE_MATRIX
int FEMTree< Dim , Real >::_setProlongedMatrixRow( const typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , Pointer( MatrixEntry< Real > ) row , node_index_type offset , const DynamicWindow< double , UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& stencil , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
int FEMTree< Dim , Real >::_setProlongedMatrixRow( const typename BaseFEMIntegrator::System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , Pointer( MatrixEntry< Real > ) row , int offset , const DynamicWindow< double , UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& stencil , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_CODE
{
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	typedef UIntPack<    BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize    ... > OverlapSizes;
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > FEMDegrees;

	int count = 0;
	const FEMTreeNode* node = neighbors.neighbors.data[ WindowIndex< OverlapSizes , OverlapRadii >::Index ];
	LocalDepth d , parentD ; LocalOffset off , parentOff;
	_localDepthAndOffset( node , d , off );
	_localDepthAndOffset( node->parent , parentD , parentOff );
	bool isInterior = _isInteriorlyOverlapped( FEMDegrees() , FEMDegrees() , node->parent );
	StaticWindow< Real , UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > > pointValues;
	memset( pointValues.data , 0 , sizeof(Real)*WindowSize< OverlapSizes >::Size );
	_addProlongedPointValues( UIntPack< FEMSigs ... >() , pointValues() , neighbors , pNeighbors , bsData , interpolationInfo ... );

#ifdef NEW_CODE
	node_index_type nodeIndex = node->nodeData.nodeIndex;
#else // !NEW_CODE
	int nodeIndex = node->nodeData.nodeIndex;
#endif // NEW_CODE

	int start[Dim] , end[Dim];
	_SetParentOverlapBounds( FEMDegrees() , FEMDegrees() , node , start , end );
	if( isInterior ) // General case, so try to make fast
	{
		WindowLoop< Dim >::Run
		(
			start , end ,
			[&]( int , int ){;} ,
			[&]( const FEMTreeNode* node , const Real& pointValue , const Real& stencilValue )
			{
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
			if( _isValidFEM1Node( node ) ) row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( node->nodeData.nodeIndex - offset ) , pointValue + stencilValue );
#else // !NEW_CODE_SPARSE_MATRIX
			if( _isValidFEM1Node( node ) ) row[count++] = MatrixEntry< Real >( (matrix_index_type)( node->nodeData.nodeIndex - offset ) , pointValue + stencilValue );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
				if( _isValidFEM1Node( node ) ) row[count++] = MatrixEntry< Real >( node->nodeData.nodeIndex - offset , pointValue + stencilValue );
#endif // NEW_CODE
			} ,
			pNeighbors.neighbors() , pointValues() , stencil()
		);
	}
	else
	{
		LocalDepth d ; LocalOffset off;
		_localDepthAndOffset( node , d , off );
		WindowLoop< Dim >::Run
		(
			start , end , 
			[&]( int , int ){;} ,
			[&]( const FEMTreeNode* node , const Real& pointValue )
			{
				if( _isValidFEM1Node( node ) )
				{
					LocalDepth d ; LocalOffset _off;
					_localDepthAndOffset( node , d , _off );
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
					row[count++] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( node->nodeData.nodeIndex - offset ) , (Real)F.pcIntegrate( _off , off ) + pointValue );
#else // !NEW_CODE_SPARSE_MATRIX
					row[count++] = MatrixEntry< Real >( (matrix_index_type)( node->nodeData.nodeIndex - offset ) , (Real)F.pcIntegrate( _off , off ) + pointValue );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
					row[count++] = MatrixEntry< Real >( node->nodeData.nodeIndex - offset , (Real)F.pcIntegrate( _off , off ) + pointValue );
#endif // NEW_CODE
				}
			} ,
			pNeighbors.neighbors() , pointValues()
		);
	}
	return count;
}

template< unsigned int Dim , class Real >
template< unsigned int FEMDegree1 , unsigned int FEMDegree2 >
void FEMTree< Dim , Real >::_SetParentOverlapBounds( const FEMTreeNode* node , int start[Dim] , int end[Dim] )
{
	const int OverlapStart = BSplineOverlapSizes< FEMDegree1 , FEMDegree2 >::OverlapStart;

	if( node->parent )
	{
		int cIdx = (int)( node - node->parent->children );
		for( int d=0 ; d<Dim ; d++ )
		{
			start[d] = BSplineOverlapSizes< FEMDegree1 , FEMDegree2 >::ParentOverlapStart[ (cIdx>>d) & 1 ] - OverlapStart;
			end  [d] = BSplineOverlapSizes< FEMDegree1 , FEMDegree2 >::ParentOverlapEnd  [ (cIdx>>d) & 1 ] - OverlapStart + 1;
		}
	}
}
template< unsigned int Dim , class Real >
template< unsigned int FEMDegree1 , unsigned int FEMDegree2 >
void FEMTree< Dim , Real >::_SetParentOverlapBounds( int cIdx , int start[Dim] , int end[Dim] )
{
	const int OverlapStart = BSplineOverlapSizes< FEMDegree1 , FEMDegree2 >::OverlapStart;

	for( int d=0 ; d<Dim ; d++ )
	{
		start[d] = BSplineOverlapSizes< FEMDegree1 , FEMDegree2 >::ParentOverlapStart[ (cIdx>>d) & 1 ] - OverlapStart;
		end  [d] = BSplineOverlapSizes< FEMDegree1 , FEMDegree2 >::ParentOverlapEnd  [ (cIdx>>d) & 1 ] - OverlapStart + 1;
	}
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , unsigned int PointD >
T FEMTree< Dim , Real >::_getInterpolationConstraintFromProlongedSolution( const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , const FEMTreeNode* node , ConstPointer( T ) prolongedSolution , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , const InterpolationInfo< T , PointD >* interpolationInfo ) const
{
	if( !interpolationInfo ) return T();
	typedef UIntPack<    BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize    ... > SupportSizes;
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	typedef UIntPack< ( -BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportStart ) ... > LeftSupportRadii;
	typedef UIntPack<    BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize    ... > OverlapSizes;
	typedef PointEvaluatorState< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > > _PointEvaluatorState;
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	T temp = {};
	if( _isValidFEM1Node( node ) )
	{
		int s[Dim];
#if defined( _WIN32 ) || defined( _WIN64 )
#pragma message ( "[WARNING] You've got me MSVC" )
		auto  UpdateFunction = [&]( int d , int i ){ s[d] = (int)SupportSizes::Values[d] - 1 - ( i - (int)OverlapRadii::Values[d] + (int)LeftSupportRadii::Values[d] ); };
		auto ProcessFunction = [&]( const FEMTreeNode* pNode )
		{
			if( _isValidSpaceNode( pNode ) )
			{
				size_t begin , end;
				interpolationInfo->range( pNode , begin , end );
				for( size_t pIndex=begin ; pIndex<end ; pIndex++ )
				{
					const DualPointInfo< Dim , Real , T , PointD > _pData = (*interpolationInfo)[ pIndex ];
					_PointEvaluatorState peState;
					Point< Real , Dim > p = _pData.position;
					LocalDepth pD ; LocalOffset pOff;
					_localDepthAndOffset( pNode , pD , pOff );
					bsData.initEvaluationState( p , pD , pOff , peState );
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Why is this necessary?" )
#endif // SHOW_WARNINGS
					const int *_off = off;
					CumulativeDerivativeValues< Real , Dim , PointD > values = peState.template dValues< Real , CumulativeDerivatives< Dim , PointD > >( _off );
					for( int d=0 ; d<CumulativeDerivatives< Dim , PointD >::Size ; d++ ) temp += _pData.dualValues[d] * values[d];
				}
			}
		};
#endif // _WIN32 || _WIN64
		WindowLoop< Dim >::Run
		(
			OverlapRadii() - LeftSupportRadii() , OverlapRadii() - LeftSupportRadii() + SupportSizes() ,
#if defined( _WIN32 ) || defined( _WIN64 )
			UpdateFunction , ProcessFunction ,
#else // !_WIN32 && !_WIN64
			[&]( int d , int i ){ s[d] = (int)SupportSizes::Values[d] - 1 - ( i - (int)OverlapRadii::Values[d] + (int)LeftSupportRadii::Values[d] ); } ,
			[&]( const FEMTreeNode* pNode )
			{
				if( _isValidSpaceNode( pNode ) )
				{
					size_t begin , end;
					interpolationInfo->range( pNode , begin , end );
					for( size_t pIndex=begin ; pIndex<end ; pIndex++ )
					{
						const DualPointInfo< Dim , Real , T , PointD > _pData = (*interpolationInfo)[ pIndex ];
						_PointEvaluatorState peState;
						Point< Real , Dim > p = _pData.position;
						LocalDepth pD ; LocalOffset pOff;
						_localDepthAndOffset( pNode , pD , pOff );
						bsData.initEvaluationState( p , pD , pOff , peState );
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Why is this necessary?" )
#endif // SHOW_WARNINGS
						const int *_off = off;
						CumulativeDerivativeValues< Real , Dim , PointD > values = peState.template dValues< Real , CumulativeDerivatives< Dim , PointD > >( _off );
						for( int d=0 ; d<CumulativeDerivatives< Dim , PointD >::Size ; d++ ) temp += _pData.dualValues[d] * values[d];
					}
				}
			} ,
#endif // _WIN32 || _WIN64
			neighbors.neighbors()
		);
	}
	return temp;
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , unsigned int ... PointDs >
T FEMTree< Dim , Real >::_getConstraintFromProlongedSolution( UIntPack< FEMSigs ... > , const BaseSystem< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& neighbors , const typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& pNeighbors , const FEMTreeNode* node , ConstPointer( T ) prolongedSolution , const DynamicWindow< double , UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >& stencil , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
{
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > FEMDegrees;

	if( _localDepth( node )<=0 ) return T();
	// This is a conservative estimate as we only need to make sure that the parent nodes don't overlap the child (not the parent itself)
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node->parent , d , off );
	bool isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( FEMDegrees() , FEMDegrees() , d , off );

	// Offset the constraints using the solution from lower resolutions.
	T constraint = {};
	static const WindowLoopData< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > > loopData( []( int c , int* start , int* end ){ BaseFEMIntegrator::ParentOverlapBounds( FEMDegrees() , FEMDegrees() , c , start , end ); } );
	int cIdx = (int)( node - node->parent->children );
	unsigned int size = loopData.size[cIdx];
	const unsigned int* indices = loopData.indices[cIdx];
	ConstPointer( double ) values = stencil.data;
	ConstPointer( FEMTreeNode * const ) nodes = pNeighbors.neighbors().data;
	if( isInterior )
	{
		for( unsigned int i=0 ; i<size ; i++ )
		{
			unsigned int idx = indices[i];
			if( _isValidFEM1Node( nodes[idx] ) ) constraint += (T)( prolongedSolution[ nodes[idx]->nodeData.nodeIndex ] * (Real)values[idx] );
		}
	}
	else
	{
		LocalDepth d ; LocalOffset off;
		_localDepthAndOffset( node , d , off );
		for( unsigned int i=0 ; i<size ; i++ )
		{
			unsigned int idx = indices[i];
			if( _isValidFEM1Node( nodes[idx] ) )
			{
				LocalDepth _d ; LocalOffset _off;
				_localDepthAndOffset( nodes[idx] , _d , _off );
				constraint += (T)( prolongedSolution[ nodes[idx]->nodeData.nodeIndex ] * (Real)F.pcIntegrate( _off , off ) );
			}

		}
	}
	return constraint + _getInterpolationConstraintFromProlongedSolution( neighbors , node , prolongedSolution , bsData , interpolationInfo... );
}

// Given the solution @( depth ) add to the met constraints @( depth-1 )
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_updateRestrictedIntegralConstraints( UIntPack< FEMSigs ... > , ThreadPool &tp , const typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth highDepth , ConstPointer( T ) fineSolution , Pointer( T ) restrictedConstraints ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_updateRestrictedIntegralConstraints( UIntPack< FEMSigs ... > , const typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth highDepth , ConstPointer( T ) fineSolution , Pointer( T ) restrictedConstraints ) const
#endif // NEW_THREADS
{
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > FEMDegrees;
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	typedef UIntPack<    BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize    ... > OverlapSizes;

	if( highDepth<=0 ) return;
	// Get the stencil describing the Laplacian relating coefficients @(highDepth) with coefficients @(highDepth-1)
	PCStencils< FEMDegrees > stencils;
	F.template setStencils< true >(  stencils );
#ifdef NEW_CODE
	node_index_type start = _sNodesBegin(highDepth) , end = _sNodesEnd(highDepth);
	node_index_type range = end-start;
	node_index_type lStart = _sNodesBegin(highDepth-1);
#else // !NEW_CODE
	size_t start = _sNodesBegin(highDepth) , end = _sNodesEnd(highDepth) , range = end-start;
	int lStart = _sNodesBegin(highDepth-1);
#endif // NEW_CODE

	// Iterate over the nodes @(highDepth)
#ifdef NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( highDepth )-1 );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(highDepth) , _sNodesEnd(highDepth) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(highDepth) , _sNodesEnd(highDepth) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(highDepth) ; i<_sNodesEnd(highDepth) ; i++ )
#else // !NEW_CODE
	for( int i=_sNodesBegin(highDepth) ; i<_sNodesEnd(highDepth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
		{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			FEMTreeNode* node = _sNodes.treeNodes[i];

			// Offset the coarser constraints using the solution from the current resolutions.
			int cIdx = (int)( node - node->parent->children );

			{
				typename FEMTreeNode::template ConstNeighbors< OverlapSizes > pNeighbors;
				neighborKey.getNeighbors( OverlapRadii() , OverlapRadii() , node->parent , pNeighbors );
				const DynamicWindow< double , OverlapSizes >& stencil = stencils.data[cIdx];

				bool isInterior = _isInteriorlyOverlapped( FEMDegrees() , FEMDegrees() , node->parent );
				LocalDepth d ; LocalOffset off;
				_localDepthAndOffset( node , d , off );

				// Offset the constraints using the solution from finer resolutions.
				int start[Dim] , end[Dim];
				_SetParentOverlapBounds( FEMDegrees() , FEMDegrees() , node , start , end );

				T solution = fineSolution[ node->nodeData.nodeIndex ];
				ConstPointer( FEMTreeNode * const ) nodes = pNeighbors.neighbors().data;
				ConstPointer( double ) stencilValues = stencil.data;
				if( isInterior )
				{
					for( int i=0 ; i<WindowSize< OverlapSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) )
					{
						AddAtomic( restrictedConstraints[ nodes[i]->nodeData.nodeIndex ] , solution * (Real)stencilValues[i] );
					}
				}
				else
				{
					for( int i=0 ; i<WindowSize< OverlapSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) )
					{
						LocalDepth _d ; LocalOffset _off;
						_localDepthAndOffset( nodes[i] , _d , _off );
						AddAtomic( restrictedConstraints[ nodes[i]->nodeData.nodeIndex ] , solution * (Real)F.pcIntegrate( _off , off ) );
					}
				}
			}
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , unsigned int PointD >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_setPointValuesFromProlongedSolution( ThreadPool &tp , LocalDepth highDepth , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , ConstPointer( T ) prolongedSolution , InterpolationInfo< T , PointD >* iInfo ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_setPointValuesFromProlongedSolution( LocalDepth highDepth , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , ConstPointer( T ) prolongedSolution , InterpolationInfo< T , PointD >* iInfo ) const
#endif // NEW_THREADS
{
	if( !iInfo ) return;
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > FEMDegrees;
	InterpolationInfo< T , PointD >& interpolationInfo = *iInfo;

	LocalDepth lowDepth = highDepth-1;
	if( lowDepth<0 ) return;
	// For every node at the current depth
#ifdef NEW_THREADS
	std::vector< ConstPointSupportKey< FEMDegrees > > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< ConstPointSupportKey< FEMDegrees > > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( lowDepth ) );

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(highDepth) , _sNodesEnd(highDepth) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(highDepth) , _sNodesEnd(highDepth) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(highDepth) ; i<_sNodesEnd(highDepth) ; i++ )
#else // !NEW_CODE
	for( int i=_sNodesBegin(highDepth) ; i<_sNodesEnd(highDepth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
		{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			ConstPointSupportKey< FEMDegrees >& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			ConstPointSupportKey< FEMDegrees >& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			ConstPointSupportKey< FEMDegrees >& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			if( _isValidSpaceNode( _sNodes.treeNodes[i] ) )
			{
				size_t begin , end;
				interpolationInfo.range( _sNodes.treeNodes[i] , begin , end );
				for( size_t pIndex=begin ; pIndex<end ; pIndex++ )
				{
					DualPointInfo< Dim , Real , T , PointD >& pData = interpolationInfo[ pIndex ];
					neighborKey.getNeighbors( _sNodes.treeNodes[i]->parent );
#ifdef _MSC_VER
					pData.dualValues = interpolationInfo( pIndex , _coarserFunctionValues< PointD , T , FEMSigs ... >( UIntPack< FEMSigs ... >() , pData.position , neighborKey , _sNodes.treeNodes[i] , bsData , prolongedSolution ) ) * pData.weight;
#else // !_MSC_VER
					pData.dualValues = interpolationInfo( pIndex , _coarserFunctionValues< PointD >( UIntPack< FEMSigs ... >() , pData.position , neighborKey , _sNodes.treeNodes[i] , bsData , prolongedSolution ) ) * pData.weight;
#endif // _MSC_VER
				}
			}
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , unsigned int PointD >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_updateRestrictedInterpolationConstraints( ThreadPool &tp , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth highDepth , ConstPointer( T ) solution , Pointer( T ) restrictedConstraints , const InterpolationInfo< T , PointD >* iInfo ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_updateRestrictedInterpolationConstraints( const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth highDepth , ConstPointer( T ) solution , Pointer( T ) restrictedConstraints , const InterpolationInfo< T , PointD >* iInfo ) const
#endif // NEW_THREADS
{
	if( !iInfo ) return;
	const InterpolationInfo< T , PointD >& interpolationInfo = *iInfo;
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > FEMDegrees;
	typedef UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > SupportSizes;

	// Note: We can't iterate over the finer point nodes as the point weights might be
	// scaled incorrectly, due to the adaptive exponent. So instead, we will iterate
	// over the coarser nodes and evaluate the finer solution at the associated points.
	LocalDepth  lowDepth = highDepth-1;
	if( lowDepth<0 ) return;

#ifdef NEW_CODE
	node_index_type start = _sNodesBegin(lowDepth) , end = _sNodesEnd(lowDepth);
#else // !NEW_CODE
	size_t start = _sNodesBegin(lowDepth) , end = _sNodesEnd(lowDepth);
#endif // NEW_CODE
#ifdef NEW_THREADS
	std::vector< ConstPointSupportKey< FEMDegrees > > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< ConstPointSupportKey< FEMDegrees > > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( lowDepth ) );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(lowDepth) , _sNodesEnd(lowDepth) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(lowDepth) , _sNodesEnd(lowDepth) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(lowDepth) ; i<_sNodesEnd(lowDepth) ; i++ )
#else // !NEW_CODE
	for( int i=_sNodesBegin(lowDepth) ; i<_sNodesEnd(lowDepth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidSpaceNode( _sNodes.treeNodes[i] ) ) if( _isValidSpaceNode( _sNodes.treeNodes[i] ) )
		{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			ConstPointSupportKey< FEMDegrees >& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			ConstPointSupportKey< FEMDegrees >& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			ConstPointSupportKey< FEMDegrees >& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			PointEvaluatorState< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > > peState;
			const FEMTreeNode* node = _sNodes.treeNodes[i];

			LocalDepth d ; LocalOffset off;
			_localDepthAndOffset( node , d , off );
			typename FEMTreeNode::template ConstNeighbors< SupportSizes >& neighbors = neighborKey.getNeighbors( node );
			size_t begin , end;
			interpolationInfo.range( node , begin , end );
			for( size_t pIndex=begin ; pIndex<end ; pIndex++ )
			{
				const DualPointInfo< Dim , Real , T , PointD >& pData = interpolationInfo[ pIndex ];
				Point< Real , Dim > p = pData.position;
				bsData.initEvaluationState( p , d , off , peState );

#ifdef _MSC_VER
				CumulativeDerivativeValues< T , Dim , PointD > dualValues = interpolationInfo( pIndex , _finerFunctionValues< PointD , T , FEMSigs ... >( UIntPack< FEMSigs ... >() , pData.position , neighborKey , node , bsData , solution ) ) * pData.weight;
#else // !_MSC_VER
				CumulativeDerivativeValues< T , Dim , PointD > dualValues = interpolationInfo( pIndex , _finerFunctionValues< PointD >( UIntPack< FEMSigs ... >() , pData.position , neighborKey , node , bsData , solution ) ) * pData.weight;
#endif // _MSC_VER
				// Update constraints for all nodes @( depth-1 ) that overlap the point
				int s[Dim];
				WindowLoop< Dim >::Run
				(
					ZeroUIntPack< Dim >() , SupportSizes() ,
					[&]( int d , int i ){ s[d] = i; } ,
					[&]( const FEMTreeNode* node )
				{
					if( _isValidFEM1Node( node ) )
					{
						LocalDepth d ; LocalOffset off;
						_localDepthAndOffset( node , d , off );
						CumulativeDerivativeValues< Real , Dim , PointD > values = peState.template dValues< Real , CumulativeDerivatives< Dim , PointD > >( off );
						T temp = {};
						for( int d=0 ; d<CumulativeDerivatives< Dim , PointD >::Size ; d++ ) temp += dualValues[d] * values[d];
						AddAtomic( restrictedConstraints[ node->nodeData.nodeIndex ] , temp );
					}
				} ,
					neighbors.neighbors()
					);
			}
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< unsigned int Dim , class Real >
template< class C , unsigned int ... FEMSigs >
#ifdef NEW_THREADS
DenseNodeData< C , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::coarseCoefficients( ThreadPool &tp , const DenseNodeData< C , UIntPack< FEMSigs ... > >& coefficients ) const
#else // !NEW_THREADS
DenseNodeData< C , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::coarseCoefficients( const DenseNodeData< C , UIntPack< FEMSigs ... > >& coefficients ) const
#endif // NEW_THREADS
{
	DenseNodeData< C , UIntPack< FEMSigs ... > > coarseCoefficients( _sNodesEnd(_maxDepth-1) );
	memset( coarseCoefficients() , 0 , sizeof(Real)*_sNodesEnd(_maxDepth-1) );
#ifdef NEW_CODE
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(0) , _sNodesEnd(_maxDepth-1) , [&]( const ThreadPool::ThreadNum & , size_t i ){ coarseCoefficients[i] = coefficients[i]; } );
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(0) , _sNodesEnd(_maxDepth-1) , [&]( unsigned int , size_t i ){ coarseCoefficients[i] = coefficients[i]; } );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
	for( node_index_type i=_sNodesBegin(0) ; i<_sNodesEnd(_maxDepth-1) ; i++ ) coarseCoefficients[i] = coefficients[i];
#endif //NEW_THREADS
#else // !NEW_CODE
#pragma omp parallel for
	for( int i=_sNodesBegin(0) ; i<_sNodesEnd(_maxDepth-1) ; i++ ) coarseCoefficients[i] = coefficients[i];
#endif // NEW_CODE
	typename FEMIntegrator::template RestrictionProlongation< UIntPack< FEMSigs ... > > rp;
#ifdef NEW_THREADS
	for( LocalDepth d=1 ; d<_maxDepth ; d++ ) _upSample( UIntPack< FEMSigs ... >() , tp , rp , d , coarseCoefficients() );
#else // !NEW_THREADS
	for( LocalDepth d=1 ; d<_maxDepth ; d++ ) _upSample( UIntPack< FEMSigs ... >() , rp , d , coarseCoefficients() );
#endif // NEW_THREADS
	return coarseCoefficients;
}

template< unsigned int Dim , class Real >
template< class C , unsigned int ... FEMSigs >
#ifdef NEW_THREADS
DenseNodeData< C , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::coarseCoefficients( ThreadPool &tp , const SparseNodeData< C , UIntPack< FEMSigs ... > >& coefficients ) const
#else // !NEW_THREADS
DenseNodeData< C , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::coarseCoefficients( const SparseNodeData< C , UIntPack< FEMSigs ... > >& coefficients ) const
#endif // NEW_THREADS
{
	DenseNodeData< C , UIntPack< FEMSigs ... > > coarseCoefficients( _sNodesEnd(_maxDepth-1) );
	memset( coarseCoefficients() , 0 , sizeof(C)*_sNodesEnd(_maxDepth-1) );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(0) , _sNodesEnd(_maxDepth-1) , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(0) , _sNodesEnd(_maxDepth-1) , [&]( unsigned int , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(0) ; i<_sNodesEnd(_maxDepth-1) ; i++ )
#else // !NEW_CODE
	for( int i=_sNodesBegin(0) ; i<_sNodesEnd(_maxDepth-1) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		const C* c = coefficients( _sNodes.treeNodes[i] );
		if( c ) coarseCoefficients[i] = *c;
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	typename FEMIntegrator::template RestrictionProlongation< UIntPack< FEMSigs ... > > rp;
#ifdef NEW_THREADS
	for( LocalDepth d=1 ; d<_maxDepth ; d++ ) _upSample( UIntPack< FEMSigs ... >() , tp , rp , d , coarseCoefficients() );
#else // !NEW_THREADS
	for( LocalDepth d=1 ; d<_maxDepth ; d++ ) _upSample( UIntPack< FEMSigs ... >() , rp , d , coarseCoefficients() );
#endif // NEW_THREADS
	return coarseCoefficients;
}

template< unsigned int Dim , class Real >
template< unsigned int PointD , typename T , unsigned int ... FEMSigs >
CumulativeDerivativeValues< T , Dim , PointD > FEMTree< Dim , Real >::_coarserFunctionValues( UIntPack< FEMSigs ... > , Point< Real , Dim > p , const ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey , const FEMTreeNode* pointNode , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , ConstPointer( T ) solution ) const
{
	typedef UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > SupportSizes;

	CumulativeDerivativeValues< T , Dim , PointD > values;
	LocalDepth depth = _localDepth( pointNode );
	if( depth<0 ) return values;
	// Iterate over all basis functions that overlap the point at the coarser resolutions
	{

		PointEvaluatorState< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > > peState;
		LocalDepth _d ; LocalOffset _off;
		_localDepthAndOffset( pointNode->parent , _d , _off );
		bsData.initEvaluationState( p , _d , _off , peState );
		const typename FEMTreeNode::template ConstNeighbors< SupportSizes >& neighbors = neighborKey.neighbors[ _localToGlobal( depth-1 ) ];
		ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;

		for( unsigned int i=0 ; i<WindowSize< SupportSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) )
		{
			LocalDepth d ; LocalOffset off;
			_localDepthAndOffset( nodes[i] , d , off );
			CumulativeDerivativeValues< Real , Dim , PointD > temp = peState.template dValues< Real , CumulativeDerivatives< Dim , PointD > >( off );
			const T& _solution = solution[ nodes[i]->nodeData.nodeIndex ];
			for( int s=0 ; s<CumulativeDerivatives< Dim , PointD >::Size ; s++ ) values[s] += _solution * temp[s];
		}
	}
	return values;
}

template< unsigned int Dim , class Real >
template< unsigned int PointD , typename T , unsigned int ... FEMSigs >
CumulativeDerivativeValues< T , Dim , PointD > FEMTree< Dim , Real >::_finerFunctionValues( UIntPack< FEMSigs ... > , Point< Real , Dim > p , const ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey , const FEMTreeNode* pointNode , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , ConstPointer( T ) solution ) const
{
	typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > > childNeighbors;
	typedef UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree  >::SupportSize ... > SupportSizes;

	CumulativeDerivativeValues< T , Dim , PointD > values;
	LocalDepth depth = _localDepth( pointNode );
	neighborKey.getChildNeighbors( _childIndex( pointNode , p ) , _localToGlobal( depth ) , childNeighbors );
	PointEvaluatorState< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > > peState;
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( pointNode , d , off );
	int cIdx = _childIndex( pointNode , p );
	d++;
	for( int dd=0 ; dd<Dim ; dd++ ) off[dd] = (off[dd]<<1) | ( (cIdx>>dd) & 1 );
	bsData.initEvaluationState( p , d , off , peState );
	int s[Dim];
	WindowLoop< Dim >::Run
	(
		ZeroUIntPack< Dim >() , SupportSizes() ,
		[&]( int d , int i ){ s[d] = i; } ,
		[&]( const FEMTreeNode* node )
		{
			if( _isValidFEM1Node( node ) )
			{
				LocalDepth d ; LocalOffset off;
				_localDepthAndOffset( node , d , off );
				CumulativeDerivativeValues< Real , Dim , PointD > dValues = peState.template dValues< Real , CumulativeDerivatives< Dim , PointD > >( off );
				const T& _solution = solution[ node->nodeData.nodeIndex ];
				for( int s=0 ; s<CumulativeDerivatives< Dim , PointD >::Size ; s++ ) values[s] += _solution * dValues[s];
			}
		} ,
		childNeighbors.neighbors()
	);
	return values;
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , unsigned int ... PointDs >
#ifdef NEW_CODE
#ifdef NEW_THREADS
int FEMTree< Dim , Real >::_getSliceMatrixAndProlongationConstraints( UIntPack< FEMSigs ... > , ThreadPool &tp , const typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size >& matrix , Pointer( Real ) diagonalR , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , node_index_type nBegin , node_index_type nEnd , ConstPointer( T ) prolongedSolution , Pointer( T ) constraints , const CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& ccStencil , const PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& pcStencils , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
int FEMTree< Dim , Real >::_getSliceMatrixAndProlongationConstraints( UIntPack< FEMSigs ... > , const typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , SparseMatrix< Real , matrix_index_type , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size >& matrix , Pointer( Real ) diagonalR , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , node_index_type nBegin , node_index_type nEnd , ConstPointer( T ) prolongedSolution , Pointer( T ) constraints , const CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& ccStencil , const PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& pcStencils , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_THREADS
#else // !NEW_CODE
int FEMTree< Dim , Real >::_getSliceMatrixAndProlongationConstraints( UIntPack< FEMSigs ... > , const typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , SparseMatrix< Real , int , WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > >::Size >& matrix , Pointer( Real ) diagonalR , const PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bsData , LocalDepth depth , int nBegin , int nEnd , ConstPointer( T ) prolongedSolution , Pointer( T ) constraints , const CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& ccStencil , const PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& pcStencils , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_CODE
{
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > FEMDegrees;
	typedef UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize ... > OverlapSizes;
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	size_t range = nEnd - nBegin;
#ifdef NEW_CODE
	matrix.resize( range );
#else // !NEW_CODE
	matrix.resize( (int)range );
#endif // NEW_CODE
#ifdef NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( depth ) );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( 0 , range , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( 0 , range , [&]( unsigned int thread , size_t i )
#endif //NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=0 ; i<(node_index_type)range ; i++ )
#else // !NEW_CODE
	for( int i=0 ; i<(int)range ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidFEM1Node( _sNodes.treeNodes[i+nBegin] ) )
		{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			FEMTreeNode* node = _sNodes.treeNodes[i+nBegin];
			// Get the matrix row size	
			typename FEMTreeNode::template ConstNeighbors< OverlapSizes > neighbors , pNeighbors;
			neighborKey.getNeighbors( OverlapRadii() , OverlapRadii() , node , pNeighbors , neighbors );
			// Set the row entries
			if( constraints ) constraints[i] = _setMatrixRowAndGetConstraintFromProlongation( UIntPack< FEMSigs ... >() , F , pNeighbors , neighbors , i , matrix , nBegin , pcStencils , ccStencil , bsData , prolongedSolution , interpolationInfo... );
			else                               _setMatrixRowAndGetConstraintFromProlongation( UIntPack< FEMSigs ... >() , F , pNeighbors , neighbors , i , matrix , nBegin , pcStencils , ccStencil , bsData , prolongedSolution , interpolationInfo... );
			if( diagonalR ) diagonalR[i] = (Real)1. / matrix[i][0].Value;
		}
		else if( constraints ) constraints[i] = T();
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Why do we care if the node is not valid?" )
#endif // SHOW_WARNINGS
#if !defined( _WIN32 ) && !defined( _WIN64 )
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] I'm not sure how expensive this system call is on non-Windows system. (You may want to comment this out.)" )
#endif // SHOW_WARNINGS
#endif // !_WIN32 && !_WIN64
	MemoryUsage();
#ifdef NEW_THREADS
#endif // NEW_THREADS
	return 1;
}

template< unsigned int Dim , class Real >
template< typename T , unsigned int ... PointDs , unsigned int ... FEMSigs >
#ifdef NEW_CODE
#ifdef NEW_THREADS
SparseMatrix< Real , matrix_index_type > FEMTree< Dim , Real >::systemMatrix( UIntPack< FEMSigs ... > , ThreadPool &tp , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth depth , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
SparseMatrix< Real , matrix_index_type > FEMTree< Dim , Real >::systemMatrix( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth depth , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif //  NEW_THREADS
#else // !NEW_CODE
SparseMatrix< Real , int > FEMTree< Dim , Real >::systemMatrix( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth depth , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_CODE
{
	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	typedef typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > > BaseSystem;
	if( depth<0 || depth>_maxDepth ) ERROR_OUT( "System depth out of bounds: 0 <= " , depth , " <= " , _maxDepth );
#ifdef NEW_CODE
	SparseMatrix< Real , matrix_index_type > matrix;
#else // !NEW_CODE
	SparseMatrix< Real , int > matrix;
#endif // NEW_CODE
	F.init( depth );
	PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > > bsData( depth );

	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	typedef UIntPack<    BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize    ... > OverlapSizes;

	CCStencil< UIntPack< FEMSignature< FEMSigs >::Degree ... > > stencil;
	PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > > stencils;

	F.template setStencil< false >( stencil );

	matrix.resize( _sNodesSize(depth) );
#ifdef NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( depth ) );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd( depth ) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd( depth ) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(depth) ; i<_sNodesEnd( depth ) ; i++ )
#else // !NEW_CODE
	for( int i=_sNodesBegin(depth) ; i<_sNodesEnd( depth ) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
		{
#ifdef NEW_CODE
			node_index_type ii = i - _sNodesBegin(depth);
#else // !NEW_CODE
			int ii = i - _sNodesBegin(depth);
#endif // NEW_CODE
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			typename FEMTreeNode::template ConstNeighbors< OverlapSizes > neighbors;
			neighborKey.getNeighbors( OverlapRadii() , OverlapRadii() , _sNodes.treeNodes[i] , neighbors );

#if defined( __GNUC__ ) && __GNUC__ < 5
			#warning "you've got me gcc version<5"
				matrix.setRowSize( ii , _getMatrixRowSize( UIntPack< FEMSigs ... >() , neighbors ) );
#else // !__GNUC__ || __GNUC__ >=5
			matrix.setRowSize( ii , _getMatrixRowSize< FEMSigs ... >( neighbors ) );
#endif // __GNUC__ || __GNUC__ < 4
			_setMatrixRowAndGetConstraintFromProlongation( UIntPack< FEMSigs ... >() , F ,  neighbors , neighbors , matrix[ii] , _sNodesBegin(depth) , stencils , stencil , bsData , ( ConstPointer( T ) )NullPointer( T ) , interpolationInfo ... );
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	return matrix;
}

template< unsigned int Dim , class Real >
template< typename T , unsigned int ... PointDs , unsigned int ... FEMSigs >
#ifdef NEW_CODE
#ifdef NEW_THREADS
SparseMatrix< Real , matrix_index_type > FEMTree< Dim , Real >::prolongedSystemMatrix( UIntPack< FEMSigs ... > , ThreadPool &tp , typename BaseFEMIntegrator::template System< UIntPack<FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth highDepth , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
SparseMatrix< Real , matrix_index_type > FEMTree< Dim , Real >::prolongedSystemMatrix( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template System< UIntPack<FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth highDepth , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_THREADS
#else // !NEW_CODE
SparseMatrix< Real , int > FEMTree< Dim , Real >::prolongedSystemMatrix( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template System< UIntPack<FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth highDepth , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_CODE
{
	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	if( highDepth<=0 || highDepth>_maxDepth ) ERROR_OUT( "System depth out of bounds: 0 < " , highDepth , " <= " , _maxDepth );

	LocalDepth lowDepth = highDepth-1;
#ifdef NEW_CODE
	SparseMatrix< Real , matrix_index_type > matrix;
#else // !NEW_CODE
	SparseMatrix< Real , int > matrix;
#endif // NEW_CODE
	F.init( highDepth );
	PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > > bsData( highDepth );
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapStart ) ... > OverlapRadii;
	typedef UIntPack<    BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree >::OverlapSize    ... > OverlapSizes;

	PCStencils< UIntPack< FEMSignature< FEMSigs >::Degree ... > > stencils;
	F.template setStencils< true >( stencils );

	matrix.resize( _sNodesSize(highDepth) );
#ifdef NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( highDepth ) );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(highDepth) , _sNodesEnd(highDepth) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(highDepth) , _sNodesEnd(highDepth) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(highDepth) ; i<_sNodesEnd(highDepth) ; i++ )
#else // !NEW_CODE
	for( int i=_sNodesBegin(highDepth) ; i<_sNodesEnd(highDepth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
		{
#ifdef NEW_CODE
			node_index_type ii = i - _sNodesBegin(highDepth);
#else // !NEW_CODE
			int ii = i - _sNodesBegin(highDepth);
#endif // NEW_CODE
			int cIdx = (int)( _sNodes.treeNodes[i]-_sNodes.treeNodes[i]->parent->children );

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			typename FEMTreeNode::template ConstNeighbors< OverlapSizes > neighbors , pNeighbors;
			neighborKey.getNeighbors( OverlapRadii() , OverlapRadii() , _sNodes.treeNodes[i] , neighbors );
			neighborKey.getNeighbors( OverlapRadii() , OverlapRadii() , _sNodes.treeNodes[i]->parent , pNeighbors );

			matrix.setRowSize( ii , _getProlongedMatrixRowSize< FEMSigs ... >( _sNodes.treeNodes[i] , pNeighbors ) );
			_setProlongedMatrixRow< Real , PointDs ... >( F , neighbors , pNeighbors , matrix[ii] , _sNodesBegin(lowDepth) , stencils.data[cIdx] , bsData , interpolationInfo... );
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	return matrix;
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
#ifdef NEW_CODE
#ifdef NEW_THREADS
SparseMatrix< Real , matrix_index_type > FEMTree< Dim , Real >::downSampleMatrix( UIntPack< FEMSigs ... > , ThreadPool &tp , LocalDepth highDepth ) const
#else // !NEW_THREADS
SparseMatrix< Real , matrix_index_type > FEMTree< Dim , Real >::downSampleMatrix( UIntPack< FEMSigs ... > , LocalDepth highDepth ) const
#endif // NEW_THREADS
#else // !NEW_CODE
SparseMatrix< Real , int > FEMTree< Dim , Real >::downSampleMatrix( UIntPack< FEMSigs ... > , LocalDepth highDepth ) const
#endif // NEW_CODE
{
#ifdef NEW_CODE
	SparseMatrix< Real , matrix_index_type > matrix;
#else // !NEW_CODE
	SparseMatrix< Real , int > matrix;
#endif // NEW_CODE
	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > FEMDegrees;
	typedef UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleSize ... > UpSampleSizes;
	typedef IntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleStart ... > UpSampleStarts;
	typedef typename FEMTreeNode::template ConstNeighborKey< UIntPack< -BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleStart ... > , UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleEnd ... > > UpSampleKey;

	LocalDepth lowDepth = highDepth-1;
	if( lowDepth<0 ) return matrix;

	matrix.resize( _sNodesSize( lowDepth ) );

	typename EvaluationData::UpSampleEvaluator* upSampleEvaluators[] = { new typename BSplineEvaluationData< FEMSigs >::UpSampleEvaluator() ... };
	for( int d=0 ; d<Dim ; d++ ) upSampleEvaluators[d]->set( lowDepth );
#ifdef NEW_THREADS
	std::vector< UpSampleKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< UpSampleKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( lowDepth ) );

	DynamicWindow< double , UpSampleSizes > upSampleStencil;
	int lowCenter = ( 1<<lowDepth )>>1;
	double value[Dim+1] ; value[0] = 1;
	WindowLoop< Dim >::Run
	(
		ZeroUIntPack< Dim >() , UpSampleSizes() ,
		[&]( int d , int i ){ value[d+1] = value[d] * upSampleEvaluators[d]->value( lowCenter , 2*lowCenter + i + UpSampleStarts::Values[d] ); } ,
		[&]( double&  stencilValue ){ stencilValue = value[Dim]; } ,
		upSampleStencil()
	);

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(lowDepth) , _sNodesEnd(lowDepth) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(lowDepth) , _sNodesEnd(lowDepth) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(lowDepth) ; i<_sNodesEnd(lowDepth) ; i++ )
#else // !NEW_CODE
	for( int i=_sNodesBegin(lowDepth) ; i<_sNodesEnd(lowDepth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
		{
#ifdef NEW_CODE
			node_index_type _i = i - _sNodesBegin(lowDepth);
#else // !NEW_CODE
			int _i = i - _sNodesBegin(lowDepth);
#endif // NEW_CODE
			FEMTreeNode* pNode = _sNodes.treeNodes[i];

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			UpSampleKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			UpSampleKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			UpSampleKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			LocalDepth d ; LocalOffset off;
			_localDepthAndOffset( pNode , d , off );
			neighborKey.getNeighbors( pNode );
			typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleSize ... > > neighbors;
			neighborKey.getChildNeighbors( 0 , _localToGlobal( d ) , neighbors );


			int rowSize = 0;
			ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
			for( int i=0 ; i<WindowSize< UpSampleSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) ) rowSize++;

			matrix.setRowSize( _i , rowSize );
			matrix.rowSizes[_i] = 0;

			// Want to make sure test if contained children are interior.
			// This is more conservative because we are test that overlapping children are interior
			bool isInterior = _isInteriorlyOverlapped( FEMDegrees() , FEMDegrees() , pNode );

			if( isInterior )
			{
				ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
				ConstPointer( double ) stencilValues = upSampleStencil().data;
				for( int i=0 ; i<WindowSize< UpSampleSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) )
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
					matrix[_i][ matrix.rowSizes[_i]++ ] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( nodes[i]->nodeData.nodeIndex - _sNodesBegin(highDepth) ) , (Real)stencilValues[i] );
#else // !NEW_CODE_SPARSE_MATRIX
					matrix[_i][ matrix.rowSizes[_i]++ ] = MatrixEntry< Real >( (matrix_index_type)( nodes[i]->nodeData.nodeIndex - _sNodesBegin(highDepth) ) , (Real)stencilValues[i] );
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
					matrix[_i][ matrix.rowSizes[_i]++ ] = MatrixEntry< Real >( nodes[i]->nodeData.nodeIndex - _sNodesBegin(highDepth) , (Real)stencilValues[i] );
#endif // NEW_CODE
			}
			else
			{
				double upSampleValues[Dim][ UpSampleSizes::Max() ];

				WindowLoop< Dim >::Run
				(
					ZeroUIntPack< Dim >() , UpSampleSizes() ,
					[&]( int d , int i ){ upSampleValues[d][i] = upSampleEvaluators[d]->value( off[d] , 2*off[d] + i + UpSampleStarts::Values[d] ); } ,
					[&]( void ){}
				);

				double values[Dim+1] ; values[0] = 1;
				WindowLoop< Dim , Dim >::Run
				(
					ZeroUIntPack< Dim >() , UpSampleSizes() ,
					[&]( int d , int i ){ values[d+1] = values[d] * upSampleValues[d][i]; } ,
#ifdef NEW_CODE
#ifdef NEW_CODE_SPARSE_MATRIX
					[&]( const FEMTreeNode* node ){ if( _isValidFEM1Node( node ) ) matrix[_i][ matrix.rowSizes[_i]++ ] = MatrixEntry< Real , matrix_index_type >( (matrix_index_type)( node->nodeData.nodeIndex - _sNodesBegin(highDepth) ) , (Real)values[Dim] ); } ,
#else // !NEW_CODE_SPARSE_MATRIX
					[&]( const FEMTreeNode* node ){ if( _isValidFEM1Node( node ) ) matrix[_i][ matrix.rowSizes[_i]++ ] = MatrixEntry< Real >( (matrix_index_type)( node->nodeData.nodeIndex - _sNodesBegin(highDepth) ) , (Real)values[Dim] ); } ,
#endif // NEW_CODE_SPARSE_MATRIX
#else // !NEW_CODE
					[&]( const FEMTreeNode* node ){ if( _isValidFEM1Node( node ) ) matrix[_i][ matrix.rowSizes[_i]++ ] = MatrixEntry< Real >( node->nodeData.nodeIndex - _sNodesBegin(highDepth) , (Real)values[Dim] ); } ,
#endif // NEW_CODE
					neighbors.neighbors()
				);
			}
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	for( int d=0 ; d<Dim ; d++ ) delete upSampleEvaluators[d];
	return matrix;
}

template< unsigned int Dim , class Real >
template< typename T , unsigned int ... PointDs , unsigned int ... FEMSigs >
#ifdef NEW_CODE
#ifdef NEW_THREADS
SparseMatrix< Real , matrix_index_type > FEMTree< Dim , Real >::fullSystemMatrix( UIntPack< FEMSigs ... > , ThreadPool &tp , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth depth , bool nonRefinableOnly , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
SparseMatrix< Real , matrix_index_type > FEMTree< Dim , Real >::fullSystemMatrix( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth depth , bool nonRefinableOnly , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_THREADS
#else // !NEW_CODE
SparseMatrix< Real , int > FEMTree< Dim , Real >::fullSystemMatrix( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , LocalDepth depth , bool nonRefinableOnly , const InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_CODE
{
#ifdef NEW_CODE
	SparseMatrix< Real , matrix_index_type > M;
	std::vector< SparseMatrix< Real , matrix_index_type > >                  systemMatrices( depth+1 );
	std::vector< SparseMatrix< Real , matrix_index_type > >         prolongedSystemMatrices( depth   );
	std::vector< std::vector< SparseMatrix< Real , matrix_index_type > > > upSampleMatrices( depth-1 );
#else // !NEW_CODE
	SparseMatrix< Real , int > M;
	std::vector< SparseMatrix< Real , int > >                  systemMatrices( depth+1 );
	std::vector< SparseMatrix< Real , int > >         prolongedSystemMatrices( depth   );
	std::vector< std::vector< SparseMatrix< Real , int > > > upSampleMatrices( depth-1 );
#endif // NEW_CODE

	for( int d=0 ; d<depth-1 ; d++ ) upSampleMatrices[d].resize( depth );
#ifdef NEW_CODE
	node_index_type size = _sNodesEnd( depth );
#else // !NEW_CODE
	unsigned int size = _sNodesEnd( depth );
#endif // NEW_CODE
	for( int d=0 ; d<=depth ; d++ )
	{
#ifdef NEW_CODE
		SparseMatrix< Real , matrix_index_type >& M = systemMatrices[d];
		M.resize( size );
#ifdef NEW_THREADS
		SparseMatrix< Real , matrix_index_type > _M = systemMatrix< Real >( UIntPack< FEMSigs ... >() , tp , F , d , interpolationInfo ... );
#else // !NEW_THREADS
		SparseMatrix< Real , matrix_index_type > _M = systemMatrix< Real >( UIntPack< FEMSigs ... >() , F , d , interpolationInfo ... );
#endif // NEW_THREADS
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( unsigned int , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
		for( node_index_type i=0 ; i<_M.rows() ; i++ )
#endif // NEW_THREADS
		{
			M.setRowSize( (matrix_index_type)( i + _sNodesBegin(d) ) , _M.rowSize(i) );
#ifdef NEW_CODE_SPARSE_MATRIX
			for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i+_sNodesBegin(d)][j] = MatrixEntry< Real , matrix_index_type >( _M[i][j].N + (matrix_index_type)_sNodesBegin(d) , _M[i][j].Value );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i+_sNodesBegin(d)][j] = MatrixEntry< Real >( _M[i][j].N + (matrix_index_type)_sNodesBegin(d) , _M[i][j].Value );
#endif // NEW_CODE_SPARSE_MATRIX
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
#else // !NEW_CODE
		SparseMatrix< Real , int >& M = systemMatrices[d];
		M.resize( size );
#ifdef NEW_THREADS
		SparseMatrix< Real , int > _M = systemMatrix< Real >( UIntPack< FEMSigs ... >() , tp , F , d , interpolationInfo ... );
#else // !NEW_THREADS
		SparseMatrix< Real , int > _M = systemMatrix< Real >( UIntPack< FEMSigs ... >() , F , d , interpolationInfo ... );
#endif // NEW_THREADS
#pragma omp parallel for
		for( int i=0 ; i<_M.rows() ; i++ )
		{
			M.setRowSize( i + _sNodesBegin(d) , _M.rowSize(i) );
			for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i+_sNodesBegin(d)][j] = MatrixEntry< Real >( _M[i][j].N + _sNodesBegin(d) , _M[i][j].Value );
		}
#endif // NEW_CODE
	}
	for( int d=0 ; d<depth ; d++ )
	{
#ifdef NEW_CODE
		SparseMatrix< Real , matrix_index_type >& M = prolongedSystemMatrices[d];
		M.resize( size );
#ifdef NEW_THREADS
		SparseMatrix< Real , matrix_index_type > _M = prolongedSystemMatrix< Real >( UIntPack< FEMSigs ... >() , tp , F , d+1 , interpolationInfo ... );
#ifdef NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( unsigned int , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
		SparseMatrix< Real , matrix_index_type > _M = prolongedSystemMatrix< Real >( UIntPack< FEMSigs ... >() , F , d+1 , interpolationInfo ... );
#pragma omp parallel for
		for( node_index_type i=0 ; i<_M.rows() ; i++ )
#endif // NEW_THREADS
		{
			M.setRowSize( i + (matrix_index_type)_sNodesBegin(d+1) , _M.rowSize(i) );
#ifdef NEW_CODE_SPARSE_MATRIX
			for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i+_sNodesBegin(d+1)][j] = MatrixEntry< Real , matrix_index_type >( _M[i][j].N + (matrix_index_type)_sNodesBegin(d) , _M[i][j].Value );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i+_sNodesBegin(d+1)][j] = MatrixEntry< Real >( _M[i][j].N + (matrix_index_type)_sNodesBegin(d) , _M[i][j].Value );
#endif // NEW_CODE_SPARSE_MATRIX
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
#else // !NEW_CODE
		SparseMatrix< Real , int >& M = prolongedSystemMatrices[d];
		M.resize( size );
		SparseMatrix< Real , int > _M = prolongedSystemMatrix< Real >( UIntPack< FEMSigs ... >() , F , d+1 , interpolationInfo ... );
#pragma omp parallel for
		for( int i=0 ; i<_M.rows() ; i++ )
		{
			M.setRowSize( i + _sNodesBegin(d+1) , _M.rowSize(i) );
			for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i+_sNodesBegin(d+1)][j] = MatrixEntry< Real >( _M[i][j].N + _sNodesBegin(d) , _M[i][j].Value );
		}
#endif // NEW_CODE
	}
	for( int d=0 ; d<depth-1 ; d++ )
	{
#ifdef NEW_CODE
		SparseMatrix< Real , matrix_index_type >& M = upSampleMatrices[d][d+1];
		M.resize( size );
#ifdef NEW_THREADS
		SparseMatrix< Real , matrix_index_type > _M = downSampleMatrix( UIntPack< FEMSigs ... >() , tp , d+1 ).transpose( _sNodesSize( d+1 ) );
#ifdef NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREAD_NUM
		tp.parallel_for( 0 , _M.rows() , [&]( unsigned int , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
		SparseMatrix< Real , matrix_index_type > _M = downSampleMatrix( UIntPack< FEMSigs ... >() , d+1 ).transpose( _sNodesSize( d+1 ) );
#pragma omp parallel for
		for( node_index_type i=0 ; i<_M.rows() ; i++ )
#endif // NEW_THREADS
		{
			M.setRowSize( i + (matrix_index_type)_sNodesBegin(d+1) , _M.rowSize(i) );
#ifdef NEW_CODE_SPARSE_MATRIX
			for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i+_sNodesBegin(d+1)][j] = MatrixEntry< Real , matrix_index_type >( _M[i][j].N + (matrix_index_type)_sNodesBegin(d) , _M[i][j].Value );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i+_sNodesBegin(d+1)][j] = MatrixEntry< Real >( _M[i][j].N + (matrix_index_type)_sNodesBegin(d) , _M[i][j].Value );
#endif // NEW_CODE_SPARSE_MATRIX
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
#else // !NEW_CODE
		SparseMatrix< Real , int >& M = upSampleMatrices[d][d+1];
		M.resize( size );
		SparseMatrix< Real , int > _M = downSampleMatrix( UIntPack< FEMSigs ... >() , d+1 ).transpose( _sNodesSize( d+1 ) );
#pragma omp parallel for
		for( int i=0 ; i<_M.rows() ; i++ )
		{
			M.setRowSize( i + _sNodesBegin(d+1) , _M.rowSize(i) );
			for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i+_sNodesBegin(d+1)][j] = MatrixEntry< Real >( _M[i][j].N + _sNodesBegin(d) , _M[i][j].Value );
		}
#endif // NEW_CODE
		for( int dd=0 ; dd<d ; dd++ ) upSampleMatrices[dd][d+1] = upSampleMatrices[d][d+1] * upSampleMatrices[dd][d];
	}

	auto Matrix = [&]( int d1 , int d2 )
	{
#ifdef NEW_CODE
		SparseMatrix< Real , matrix_index_type > _M;
#else // !NEW_CODE
		SparseMatrix< Real , int > _M;
#endif // NEW_CODE
		int _d1 = d1<d2 ? d1 : d2 , _d2 = d2<d1 ? d1 : d2;
		if     ( _d1==_d2   ) _M =          systemMatrices[_d1];
		else if( _d2==_d1+1 ) _M = prolongedSystemMatrices[_d2-1];
		else                  _M = prolongedSystemMatrices[_d2-1] * upSampleMatrices[_d1][_d2-1];
		if( d2<d1 ) return _M.transpose( size );
		else        return _M;
	};

	for( int d1=0 ; d1<=depth ; d1++ )
	{
		M += Matrix( d1 , d1 );
		for( int d2=0 ; d2<=depth ; d2++ ) if( d1!=d2 )
		{
#ifdef NEW_CODE
			SparseMatrix< Real , matrix_index_type > _M = Matrix( d1 , d2 );
#else // !NEW_CODE
			SparseMatrix< Real , int > _M = Matrix( d1 , d2 );
#endif // NEW_CODE
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			tp.parallel_for( 0 , _M.rows() , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( 0 , _M.rows() , [&]( unsigned int , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
			for( matrix_index_type i=0 ; i<(matrix_index_type)_M.rows() ; i++ )
#else // !NEW_CODE
			for( int i=0 ; i<_M.rows() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
				if( _M.rowSize(i) )
				{
					size_t oldSize = M.rowSize(i);
					M.resetRowSize( i , oldSize + _M.rowSize(i) );
					for( int j=0 ; j<_M.rowSize(i) ; j++ ) M[i][oldSize+j] = _M[i][j];
				}
			}
#ifdef NEW_THREADS
			);
#endif // NEW_THREADS
		}
	}
	if( nonRefinableOnly )
	{
#ifdef NEW_THREADS
		_setRefinabilityFlags( UIntPack< FEMSigs ... >() , tp );
#else // !NEW_THREADS
		_setRefinabilityFlags( UIntPack< FEMSigs ... >() );
#endif // NEW_THREADS
		_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
		tp.parallel_for( 0 , M.rows() , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREAD_NUM
		tp.parallel_for( 0 , M.rows() , [&]( unsigned int , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( matrix_index_type i=0 ; i<(matrix_index_type)M.rows() ; i++ )
#else // !NEW_CODE
		for( int i=0 ; i<M.rows() ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( ( _isRefinableNode( _sNodes.treeNodes[i] ) && _localDepth( _sNodes.treeNodes[i] )<depth ) || !_isValidFEM1Node( _sNodes.treeNodes[i] ) )
			{
				// Setting this to the (local) identity so it doesn't make the system singular.
				M.resetRowSize( i , 1 );
#ifdef NEW_CODE_SPARSE_MATRIX
				M[i][0] = MatrixEntry< Real , matrix_index_type >( i , (Real)1. );
#else // !NEW_CODE_SPARSE_MATRIX
				M[i][0] = MatrixEntry< Real >( i , (Real)1. );
#endif // NEW_CODE_SPARSE_MATRIX
			}
			else
			{
				int jj=0;
				for( int j=0 ; j<M.rowSize(i) ; j++ ) if( !( _isRefinableNode( _sNodes.treeNodes[ M[i][j].N ] ) && _localDepth( _sNodes.treeNodes[ M[i][j].N ] )<depth ) && _isValidFEM1Node( _sNodes.treeNodes[i] ) ) M[i][jj++] = M[i][j];
				if( jj!=M.rowSize(i) ) M.resetRowSize( i , jj );
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
	return M;
}
template< unsigned int Dim , class Real >
template< class C , unsigned int ... Degrees , unsigned int ... FEMSigs >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_downSample( UIntPack< FEMSigs ... > , ThreadPool &tp , typename BaseFEMIntegrator::template RestrictionProlongation< UIntPack< Degrees ... > >& rp , LocalDepth highDepth , Pointer( C ) constraints ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_downSample( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template RestrictionProlongation< UIntPack< Degrees ... > >& rp , LocalDepth highDepth , Pointer( C ) constraints ) const
#endif // NEW_THREADS
{
	LocalDepth lowDepth = highDepth-1;
	if( lowDepth<0 ) return;

	typedef typename BaseFEMIntegrator::RestrictionProlongation< UIntPack< Degrees ... > > BaseRestrictionProlongation;
	typedef typename FEMTreeNode::template ConstNeighborKey< UIntPack< ( - BSplineSupportSizes< Degrees >::UpSampleStart ) ... > , UIntPack< BSplineSupportSizes< Degrees >::UpSampleEnd ... > > UpSampleKey;
	typedef typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineSupportSizes< Degrees >::UpSampleSize ... > > UpSampleNeighbors;
	typedef UIntPack< BSplineSupportSizes< Degrees >::UpSampleSize ... > UpSampleSizes;

#ifdef NEW_THREADS
	std::vector< UpSampleKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< UpSampleKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( lowDepth ) );

	( ( BaseRestrictionProlongation& )rp ).init( highDepth );
	typename BaseRestrictionProlongation::UpSampleStencil upSampleStencil;
	rp.setStencil( upSampleStencil );

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(lowDepth) , _sNodesEnd(lowDepth) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(lowDepth) , _sNodesEnd(lowDepth) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(lowDepth) ; i<_sNodesEnd(lowDepth) ; i++ )
#else // !NEW_CODE
	for( int i=_sNodesBegin(lowDepth) ; i<_sNodesEnd(lowDepth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
		{
			FEMTreeNode* pNode = _sNodes.treeNodes[i];
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			UpSampleKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			UpSampleKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			UpSampleKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			LocalDepth d ; LocalOffset off;
			_localDepthAndOffset( pNode , d , off );

			neighborKey.getNeighbors( pNode );
			UpSampleNeighbors neighbors;
			neighborKey.getChildNeighbors( 0 , _localToGlobal( d ) , neighbors );

			C& coarseConstraint = constraints[i];

			// Want to make sure test if contained children are interior.
			// This is more conservative because we are test that overlapping children are interior
			bool isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< Degrees ... >() , UIntPack< Degrees ... >() , d , off );
			if( isInterior )
			{
				Pointer( const FEMTreeNode* ) nodes = neighbors.neighbors().data;
				Pointer( double ) stencilValues = upSampleStencil.data;
				for( unsigned int i=0 ; i<WindowSize< UpSampleSizes >::Size ; i++ )
					if( _isValidFEM1Node( nodes[i] ) ) coarseConstraint += (C)( constraints[ nodes[i]->nodeData.nodeIndex ] * (Real)stencilValues[i] );
			}
			else
			{
				ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
				for( int i=0 ; i<WindowSize< UpSampleSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) )
				{
					LocalDepth _d ; LocalOffset _off;
					_localDepthAndOffset( nodes[i] , _d , _off );
					coarseConstraint += (C)( constraints[ nodes[i]->nodeData.nodeIndex ] * (Real)rp.upSampleCoefficient( off , _off ) );
				}
			}
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
#ifdef NEW_THREADS
DenseNodeData< Real , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::supportWeights( UIntPack< FEMSigs ... > , ThreadPool &tp ) const
#else // !NEW_THREADS
DenseNodeData< Real , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::supportWeights( UIntPack< FEMSigs ... > ) const
#endif // NEW_THREADS
{
	typedef typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > > BaseSystem;
	typedef typename BaseFEMIntegrator::template Constraint< UIntPack< FEMSignature< FEMSigs >::Degree ... > , IsotropicUIntPack< Dim , 0 > , 1 > BaseConstraint;
	typedef UIntPack< (  BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree , 0 >::OverlapSize  ) ... >          OverlapSizes;
	typedef UIntPack< ( -BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree , 0 >::OverlapStart ) ... >  LeftFEMCOverlapRadii;
	typedef UIntPack< (  BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree , 0 >::OverlapEnd   ) ... > RightFEMCOverlapRadii;
	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	typename FEMIntegrator::template ScalarConstraint< UIntPack< FEMSigs ... > , ZeroUIntPack< Dim > , IsotropicUIntPack< Dim , FEMTrivialSignature > , ZeroUIntPack< Dim > > F( {1.} );
	DenseNodeData< Real , UIntPack< FEMSigs ... > > weights = initDenseNodeData( UIntPack< FEMSigs ... >() );
	typename BaseConstraint::CCStencil stencil;
#ifdef NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( int d=0 ; d<=_maxDepth ; d++ )
	{
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( d ) );
		F.init( d );
		F.template setStencil< false >( stencil );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
		tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
		tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#else // !NEW_CODE
		for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
			{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
				ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
				ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
				ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS

				FEMTreeNode* node = _sNodes.treeNodes[i];
				typename FEMTreeNode::template ConstNeighbors< OverlapSizes > neighbors;
				LocalOffset off;
				{
					LocalDepth d ; _localDepthAndOffset( node , d , off );
				}
				neighborKey.getNeighbors( LeftFEMCOverlapRadii() , RightFEMCOverlapRadii() , node , neighbors );
				bool isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< FEMSignature< FEMSigs >::Degree ... >() , ZeroUIntPack< Dim >() , d , off );
				double sum=0 , totalSum=0;
				if( isInterior )
				{
					ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
					ConstPointer( Point< double , 1 > ) stencilValues = stencil.data;
					for( int i=0 ; i<WindowSize< UIntPack< BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree , 0 >::OverlapSize ... > >::Size ; i++ )
					{
						double s = stencilValues[i][0];
						totalSum += s;
						if( isValidSpaceNode( nodes[i] ) ) sum += s;
					}
				}
				else
				{
					static const int OverlapStart [] = { BSplineOverlapSizes< FEMSignature< FEMSigs >::Degree , 0 >::OverlapStart ... };
					LocalOffset _off;
					WindowLoop< Dim >::Run
					(
						IsotropicUIntPack< Dim , 0 >() , OverlapSizes() ,
						[&]( int d , int i ){ _off[d] = off[d]+i+OverlapStart[d]; } ,
						[&]( const FEMTreeNode* node )
					{
						double s = F.ccIntegrate( off , _off )[0];
						totalSum += s;
						if( isValidSpaceNode( node ) ) sum += s;
					} ,
						neighbors.neighbors()
						);
				}
				weights[i] = (Real)( sum / totalSum );
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
	return weights;
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
#ifdef NEW_THREADS
DenseNodeData< Real , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::prolongationWeights( UIntPack< FEMSigs ... > , ThreadPool &tp , bool prolongToChildren ) const
#else // !NEW_THREADS
DenseNodeData< Real , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::prolongationWeights( UIntPack< FEMSigs ... > , bool prolongToChildren ) const
#endif // NEW_THREADS
{
	DenseNodeData< Real , UIntPack< FEMSigs ... > > weights = initDenseNodeData( UIntPack< FEMSigs ... >() );

	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	typedef typename BaseFEMIntegrator::RestrictionProlongation< UIntPack< FEMSignature< FEMSigs >::Degree ... > > BaseRestrictionProlongation;
	typedef typename     FEMIntegrator::template RestrictionProlongation< UIntPack< FEMSigs ... > > RestrictionProlongation;

	typename BaseRestrictionProlongation::DownSampleStencils downSampleStencils;
	RestrictionProlongation rp;

	typedef typename FEMTreeNode::template ConstNeighborKey< UIntPack< ( - BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleStart ) ... > , UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleEnd ... > > UpSampleKey;
	typedef typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleSize ... > > UpSampleNeighbors;
	typedef typename FEMTreeNode::template ConstNeighborKey< UIntPack< - BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Start ... > , UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1End ... > > DownSampleKey;
	typedef typename FEMTreeNode::template ConstNeighbors< UIntPack< ( - BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Start + BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1End + 1 ) ... > > DownSampleNeighbors;
	const int      UpSampleStart[] =   { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleStart    ... };
	const int DownSampleStart[2][Dim] = { { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Start ... } , { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1Start ... } };
	const int   DownSampleEnd[2][Dim] = { { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0End   ... } , { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1End   ... } };

#ifdef NEW_THREADS
	std::vector< UpSampleKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< UpSampleKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( _maxDepth-1 ) );

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(_maxDepth) , _sNodesEnd(_maxDepth) , [&]( const ThreadPool::ThreadNum & , size_t i ){ weights[i] = (Real)0.; } );
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(_maxDepth) , _sNodesEnd(_maxDepth) , [&]( unsigned int , size_t i ){ weights[i] = (Real)0.; } );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(_maxDepth) ; i<_sNodesEnd(_maxDepth) ; i++ ) weights[i] = (Real)0.;
#else // !NEW_CODE
	for( int i=_sNodesBegin(_maxDepth) ; i<_sNodesEnd(_maxDepth) ; i++ ) weights[i] = (Real)0.;
#endif // NEW_CODE
#endif // NEW_THREADS

	for( int lowDepth=0 ; lowDepth<_maxDepth ; lowDepth++ )
	{
		( ( BaseRestrictionProlongation& )rp ).init( lowDepth+1 );
		typename BaseRestrictionProlongation::UpSampleStencil upSampleStencil;
		rp.setStencil( upSampleStencil );

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
		tp.parallel_for( _sNodesBegin(lowDepth) , _sNodesEnd(lowDepth) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
		tp.parallel_for( _sNodesBegin(lowDepth) , _sNodesEnd(lowDepth) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=_sNodesBegin(lowDepth) ; i<_sNodesEnd(lowDepth) ; i++ )
#else // !NEW_CODE
		for( int i=_sNodesBegin(lowDepth) ; i<_sNodesEnd(lowDepth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
			{
				FEMTreeNode* pNode = _sNodes.treeNodes[i];
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
				UpSampleKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
				UpSampleKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
				UpSampleKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
				LocalDepth d ; LocalOffset pOff;
				_localDepthAndOffset( pNode , d , pOff );

				neighborKey.getNeighbors( pNode );
				UpSampleNeighbors neighbors;
				neighborKey.getChildNeighbors( 0 , _localToGlobal( d ) , neighbors );

				double partialSum = 0 , totalSum = 0;

				// Want to make sure test if contained children are interior.
				// This is more conservative because we are test that overlapping children are interior
				bool isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< FEMSignature< FEMSigs >::Degree ... >() , UIntPack< FEMSignature< FEMSigs >::Degree ... >() , d , pOff );

				LocalOffset cOff;
				if( isInterior )
				{
					WindowLoop< Dim >::Run
					(
						IsotropicUIntPack< Dim , 0 >() , UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleSize ... >() ,
						[&]( int d , int i ){ cOff[d] = UpSampleStart[d] + pOff[d]*2 + i; } ,
						[&]( const FEMTreeNode* node , double stencilValue )
					{
						if( FEMIntegrator::IsValidFEMNode( UIntPack< FEMSigs ... >() , lowDepth+1 , cOff ) )
						{
							totalSum += stencilValue;
							if( _isValidFEM1Node( node ) ) partialSum += stencilValue;
						}
					} ,
						neighbors.neighbors() , upSampleStencil()
						);
				}
				else
				{
					WindowLoop< Dim >::Run
					(
						IsotropicUIntPack< Dim , 0 >() , UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleSize ... >() ,
						[&]( int d , int i ){ cOff[d] = UpSampleStart[d] + pOff[d]*2 + i; } ,
						[&]( const FEMTreeNode* node )
					{
						if( FEMIntegrator::IsValidFEMNode( UIntPack< FEMSigs ... >() , lowDepth+1 , cOff ) )
						{
							double stencilValue = rp.upSampleCoefficient( pOff , cOff );
							totalSum += stencilValue;
							if( _isValidFEM1Node( node ) ) partialSum += stencilValue;
						}
					} ,
						neighbors.neighbors() 
						);
				}
				weights[i] = (Real)( partialSum / totalSum );
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
	if( prolongToChildren )
	{
#ifdef NEW_THREADS
		std::vector< DownSampleKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
		std::vector< DownSampleKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( _maxDepth-1 ) );

		for( int lowDepth=_maxDepth-1 ; lowDepth>=0 ; lowDepth-- )
		{
			( ( BaseRestrictionProlongation& )rp ).init( lowDepth+1 );
			typename BaseRestrictionProlongation::DownSampleStencils downSampleStencils;
			rp.setStencils( downSampleStencils );

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(lowDepth+1) , _sNodesEnd(lowDepth+1) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(lowDepth+1) , _sNodesEnd(lowDepth+1) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(lowDepth+1) ; i<_sNodesEnd(lowDepth+1) ; i++ )
#else // !NEW_CODE
			for( int i=_sNodesBegin(lowDepth+1) ; i<_sNodesEnd(lowDepth+1) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
				if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
				{
					FEMTreeNode *cNode = _sNodes.treeNodes[i];
					int c = (int)( cNode-cNode->parent->children );

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
					DownSampleKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
					DownSampleKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
					DownSampleKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
					LocalDepth d ; LocalOffset cOff;
					_localDepthAndOffset( cNode , d , cOff );
					DownSampleNeighbors neighbors = neighborKey.getNeighbors( cNode->parent );
					// Want to make sure test if contained children are interior.
					// This is more conservative because we are test that overlapping children are interior
					bool isInterior;
					{
						LocalDepth d ; LocalOffset pOff;
						_localDepthAndOffset( cNode->parent , d , pOff );
						isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< FEMSignature< FEMSigs >::Degree ... >() , UIntPack< FEMSignature< FEMSigs >::Degree ... >() , d , pOff );
					}

					typename BaseRestrictionProlongation::DownSampleStencil& downSampleStencil = downSampleStencils.data[c];
					int start[Dim] , end[Dim];
					for( int d=0 ; d<Dim ; d++ ) start[d] = DownSampleStart[(c>>d)&1][d] - DownSampleStart[0][d] , end[d] = - DownSampleStart[0][d] + DownSampleEnd[(c>>d)&1][d] + 1;

					double partialSum = 0 , totalSum = 0;
					if( isInterior )
					{
						WindowLoop< Dim >::Run
						(
							start , end ,
							[&]( int , int ){ } ,
							[&]( const FEMTreeNode* node , double stencilValue ){ if( _isValidFEM1Node( node ) ) totalSum += stencilValue , partialSum += weights[ node->nodeData.nodeIndex ] * stencilValue; } ,
							neighbors.neighbors() , downSampleStencil()
						);
					}
					else
					{
						WindowLoop< Dim >::Run
						(
							start , end ,
							[&]( int , int ){ } ,
							[&]( const FEMTreeNode* node )
						{
							if( _isValidFEM1Node( node ) )
							{
								LocalDepth d ; LocalOffset pOff;
								_localDepthAndOffset( node , d , pOff );
								double stencilValue = rp.upSampleCoefficient( pOff , cOff );
								totalSum += stencilValue , partialSum += weights[ node->nodeData.nodeIndex ] * stencilValue;
							}
						} ,
							neighbors.neighbors()
							);
					}
					weights[i] = (Real)( partialSum / totalSum );
				}
			}
#ifdef NEW_THREADS
			);
#endif // NEW_THREADS
		}
	}
	return weights;
}

template< unsigned int Dim , class Real >
template< class C , unsigned int ... Degrees , unsigned int ... FEMSigs >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_upSample( UIntPack< FEMSigs ... > , ThreadPool &tp , typename BaseFEMIntegrator::template RestrictionProlongation< UIntPack< Degrees ... > >& rp , LocalDepth highDepth , Pointer( C ) coefficients ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_upSample( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template RestrictionProlongation< UIntPack< Degrees ... > >& rp , LocalDepth highDepth , Pointer( C ) coefficients ) const
#endif // NEW_THREADS
{
	LocalDepth lowDepth = highDepth-1;
	if( lowDepth<0 ) return;
	typedef typename BaseFEMIntegrator::RestrictionProlongation< UIntPack< Degrees ... > > BaseRestrictionProlongation;
	typedef typename FEMTreeNode::template ConstNeighborKey< UIntPack< - BSplineSupportSizes< Degrees >::DownSample0Start ... > , UIntPack< BSplineSupportSizes< Degrees >::DownSample1End ... > > DownSampleKey;
	typedef typename FEMTreeNode::template ConstNeighbors< UIntPack< ( - BSplineSupportSizes< Degrees >::DownSample0Start + BSplineSupportSizes< Degrees >::DownSample1End + 1 ) ... > > DownSampleNeighbors;
	typedef UIntPack< ( - BSplineSupportSizes< Degrees >::DownSample0Start + BSplineSupportSizes< Degrees >::DownSample1End + 1 ) ... > DownSampleSizes;

#ifdef NEW_THREADS
	std::vector< DownSampleKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< DownSampleKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( lowDepth ) );

	( ( BaseRestrictionProlongation& )rp ).init( highDepth );
	typename BaseRestrictionProlongation::DownSampleStencils downSampleStencils;
	rp.setStencils( downSampleStencils );

	const int Start[2][Dim] = { { BSplineSupportSizes< Degrees >::DownSample0Start ... } , { BSplineSupportSizes< Degrees >::DownSample1Start ... } };
	const int   End[2][Dim] = { { BSplineSupportSizes< Degrees >::DownSample0End   ... } , { BSplineSupportSizes< Degrees >::DownSample1End   ... } };

	static const WindowLoopData< UIntPack< ( - BSplineSupportSizes< Degrees >::DownSample0Start + BSplineSupportSizes< Degrees >::DownSample1End + 1 ) ... > > loopData
	( []( int c , int* start , int* end )
	{
		const int Start[2][Dim] = { { BSplineSupportSizes< Degrees >::DownSample0Start ... } , { BSplineSupportSizes< Degrees >::DownSample1Start ... } };
		const int   End[2][Dim] = { { BSplineSupportSizes< Degrees >::DownSample0End   ... } , { BSplineSupportSizes< Degrees >::DownSample1End   ... } };
		for( int d=0 ; d<Dim ; d++ ) start[d] = Start[(c>>d)&1][d] - Start[0][d] , end[d] = - Start[0][d] + End[(c>>d)&1][d] + 1;
	} 
	);
	// For Dirichlet constraints, can't get to all children from parents because boundary nodes are invalid
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(highDepth) , _sNodesEnd(highDepth) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
	tp.parallel_for( _sNodesBegin(highDepth) , _sNodesEnd(highDepth) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
	for( node_index_type i=_sNodesBegin(highDepth) ; i<_sNodesEnd(highDepth) ; i++ )
#else // !NEW_CODE
	for( int i=_sNodesBegin(highDepth) ; i<_sNodesEnd(highDepth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
	{
		if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
		{
			FEMTreeNode *cNode = _sNodes.treeNodes[i];
			int c = (int)( cNode-cNode->parent->children );

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			DownSampleKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			DownSampleKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			DownSampleKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			DownSampleNeighbors neighbors = neighborKey.getNeighbors( cNode->parent );
			// Want to make sure test if contained children are interior.
			// This is more conservative because we are test that overlapping children are interior
			bool isInterior;
			{
				LocalDepth d ; LocalOffset off;
				_localDepthAndOffset( cNode->parent , d , off );
				isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< Degrees ... >() , UIntPack< Degrees ... >() , d , off );
			}

			C& fineCoefficient = coefficients[ cNode->nodeData.nodeIndex ];

			typename BaseRestrictionProlongation::DownSampleStencil& downSampleStencil = downSampleStencils.data[c];
			unsigned int size = loopData.size[c];
			const unsigned int* indices = loopData.indices[c];
			Pointer( const FEMTreeNode* ) nodes = neighbors.neighbors().data;
			Pointer( double ) downSampleValues = downSampleStencil.data;
			if( isInterior )
			{
				for( unsigned int i=0 ; i<size ; i++ )
				{
					unsigned int idx = indices[i];
					if( _isValidFEM1Node( nodes[idx] ) ) fineCoefficient += (C)( coefficients[ nodes[idx]->nodeData.nodeIndex ] * (Real)downSampleValues[idx] );
				}
			}
			else
			{
				LocalDepth d ; LocalOffset off;
				_localDepthAndOffset( cNode , d , off );
				for( unsigned int i=0 ; i<size ; i++ )
				{
					unsigned int idx = indices[i];
					if( _isValidFEM1Node( nodes[idx] ) )
					{
						LocalDepth _d ; LocalOffset _off;
						_localDepthAndOffset( nodes[idx] , _d , _off );
						fineCoefficient += (C)( coefficients[ nodes[idx]->nodeData.nodeIndex ] * (Real)rp.upSampleCoefficient( _off , off ) );
					}
				}
			}
		}
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< unsigned int Dim , class Real >
template< bool XMajor , class C , unsigned int ... FEMSigs >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_RegularGridUpSample( UIntPack< FEMSigs ... > , ThreadPool &tp , LocalDepth highDepth , ConstPointer( C ) lowCoefficients , Pointer( C ) highCoefficients )
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_RegularGridUpSample( UIntPack< FEMSigs ... > , LocalDepth highDepth , ConstPointer( C ) lowCoefficients , Pointer( C ) highCoefficients )
#endif // NEW_THREADS
{
	LocalDepth lowDepth = highDepth - 1;
	if( lowDepth<0 ) return;

	int lowBegin[Dim] , lowEnd[Dim] , highBegin[Dim] , highEnd[Dim];
	FEMIntegrator::BSplineBegin( UIntPack< FEMSigs ... >() ,  lowDepth ,  lowBegin );
	FEMIntegrator::BSplineEnd  ( UIntPack< FEMSigs ... >() ,  lowDepth ,  lowEnd   );
	FEMIntegrator::BSplineBegin( UIntPack< FEMSigs ... >() , highDepth , highBegin );
	FEMIntegrator::BSplineEnd  ( UIntPack< FEMSigs ... >() , highDepth , highEnd   );

#ifdef NEW_THREADS
	_RegularGridUpSample< XMajor >( UIntPack< FEMSigs ... >() , tp , lowBegin , lowEnd , highBegin , highEnd , highDepth , lowCoefficients , highCoefficients );
#else // !NEW_THREADS
	_RegularGridUpSample< XMajor >( UIntPack< FEMSigs ... >() , lowBegin , lowEnd , highBegin , highEnd , highDepth , lowCoefficients , highCoefficients );
#endif // NEW_THREADS
}
template< unsigned int Dim , class Real >
template< bool XMajor , class C , unsigned int ... FEMSigs >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_RegularGridUpSample( UIntPack< FEMSigs ... > , ThreadPool &tp , const int lowBegin[] , const int lowEnd[] , const int highBegin[] , const int highEnd[] , LocalDepth highDepth , ConstPointer( C ) lowCoefficients , Pointer( C ) highCoefficients )
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_RegularGridUpSample( UIntPack< FEMSigs ... > , const int lowBegin[] , const int lowEnd[] , const int highBegin[] , const int highEnd[] , LocalDepth highDepth , ConstPointer( C ) lowCoefficients , Pointer( C ) highCoefficients )
#endif // NEW_THREADS
{
	// Note: In contrast to the standard grid indexing, where x is the major index in (x,y,z,...)
	//       For our representation of the grid, x is the minor index
	LocalDepth lowDepth = highDepth - 1;
	if( lowDepth<0 ) return;

	static const          int LeftDownSampleRadii[] = { -( ( BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Start < BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1Start ) ? BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Start : BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1Start ) ... };
	static const          int DownSampleStart[][ sizeof...(FEMSigs) ] = { { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Start ... } , { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1Start ... } };
	static const unsigned int DownSampleSize [][ sizeof...(FEMSigs) ] = { { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Size  ... } , { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1Size  ... } };
	typedef UIntPack< ( - BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Start + BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1End + 1 ) ... > DownSampleSizes;
	typedef typename     FEMIntegrator::template RestrictionProlongation< UIntPack< FEMSigs ...                         > >     RestrictionProlongation;
	typedef typename BaseFEMIntegrator::template RestrictionProlongation< UIntPack< FEMSignature< FEMSigs >::Degree ... > > BaseRestrictionProlongation;

	RestrictionProlongation rp;
	typename BaseRestrictionProlongation::DownSampleStencils downSampleStencils;
	rp.init( highDepth );
	rp.setStencils( downSampleStencils );

	struct LoopData
	{
		unsigned int size[1<<Dim];
		unsigned int indices[1<<Dim][ WindowSize< UIntPack< ( - BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Start + BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1End + 1 ) ... > >::Size ];
		long long offsets[1<<Dim][ WindowSize< UIntPack< ( - BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample0Start + BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::DownSample1End + 1 ) ... > >::Size ];
		LoopData( const int lowBegin[] , const int lowEnd[] , const int highBegin[] , const int highEnd[] )
		{
			int start[Dim] , end[Dim] , lowDim[Dim] , highDim[Dim];
			for( int d=0 ; d<Dim ; d++ ) lowDim[d] = lowEnd[d] - lowBegin[d] , highDim[d] = highEnd[d] - highBegin[d];

			int lowDimMultiplier[Dim] , highDimMultiplier[Dim];
			if( XMajor )
			{
				lowDimMultiplier[0] = highDimMultiplier[0] = 1;
				for( int d=1 ; d<Dim ; d++ ) lowDimMultiplier[d] = lowDimMultiplier[d-1] * (lowEnd[d-1]-lowBegin[d-1]) , highDimMultiplier[d] = highDimMultiplier[d-1] * (highEnd[d-1]-highBegin[d-1]);
			}
			else
			{
				lowDimMultiplier[Dim-1] = highDimMultiplier[Dim-1] = 1;
				for( int d=Dim-2 ; d>=0 ; d-- ) lowDimMultiplier[d] = lowDimMultiplier[d+1] * (lowEnd[d+1]-lowBegin[d+1]) , highDimMultiplier[d] = highDimMultiplier[d+1] * (highEnd[d+1]-highBegin[d+1]);
			}

			for( int c=0 ; c<(1<<Dim) ; c++ )
			{
				size[c] = 0;
				for( int d=0 ; d<Dim ; d++ ) start[d] = DownSampleStart[(c>>d)&1][d] + LeftDownSampleRadii[d] , end[d] = start[d] + DownSampleSize[(c>>d)&1][d];

				unsigned int idx[Dim];
				long long off[Dim+1];
				off[0] = 0;
				WindowLoop< Dim >::Run
				(
					start , end ,
					[&]( int d , int i ){ idx[d] = i ; off[d+1] = off[d] + ( i - LeftDownSampleRadii[d] - lowBegin[d] ) * lowDimMultiplier[d]; } ,
					[&]( void ){ indices[c][ size[c] ] = GetWindowIndex( DownSampleSizes() , idx ) , offsets[c][ size[c] ] = off[Dim] ; size[c]++; }
				);
			}
		}
	};
	const LoopData loopData( lowBegin , lowEnd , highBegin , highEnd );
	int lowDim[Dim] , highDim[Dim];
	for( int d=0 ; d<Dim ; d++ ) lowDim[d] = lowEnd[d] - lowBegin[d] , highDim[d] = highEnd[d] - highBegin[d];
	int Zero[Dim];
	for( int d=0 ; d<Dim ; d++ ) Zero[d] = 0;
	int lowDimMultiplier[Dim] , highDimMultiplier[Dim];
	if( XMajor )
	{
		lowDimMultiplier[0] = highDimMultiplier[0] = 1;
		for( int d=1 ; d<Dim ; d++ ) lowDimMultiplier[d] = lowDimMultiplier[d-1] * (lowEnd[d-1]-lowBegin[d-1]) , highDimMultiplier[d] = highDimMultiplier[d-1] * (highEnd[d-1]-highBegin[d-1]);
	}
	else
	{
		lowDimMultiplier[Dim-1] = highDimMultiplier[Dim-1] = 1;
		for( int d=Dim-2 ; d>=0 ; d-- ) lowDimMultiplier[d] = lowDimMultiplier[d+1] * (lowEnd[d+1]-lowBegin[d+1]) , highDimMultiplier[d] = highDimMultiplier[d+1] * (highEnd[d+1]-highBegin[d+1]);
	}


	struct UpdateData
	{
		typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;
		LocalOffset pOff , cOff;
		int c;
		long long lowIndex[Dim+1] , highIndex[Dim+1];
		bool isInterior[Dim+1];
		int start[Dim] , end[Dim];
		void init( int lowDepth , const int lowBegin[] , const int lowEnd[] , const int highBegin[] , const int highEnd[] )
		{
			c = 0;
			lowIndex[0] = highIndex[0] = 0;
			isInterior[0] = true;
			this->lowBegin = lowBegin , this->lowEnd = lowEnd , this->highBegin = highBegin , this->highEnd = highEnd;
			if( XMajor )
			{
				_lowDim[0] = _highDim[0] = 1;
				for( int d=1 ; d<Dim ; d++ ) _lowDim[d] = _lowDim[d-1] * (lowEnd[d-1]-lowBegin[d-1]) , _highDim[d] = _highDim[d-1] * (highEnd[d-1]-highBegin[d-1]);
			}
			else
			{
				_lowDim[Dim-1] = _highDim[Dim-1] = 1;
				for( int d=Dim-2 ; d>=0 ; d-- ) _lowDim[d] = _lowDim[d+1] * (lowEnd[d+1]-lowBegin[d+1]) , _highDim[d] = _highDim[d+1] * (highEnd[d+1]-highBegin[d+1]);
			}
			BaseFEMIntegrator::InteriorOverlappedSpan( Degrees() , Degrees() , lowDepth , _begin , _end );
		}
		void set( int d , int i )
		{
			int ii = i + highBegin[d];
			cOff[d] = ii;
			pOff[d] = (ii>>1);
			c = ( c & ( ~(1<<d) ) ) | (ii&1)<<d;
			lowIndex[d+1] = lowIndex[d] + pOff[d] * _lowDim[d];
			highIndex[d+1] = highIndex[d] + i * _highDim[d];
			start[d] = DownSampleStart[(c>>d)&1][d] + LeftDownSampleRadii[d] , end[d] = start[d] + DownSampleSize[(c>>d)&1][d];
			isInterior[d+1] = isInterior[d] && ( pOff[d] + start[d] - LeftDownSampleRadii[d] )>=lowBegin[d] && ( pOff[d] + end[d] - LeftDownSampleRadii[d] )<lowEnd[d] && pOff[d]>=_begin[d] && pOff[d]<_end[d];
		}
	protected:
		const int *lowBegin , *lowEnd , *highBegin , *highEnd;
		int _lowDim[Dim] , _highDim[Dim] , _begin[Dim] , _end[Dim];
	};
#ifdef NEW_THREADS
	std::vector< UpdateData > updateData( tp.threadNum() );
#else // !NEW_THREADS
	std::vector< UpdateData > updateData( omp_get_max_threads() );
#endif // NEW_THREADS
	for( int i=0 ; i<updateData.size() ; i++ ) updateData[i].init( lowDepth , lowBegin , lowEnd , highBegin , highEnd );
	WindowLoop< Dim >::RunParallel
	(
#ifdef NEW_THREADS
		tp ,
#endif // NEW_THREADS
		Zero , highDim ,
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
		[&]( const ThreadPool::ThreadNum &t , int d , size_t i ){ updateData[t()].set( d , (int)i ); } ,
		[&]( const ThreadPool::ThreadNum &t )
#else // !NEW_THREAD_NUM
		[&]( unsigned int t , int d , size_t i ){ updateData[t].set( d , (int)i ); } ,
		[&]( unsigned int t )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
		[&]( int t , int d , int i ){ updateData[t].set( d , i ); } ,
		[&]( int t )
#endif // NEW_THREADS
		{
#ifdef NEW_THREAD_NUM
			const UpdateData& data = updateData[t()];
#else // !NEW_THREAD_NUM
			const UpdateData& data = updateData[t];
#endif // NEW_THREAD_NUM
			const long long highIdx = data.highIndex[Dim] , lowIndex = data.lowIndex[Dim];
			const int c = data.c;
			const bool isInterior = data.isInterior[Dim];

			C highCoefficient = {};

			if( isInterior )
			{
				typename BaseRestrictionProlongation::DownSampleStencil& downSampleStencil = downSampleStencils.data[c];
				const unsigned int size = loopData.size[c];
				const unsigned int* idx = loopData.indices[c];
				const long long* off = loopData.offsets[c];
				ConstPointer( double ) stencilValues = downSampleStencil.data;
				ConstPointer( C ) _lowCoefficients = lowCoefficients + lowIndex;
				for( unsigned int i=0 ; i<size ; i++ ) highCoefficient += (C)( _lowCoefficients[ off[i] ] * (Real)stencilValues[ idx[i] ] );
			}
			else
			{
				const LocalOffset& pOff = data.pOff;
				const LocalOffset& cOff = data.cOff;
				const int* start = data.start;
				const int* end = data.end;
				long long lowIdx[ Dim+1 ] ; lowIdx[0] = 0;
				bool isValid[Dim+1] ; isValid[0] = true;
				int _pOff[Dim];

				WindowLoop< Dim >::Run
				(
					start , end ,
					[&]( int d , int i )
					{
						_pOff[d] = pOff[d] + i - LeftDownSampleRadii[d];
						lowIdx[d+1] = lowIdx[d] + lowDimMultiplier[d] * ( _pOff[d] - lowBegin[d] );
						isValid[d+1] = isValid[d] && ( _pOff[d]>=lowBegin[d] && _pOff[d]<lowEnd[d] );
					} ,
					[&]( void ){ if( isValid[Dim] ) highCoefficient += (C)( lowCoefficients[ lowIdx[Dim] ] * (Real)rp.upSampleCoefficient( _pOff , cOff ) ); }
				);
			}
			highCoefficients[ highIdx ] += highCoefficient;
		}
	);
}



template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , typename TDotT , unsigned int ... PointDs >
#ifdef NEW_THREADS
DenseNodeData< T , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::solveSystem( UIntPack< FEMSigs ... > , ThreadPool &tp , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const DenseNodeData< T , UIntPack< FEMSigs ... > >& constraints , TDotT Dot , LocalDepth maxSolveDepth , const typename FEMTree< Dim , Real >::SolverInfo& solverInfo , InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
DenseNodeData< T , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::solveSystem( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const DenseNodeData< T , UIntPack< FEMSigs ... > >& constraints , TDotT Dot , LocalDepth maxSolveDepth , const typename FEMTree< Dim , Real >::SolverInfo& solverInfo , InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_THREADS
{
	DenseNodeData< T , UIntPack< FEMSigs ... > > solution;
#ifdef NEW_THREADS
	solveSystem( UIntPack< FEMSigs ... >() , tp , F , constraints , solution , Dot , maxSolveDepth , solverInfo , interpolationInfo... );
#else // !NEW_THREADS
	solveSystem( UIntPack< FEMSigs ... >() , F , constraints , solution , Dot , maxSolveDepth , solverInfo , interpolationInfo... );
#endif // NEW_THREADS
	return solution;
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , typename T , typename TDotT , unsigned int ... PointDs >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::solveSystem( UIntPack< FEMSigs ... > , ThreadPool &tp , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const DenseNodeData< T , UIntPack< FEMSigs ... > >& constraints , DenseNodeData< T , UIntPack< FEMSigs ... > >& solution , TDotT Dot , LocalDepth maxSolveDepth , const typename FEMTree< Dim , Real >::SolverInfo& solverInfo , InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::solveSystem( UIntPack< FEMSigs ... > , typename BaseFEMIntegrator::template System< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& F , const DenseNodeData< T , UIntPack< FEMSigs ... > >& constraints , DenseNodeData< T , UIntPack< FEMSigs ... > >& solution , TDotT Dot , LocalDepth maxSolveDepth , const typename FEMTree< Dim , Real >::SolverInfo& solverInfo , InterpolationInfo< T , PointDs >* ... interpolationInfo ) const
#endif // NEW_THREADS
{
	int baseDepth = solverInfo.baseDepth;
	if( baseDepth>getFullDepth( UIntPack< FEMSignature< FEMSigs >::Degree ... >() ) ) ERROR_OUT( "Base depth cannot excceed full depth: " , baseDepth , " <= " , getFullDepth( UIntPack< FEMSignature< FEMSigs >::Degree ... >() ) );

	static_assert( Dim==sizeof ... ( FEMSigs ) , "[ERROR] FEMTree:solveSystem: Dimensions and number of signatures don't match" );
	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	PointEvaluator< UIntPack< FEMSigs ... > , UIntPack< FEMSignature< FEMSigs >::Degree ... > > bsData( sizeof...(PointDs)==0 ? 0 : maxSolveDepth );

	maxSolveDepth = std::min< LocalDepth >( maxSolveDepth , _maxDepth );

	bool clearSolution = solution.size()!=_sNodesEnd( _maxDepth );
	if( clearSolution ) solution = initDenseNodeData< T >( UIntPack< FEMSigs ... >() ) , clearSolution = true;
	bool simpleSolve = clearSolution && solverInfo.vCycles==1 && solverInfo.cascadic;

	// The initial estimate of the solution (may be empty or may come in with an initial guess)
	Pointer( T ) _solution = solution();
	// The constraints
	ConstPointer( T ) _constraints = constraints();

	// _residualConstraints:
	// -- stores the difference between the initial constraints and the constraints met by the current solution at all _other_ levels
	// **** This could implemented in one of two ways:
	// **** (1) Repeatedly computing the difference using the entire solution
	// **** (2) Iteratively updating using the change in the solution
	// **** We have opted for #1 to avoid having to compute/store the change in the solution after each solve
	Pointer( T ) _residualConstraints = AllocPointer< T >( _sNodesEnd( _maxDepth-1 ) );
	// The constraints met during the restriction phase
	Pointer( T ) _restrictedConstraints = NullPointer( T );
	// The solution obtained during the prolongation phase
	Pointer( T ) _prolongedSolution = AllocPointer< T >( _sNodesEnd( _maxDepth-1 ) );

	memset( _prolongedSolution , 0 , sizeof(T) * _sNodesEnd( _maxDepth-1 ) );
	if( !( clearSolution && solverInfo.vCycles==1 && solverInfo.cascadic ) )
	{
		_restrictedConstraints = AllocPointer< T >( _sNodesEnd( _maxDepth-1 ) );
		memset( _restrictedConstraints , 0 , sizeof(T) * _sNodesEnd( _maxDepth-1 ) );
	}

	Pointer( double ) _bNorm2 = NullPointer( double );
	if( solverInfo.showGlobalResidual!=SHOW_GLOBAL_RESIDUAL_NONE )
	{
		_bNorm2 = AllocPointer< double >( _maxDepth+1 );
		memset( _bNorm2 , 0 , sizeof(double) * ( _maxDepth+1 ) );
#ifdef NEW_CODE
		for( LocalDepth d=baseDepth ; d<=maxSolveDepth ; d++ ) for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ ) _bNorm2[d] += Dot( _constraints[i] , _constraints[i] );
#else // !NEW_CODE
		for( LocalDepth d=baseDepth ; d<=maxSolveDepth ; d++ ) for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ ) _bNorm2[d] += Dot( _constraints[i] , _constraints[i] );
#endif // NEW_CODE
	}

	auto UpdateProlongation = [&] ( int depth )
	{
		if( depth<_maxDepth && _prolongedSolution )
		{
			memset( _prolongedSolution + _sNodesBegin( depth ) , 0 , sizeof( T ) * _sNodesSize( depth ) );
			// Up-sample the prolonged solution @(depth-1) into the prolonged solution @(depth)
			F.init( depth );
#ifdef NEW_THREADS
			if( depth>baseDepth ) _upSample( UIntPack< FEMSigs ... >() , tp , F.restrictionProlongation() , depth , _prolongedSolution );
			// Add in the solution @(depth) to the prolonged solution
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd(depth) , [&]( const ThreadPool::ThreadNum & , size_t i ){ _prolongedSolution[i] += solution[i]; } );
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd(depth) , [&]( unsigned int , size_t i ){ _prolongedSolution[i] += solution[i]; } );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			if( depth>baseDepth ) _upSample( UIntPack< FEMSigs ... >() , F.restrictionProlongation() , depth , _prolongedSolution );
			// Add in the solution @(depth) to the prolonged solution
#pragma omp parallel for
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(depth) ; i<_sNodesEnd(depth) ; i++ ) _prolongedSolution[i] += solution[i];
#else // !NEW_CODE
			for( int i=_sNodesBegin(depth) ; i<_sNodesEnd(depth) ; i++ ) _prolongedSolution[i] += solution[i];
#endif // NEW_CODE
#endif // NEW_THREADS
		}
	};
	auto UpdateRestriction = [&]( int depth , InterpolationInfo< T , PointDs >* ... interpolationInfo )
	{
		if( depth>baseDepth && _restrictedConstraints )
		{
			memset( _restrictedConstraints + _sNodesBegin( depth-1 ) , 0 , sizeof( T ) * _sNodesSize( depth-1 ) );
			// Update the restricted constraints @(depth-1) based on the solution @(depth)
			F.init( depth );
#ifdef NEW_THREADS
			_updateRestrictedIntegralConstraints( UIntPack< FEMSigs ... >() , tp , F , depth , ( ConstPointer(T) )_solution , _restrictedConstraints );
			_updateRestrictedInterpolationConstraints( tp , bsData , depth , ( ConstPointer( T ) )_solution , _restrictedConstraints , interpolationInfo... );
			// Down-sample the restricted constraints @(depth) into the restricted constraints @(depth-1)
			if( depth<_maxDepth ) _downSample( UIntPack< FEMSigs ... >() , tp , F.restrictionProlongation() , depth , _restrictedConstraints );
#else // !NEW_THREADS
			_updateRestrictedIntegralConstraints( UIntPack< FEMSigs ... >() , F , depth , ( ConstPointer(T) )_solution , _restrictedConstraints );
			_updateRestrictedInterpolationConstraints( bsData , depth , ( ConstPointer( T ) )_solution , _restrictedConstraints , interpolationInfo... );
			// Down-sample the restricted constraints @(depth) into the restricted constraints @(depth-1)
			if( depth<_maxDepth ) _downSample( UIntPack< FEMSigs ... >() , F.restrictionProlongation() , depth , _restrictedConstraints );
#endif // NEW_THREADS
		}
	};
	auto SetResidualConstraints = [&]( int depth , InterpolationInfo< T , PointDs >* ... interpolationInfo )
	{
		// Copy the constraints
		if( depth<_maxDepth ) memcpy( _residualConstraints + _sNodesBegin(depth) , _constraints + _sNodesBegin(depth) , sizeof( T ) * _sNodesSize(depth) );

		// Update the constraints @(depth) using the prolonged solution @(depth-1)
#ifdef NEW_THREADS
		if( depth>baseDepth && _prolongedSolution ) _setPointValuesFromProlongedSolution( tp , depth , bsData , ( ConstPointer( T ) )_prolongedSolution , interpolationInfo... );
#else // !NEW_THREADS
		if( depth>baseDepth && _prolongedSolution ) _setPointValuesFromProlongedSolution( depth , bsData , ( ConstPointer( T ) )_prolongedSolution , interpolationInfo... );
#endif // NEW_THREADS
		// Update the constraints @(depth) using the restriced residual @(depth)
		if( depth<_maxDepth && _restrictedConstraints )
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd(depth) , [&]( const ThreadPool::ThreadNum & , size_t i ){ _residualConstraints[i] -= _restrictedConstraints[i]; } );
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(depth) , _sNodesEnd(depth) , [&]( unsigned int , size_t i ){ _residualConstraints[i] -= _restrictedConstraints[i]; } );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(depth) ; i<_sNodesEnd(depth) ; i++ ) _residualConstraints[i] -= _restrictedConstraints[i];
#else // !NEW_CODE
			for( int i=_sNodesBegin(depth) ; i<_sNodesEnd(depth) ; i++ ) _residualConstraints[i] -= _restrictedConstraints[i];
#endif // NEW_CODE
#endif // NEW_THREADS
	};
	auto OutputSolverStats = [&] ( int cycle , int depth , const _SolverStats& sStats , bool showResidual , int actualIters )
	{
		if( solverInfo.verbose )
		{
#ifdef NEW_CODE
			node_index_type femNodes = (node_index_type)validFEMNodes( UIntPack< FEMSigs ... >() , depth );
#else // !NEW_CODE
			int femNodes = (int)validFEMNodes( UIntPack< FEMSigs ... >() , depth );
#endif // NEW_CODE
			if( maxSolveDepth<10 )
				if( solverInfo.vCycles<10 ) printf( "Cycle[%d] Depth[%d/%d]:\t" , cycle , depth , maxSolveDepth );
				else                        printf( "Cycle[%2d] Depth[%d/%d]:\t" , cycle , depth , maxSolveDepth );
			else 
				if( solverInfo.vCycles<10 ) printf( "Cycle[%d] Depth[%2d/%d]:\t" , cycle , depth , maxSolveDepth );
				else                        printf( "Cycle[%2d] Depth[%2d/%d]:\t" , cycle , depth , maxSolveDepth );
#ifdef NEW_CODE
			printf( "Updated constraints / Got system / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\tNodes: %llu\n" , sStats.constraintUpdateTime , sStats.systemTime , sStats.solveTime , _LocalMemoryUsage , (unsigned long long)femNodes );
#else // !NEW_CODE
			printf( "Updated constraints / Got system / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\tNodes: %d\n" , sStats.constraintUpdateTime , sStats.systemTime , sStats.solveTime , _LocalMemoryUsage , femNodes );
#endif // NEW_CODE
		}
		if( solverInfo.showResidual && showResidual )
		{
			for( int d=baseDepth ; d<depth ; d++ ) printf( "  " );
			printf( "%s: %.4e -> %.4e -> %.4e (%.1e) [%d]\n" , depth<=solverInfo.cgDepth ? "CG" : "GS" , sqrt( sStats.bNorm2 ) , sqrt( sStats.inRNorm2 ) , sqrt( sStats.outRNorm2 ) , sqrt( sStats.outRNorm2  / sStats.inRNorm2 ) , actualIters );
		}
	};

	// Set the cumulative solution
	if( !clearSolution ) for( LocalDepth d=baseDepth ; d<maxSolveDepth ; d++ ) UpdateProlongation( d );

	_SolverStats sStats;
	bool showResidual;
	int actualIters;
	double t;

#ifdef NEW_CODE
	struct TrivialSORWeights{ Real operator[] ( node_index_type idx ) const { return (Real)1; } };
#else // !NEW_CODE
	struct TrivialSORWeights{ Real operator[] ( int idx ) const { return (Real)1; } };
#endif // NEW_CODE
	struct SORWeights
	{
		DenseNodeData< Real,  UIntPack< FEMSigs ... > > supportWeights , prolongationSupportWeights;
		std::function< Real (Real,Real) > sorFunction;
#ifdef NEW_CODE
		Real operator[] ( node_index_type idx ) const
#else // !NEW_CODE
		Real operator[] ( int idx ) const
#endif // NEW_CODE
		{
			if     ( supportWeights() && prolongationSupportWeights() ) return sorFunction( supportWeights[idx] , prolongationSupportWeights[idx] );
			else if( supportWeights()                                 ) return sorFunction( supportWeights[idx] , 1                               );
			else if(                     prolongationSupportWeights() ) return sorFunction( 1                   , prolongationSupportWeights[idx] );
			else                                                        return sorFunction( 1                   , 1                               );
		}
	};
	SORWeights sorWeights;
#ifdef NEW_THREADS
	if( solverInfo.useSupportWeights ) sorWeights.supportWeights = supportWeights( UIntPack< FEMSigs ... >() , tp );
	if( solverInfo.useProlongationSupportWeights ) sorWeights.prolongationSupportWeights = prolongationWeights( UIntPack< FEMSigs ... >() , tp , false );
#else // !NEW_THREADS
	if( solverInfo.useSupportWeights ) sorWeights.supportWeights = supportWeights( UIntPack< FEMSigs ... >() );
	if( solverInfo.useProlongationSupportWeights ) sorWeights.prolongationSupportWeights = prolongationWeights( UIntPack< FEMSigs ... >() , false );
#endif // NEW_THREADS

	auto SolveRestriction = [&]( int v , int depth , InterpolationInfo< T , PointDs >* ... interpolationInfo )
	{
		sorWeights.sorFunction = solverInfo.sorRestrictionFunction;
		// The restriction phase
		if( solverInfo.cascadic )
		{
			showResidual = false;
			if( !clearSolution || v>0 ) for( LocalDepth d=depth ; d>=baseDepth ; d-- ) { F.init( d ) ; UpdateRestriction( d , interpolationInfo ... ); }
		}
		else
		{
			bool coarseToFine = false;
			for( LocalDepth d=depth ; d>=baseDepth ; d-- )
			{
				sStats.constraintUpdateTime = 0;
				showResidual = ( d!=baseDepth );
				int iters = solverInfo.iters( v , true , d );
				t = Time();
				F.init( d );
				SetResidualConstraints( d , interpolationInfo... );
				sStats.constraintUpdateTime += Time()-t;
				// In the restriction phase we do not solve at the coarsest resolution since we will do so in the prolongation phase
#ifdef NEW_THREADS
				if( d==baseDepth ) _solveRegularMG( UIntPack< FEMSigs ... >() , tp , F , bsData , d , _solution , d==_maxDepth ? _constraints : _residualConstraints , Dot , solverInfo.baseVCycles , iters , sStats , solverInfo.showResidual , solverInfo.cgAccuracy , interpolationInfo... );
#else // !NEW_THREADS
				if( d==baseDepth ) _solveRegularMG( UIntPack< FEMSigs ... >() , F , bsData , d , _solution , d==_maxDepth ? _constraints : _residualConstraints , Dot , solverInfo.baseVCycles , iters , sStats , solverInfo.showResidual , solverInfo.cgAccuracy , interpolationInfo... );
#endif // NEW_THREADS
				else
				{
#ifdef NEW_THREADS
					if( d>solverInfo.cgDepth ) actualIters = _solveSystemGS( UIntPack< FEMSigs ... >() , tp , Dim!=1 , F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , iters , coarseToFine , solverInfo.sliceBlockSize , sorWeights , sStats , solverInfo.showResidual ,                         interpolationInfo... );
					else                       actualIters = _solveSystemCG( UIntPack< FEMSigs ... >() , tp ,          F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , iters , coarseToFine ,                                          sStats , solverInfo.showResidual , solverInfo.cgAccuracy , interpolationInfo... );
#else // !NEW_THREADS
					if( d>solverInfo.cgDepth ) actualIters = _solveSystemGS( UIntPack< FEMSigs ... >() , Dim!=1 , F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , iters , coarseToFine , solverInfo.sliceBlockSize , sorWeights , sStats , solverInfo.showResidual ,                         interpolationInfo... );
					else                       actualIters = _solveSystemCG( UIntPack< FEMSigs ... >() ,          F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , iters , coarseToFine ,                                          sStats , solverInfo.showResidual , solverInfo.cgAccuracy , interpolationInfo... );
#endif // NEW_THREADS
				}
				t = Time();
				UpdateRestriction( d , interpolationInfo... );
				sStats.constraintUpdateTime += Time()-t;
				OutputSolverStats( v , d , sStats , showResidual , actualIters );
			}
		}
	};
	auto SolveProlongation = [&]( int v , int depth , InterpolationInfo< T , PointDs >* ... interpolationInfo )
	{
		sorWeights.sorFunction = solverInfo.sorProlongationFunction;
		showResidual = true;
		bool coarseToFine = true;
		for( LocalDepth d=baseDepth ; d<=depth ; d++ )
		{
			sStats.constraintUpdateTime = 0;
			int iters = solverInfo.iters( v , false , d );
			t = Time();
			F.init( d );
			SetResidualConstraints( d , interpolationInfo... );
			sStats.constraintUpdateTime += Time()-t;
#ifdef NEW_THREADS
			if( d==baseDepth ) _solveRegularMG( UIntPack< FEMSigs ... >() , tp , F , bsData , d , _solution , d==_maxDepth ? _constraints : _residualConstraints , Dot , solverInfo.baseVCycles , iters , sStats , solverInfo.showResidual , solverInfo.cgAccuracy , interpolationInfo... );
#else // !NEW_THREADS
			if( d==baseDepth ) _solveRegularMG( UIntPack< FEMSigs ... >() , F , bsData , d , _solution , d==_maxDepth ? _constraints : _residualConstraints , Dot , solverInfo.baseVCycles , iters , sStats , solverInfo.showResidual , solverInfo.cgAccuracy , interpolationInfo... );
#endif // NEW_THREADS
			else
			{
#ifdef NEW_THREADS
				if( d>solverInfo.cgDepth ) actualIters = _solveSystemGS( UIntPack< FEMSigs ... >() , tp , Dim!=1 , F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , iters , coarseToFine , solverInfo.sliceBlockSize , sorWeights , sStats , solverInfo.showResidual , interpolationInfo... );
				else                       actualIters = _solveSystemCG( UIntPack< FEMSigs ... >() , tp ,          F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , iters , coarseToFine , sStats , solverInfo.showResidual , solverInfo.cgAccuracy , interpolationInfo... );
#else // !NEW_THREADS
				if( d>solverInfo.cgDepth ) actualIters = _solveSystemGS( UIntPack< FEMSigs ... >() , Dim!=1 , F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , iters , coarseToFine , solverInfo.sliceBlockSize , sorWeights , sStats , solverInfo.showResidual , interpolationInfo... );
				else                       actualIters = _solveSystemCG( UIntPack< FEMSigs ... >() ,          F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , iters , coarseToFine , sStats , solverInfo.showResidual , solverInfo.cgAccuracy , interpolationInfo... );
#endif // NEW_THREADS
			}
			t = Time();
			UpdateProlongation( d );
			sStats.constraintUpdateTime += Time()-t;
			OutputSolverStats( v , d , sStats , showResidual , actualIters );
		}
	};

	for( int v=0 ; v<solverInfo.vCycles ; v++ )
	{
		if( solverInfo.wCycle )
		{
			for( int d=maxSolveDepth ; d>baseDepth ; d-- )
			{
				SolveRestriction ( v , d   , interpolationInfo ... );
				SolveProlongation( v , d-1 , interpolationInfo ... );
			}
			for( int d=baseDepth+1 ; d<=maxSolveDepth ; d++ )
			{
				SolveRestriction ( v , d-1 , interpolationInfo ... );
				SolveProlongation( v , d   , interpolationInfo ... );
			}
		}
		else
		{
			SolveRestriction ( v , maxSolveDepth , interpolationInfo ... );
			SolveProlongation( v , maxSolveDepth , interpolationInfo ... );
		}
		if( solverInfo.showGlobalResidual==SHOW_GLOBAL_RESIDUAL_ALL || ( solverInfo.showGlobalResidual==SHOW_GLOBAL_RESIDUAL_LAST && v==solverInfo.vCycles-1 ) )
		{
			bool coarseToFine = false;
			std::vector< double > rNorms( maxSolveDepth+1 );
			for( LocalDepth d=maxSolveDepth ; d>=baseDepth ; d-- )
			{
				F.init( d );
				SetResidualConstraints( d , interpolationInfo... );
#ifdef NEW_THREADS
				_solveSystemGS( UIntPack< FEMSigs ... >() , tp , Dim!=1 , F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , 0 , coarseToFine , solverInfo.sliceBlockSize , TrivialSORWeights() , sStats , true , interpolationInfo... );
#else // !NEW_THREADS
				_solveSystemGS( UIntPack< FEMSigs ... >() , Dim!=1 , F , bsData , d , _solution , ( ConstPointer( T ) )_prolongedSolution , d==_maxDepth ? _constraints : _residualConstraints , Dot , 0 , coarseToFine , solverInfo.sliceBlockSize , TrivialSORWeights() , sStats , true , interpolationInfo... );
#endif // NEW_THREADS
				UpdateRestriction( d , interpolationInfo... );
				rNorms[d] = sqrt( sStats.outRNorm2 / _bNorm2[d] );
			}
			printf( "%3d" , v+1 );
			for( int d=baseDepth ; d<=maxSolveDepth ; d++ ) printf( "\t%.4e" , rNorms[d] );
			printf( "\n" );
		}
	}
	MemoryUsage();

	FreePointer( _residualConstraints );
	FreePointer( _restrictedConstraints );
	FreePointer( _prolongedSolution );
	FreePointer( _bNorm2 );
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
DenseNodeData< Real , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::initDenseNodeData( UIntPack< FEMSigs ... > ) const
{
	DenseNodeData< Real , UIntPack< FEMSigs ... > > constraints( _sNodes.size() );
	memset( constraints() , 0 , sizeof(Real)*_sNodes.size() );
	return constraints;
}
template< unsigned int Dim , class Real >
template< class Data , unsigned int ... FEMSigs >
DenseNodeData< Data , UIntPack< FEMSigs ... > > FEMTree< Dim , Real >::initDenseNodeData( UIntPack< FEMSigs ... > ) const
{
	DenseNodeData< Data , UIntPack< FEMSigs ... > > constraints( _sNodes.size() );
	memset( constraints() , 0 , sizeof(Data)*_sNodes.size() );
	return constraints;
}

template< unsigned int Dim , class Real > template< class SReal , class Data , unsigned int _Dim > Data FEMTree< Dim , Real >::_StencilDot( Point< SReal , _Dim > p1 , Point< Data , _Dim > p2 ){ Data dot={} ; for( int d=0 ; d<_Dim ; d++ ) dot += p2[d] * (Real)p1[d] ; return dot; }
template< unsigned int Dim , class Real > template< class SReal , class Data                     > Data FEMTree< Dim , Real >::_StencilDot( Point< SReal , 1 >    p1 , Point< Data , 1 >    p2 ){ return p2[0] * (Real)p1[0]; }
template< unsigned int Dim , class Real > template< class SReal , class Data                     > Data FEMTree< Dim , Real >::_StencilDot( SReal                 p1 , Point< Data , 1 >    p2 ){ return p2[0] * (Real)p1; }
template< unsigned int Dim , class Real > template< class SReal , class Data                     > Data FEMTree< Dim , Real >::_StencilDot( Point< SReal , 1 >    p1 , Data                 p2 ){ return p2 * (Real)p1[0]; }
template< unsigned int Dim , class Real > template< class SReal , class Data                     > Data FEMTree< Dim , Real >::_StencilDot( SReal                 p1 , Data                 p2 ){ return p2*(Real)p1; }
template< unsigned int Dim , class Real > template< class Real1 , unsigned int _Dim > bool FEMTree< Dim , Real >::_IsZero( Point< Real1 , _Dim > p ){ for( int d=0 ; d<_Dim ; d++ ) if( !_IsZero( p[d] ) ) return false ; return true; }
template< unsigned int Dim , class Real > template< class Real1 > bool FEMTree< Dim , Real >::_IsZero( Real1 p ){ return p==0; }

template< unsigned int Dim , class Real >
template< typename T , unsigned int ... FEMSigs , unsigned int ... CSigs , unsigned int ... FEMDegrees , unsigned int ... CDegrees , unsigned int CDim , class Coefficients >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::_addFEMConstraints( UIntPack< FEMSigs ... > , UIntPack< CSigs ... > , ThreadPool &tp , typename BaseFEMIntegrator::template Constraint< UIntPack< FEMDegrees ... > , UIntPack< CDegrees ... > , CDim >& F , const Coefficients& coefficients , Pointer( T ) constraints , LocalDepth maxDepth ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::_addFEMConstraints( UIntPack< FEMSigs ... > , UIntPack< CSigs ... > , typename BaseFEMIntegrator::template Constraint< UIntPack< FEMDegrees ... > , UIntPack< CDegrees ... > , CDim >& F , const Coefficients& coefficients , Pointer( T ) constraints , LocalDepth maxDepth ) const
#endif // NEW_THREADS
{
	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	_setFEM2ValidityFlags( UIntPack<   CSigs ... >() );
	typedef typename BaseFEMIntegrator::template Constraint< UIntPack< FEMDegrees ... > , UIntPack< CDegrees ... > , CDim > BaseConstraint;
	typedef typename Coefficients::data_type D;
	typedef UIntPack< (  BSplineOverlapSizes< CDegrees , FEMDegrees >::OverlapSize  ) ... >          OverlapSizes;
	typedef UIntPack< ( -BSplineOverlapSizes< CDegrees , FEMDegrees >::OverlapStart ) ... >  LeftCFEMOverlapRadii;
	typedef UIntPack< (  BSplineOverlapSizes< CDegrees , FEMDegrees >::OverlapEnd   ) ... > RightCFEMOverlapRadii;
	typedef UIntPack< ( -BSplineOverlapSizes< FEMDegrees , CDegrees >::OverlapStart ) ... >  LeftFEMCOverlapRadii;
	typedef UIntPack< (  BSplineOverlapSizes< FEMDegrees , CDegrees >::OverlapEnd   ) ... > RightFEMCOverlapRadii;

	// To set the constraints, we iterate over the splatted normals and compute the dot-product of the divergence of the normal field with all the basis functions.
	// Within the same depth: set directly as a gather 
	// Coarser depths 
	maxDepth = std::min< LocalDepth >( maxDepth , _maxDepth );
	Pointer( T ) _constraints = AllocPointer< T >( _sNodesEnd( maxDepth-1 ) );
	memset( _constraints , 0 , sizeof(T)*( _sNodesEnd(maxDepth-1) ) );
	MemoryUsage();

	static const WindowLoopData< UIntPack< BSplineOverlapSizes< CDegrees , FEMDegrees >::OverlapSize ... > > cfemLoopData( []( int c , int* start , int* end ){ BaseFEMIntegrator::ParentOverlapBounds( UIntPack< CDegrees ... >() , UIntPack< FEMDegrees ... >() , c , start , end ); } );
	static const WindowLoopData< UIntPack< BSplineOverlapSizes< FEMDegrees , CDegrees >::OverlapSize ... > > femcLoopData( []( int c , int* start , int* end ){ BaseFEMIntegrator::ParentOverlapBounds( UIntPack< FEMDegrees ... >() , UIntPack< CDegrees ... >() , c , start , end ); } );

	bool hasCoarserCoefficients = false;
	// Iterate from fine to coarse, setting the constraints @(depth) and the cumulative constraints @(depth-1)
	for( LocalDepth d=maxDepth ; d>=0 ; d-- )
	{
		typename BaseConstraint::CCStencil  stencil;
		typename BaseConstraint::PCStencils stencils;
		F.init( d );
		F.template setStencil < false >( stencil  );
		F.template setStencils< true  >( stencils );
#ifdef NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( d ) );
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
		tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d)  , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
		tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d)  , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#else // !NEW_CODE
		for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( d<maxDepth ) constraints[i] += _constraints[i];
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
			FEMTreeNode* node = _sNodes.treeNodes[i];
			int start[Dim] , end[] = { BSplineOverlapSizes< CDegrees , FEMDegrees >::OverlapSize ... };
			memset( start , 0 , sizeof( start ) );
			typename FEMTreeNode::template ConstNeighbors< OverlapSizes > neighbors;
			neighborKey.getNeighbors( LeftFEMCOverlapRadii() , RightFEMCOverlapRadii() , node , neighbors );
			bool isInterior , isInterior2;
			{
				LocalDepth d ; LocalOffset off ; _localDepthAndOffset( node , d , off );
				isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< FEMDegrees ... >() , UIntPack< CDegrees ... >() , d , off );
			}
			{
				LocalDepth d ; LocalOffset off ; _localDepthAndOffset( node->parent , d , off );
				isInterior2 = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< CDegrees ... >() , UIntPack< FEMDegrees ... >() , d , off );
			}

			LocalDepth d ; LocalOffset off;
			_localDepthAndOffset( node , d , off );

			// Set constraints from current depth
			// Gather the constraints from _node into the constraint stored with node
			if( _isValidFEM1Node( node ) )
			{
				if( isInterior )
				{
					unsigned int size = neighbors.neighbors.Size;
					Pointer( const FEMTreeNode* ) nodes = neighbors.neighbors().data;
					Pointer( Point< double , CDim > ) stencilValues = stencil.data;
					for( unsigned int j=0 ; j<size ; j++ )
					{
						if( _isValidFEM2Node( nodes[j] ) )
						{
							const D* _data = coefficients( nodes[j] );
							if( _data ) constraints[i] += _StencilDot( stencilValues[j] , *_data );
						}
					}
				}
				else
				{
					unsigned int size = neighbors.neighbors.Size;
					Pointer( const FEMTreeNode* ) nodes = neighbors.neighbors().data;
					for( unsigned int j=0 ; j<size ; j++ )
					{
						if( _isValidFEM2Node( nodes[j] ) )
						{
							const D* _data = coefficients( nodes[j] );
							if( _data )
							{
								LocalDepth _d ; LocalOffset _off ; _localDepthAndOffset( nodes[j] , _d , _off );
								constraints[i] += _StencilDot( F.ccIntegrate( off , _off ) , *_data );
							}
						}
					}
				}
				BaseFEMIntegrator::ParentOverlapBounds( UIntPack< CDegrees ... >() , UIntPack< FEMDegrees ... >() , d , off , start , end );
			}
#ifdef NEW_THREADS
			if( !_isValidFEM2Node( node ) ) return;
#else // !NEW_THREADS
			if( !_isValidFEM2Node( node ) ) continue;
#endif // NEW_THREADS
			const D* _data = coefficients( node );
#ifdef NEW_THREADS
			if( !_data ) return;
#else // !NEW_THREADS
			if( !_data ) continue;
#endif // NEW_THREADS
			else if( d<maxDepth ) hasCoarserCoefficients = true;
			const D& data = *_data;
#ifdef NEW_THREADS
			if( _IsZero( data ) ) return;
#else // !NEW_THREADS
			if( _IsZero( data ) ) continue;
#endif // NEW_THREADS

			// Set the _constraints for the parents
			if( d>0 )
			{
				int cIdx = (int)( node - node->parent->children );
				const typename BaseConstraint::CCStencil& _stencil = stencils.data[cIdx];
				neighborKey.getNeighbors( LeftCFEMOverlapRadii() , RightCFEMOverlapRadii() , node->parent , neighbors );

				unsigned int size = cfemLoopData.size[cIdx];
				const unsigned int* indices = cfemLoopData.indices[cIdx];
				ConstPointer( Point< double , CDim > ) stencilValues = _stencil.data;
				Pointer( const FEMTreeNode* ) nodes = neighbors.neighbors().data;
				if( isInterior2 )
				{
					for( unsigned int i=0 ; i<size ; i++ )
					{
						unsigned int idx = indices[i];
						if( nodes[idx] )
						{
							AddAtomic( _constraints[ nodes[idx]->nodeData.nodeIndex ] , _StencilDot( stencilValues[idx] , data ) );
						}
					}
				}
				else
				{
					for( unsigned int i=0 ; i<size ; i++ )
					{
						unsigned int idx = indices[i];
						if( nodes[idx] )
						{
							LocalDepth _d ; LocalOffset _off ; _localDepthAndOffset( nodes[idx] , _d , _off );
							AddAtomic( _constraints[ nodes[idx]->nodeData.nodeIndex ] , _StencilDot( F.pcIntegrate( _off , off ) , data ) );
						}
					}
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
#ifdef NEW_THREADS
		if( d>0 && d<maxDepth ) _downSample( UIntPack< FEMSigs ... >() , tp , F.tRestrictionProlongation() , d , _constraints );
#else // !NEW_THREADS
		if( d>0 && d<maxDepth ) _downSample( UIntPack< FEMSigs ... >() , F.tRestrictionProlongation() , d , _constraints );
#endif // NEW_THREADS
		MemoryUsage();
	}
	FreePointer( _constraints );
	if( hasCoarserCoefficients )
	{
		Pointer( D ) _coefficients = AllocPointer< D >( _sNodesEnd( maxDepth-1 ) );
		memset( _coefficients , 0 , sizeof(D) * _sNodesEnd(maxDepth-1) );
		for( LocalDepth d=maxDepth-1 ; d>=0 ; d-- )
		{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( const ThreadPool::ThreadNum & , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ ) if( _isValidFEM2Node( _sNodes.treeNodes[i] ) )
#else // !NEW_CODE
			for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ ) if( _isValidFEM2Node( _sNodes.treeNodes[i] ) )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
				const D* d = coefficients( _sNodes.treeNodes[i] );
				if( d ) _coefficients[i] += *d;
			}
#ifdef NEW_THREADS
			);
#endif // NEW_THREADS
		}

		// Coarse-to-fine up-sampling of coefficients
#ifdef NEW_THREADS
		for( LocalDepth d=1 ; d<maxDepth ; d++ ) _upSample( UIntPack< FEMSigs ... >() , tp , F.tRestrictionProlongation() , d , _coefficients );
#else // !NEW_THREADS
		for( LocalDepth d=1 ; d<maxDepth ; d++ ) _upSample( UIntPack< FEMSigs ... >() , F.tRestrictionProlongation() , d , _coefficients );
#endif // NEW_THREADS
		// Compute the contribution from all coarser depths
		for( LocalDepth d=1 ; d<=maxDepth ; d++ )
		{
#ifdef NEW_CODE
			node_index_type start = _sNodesBegin( d ) , end = _sNodesEnd( d );
			size_t range = end - start;
#else // !NEW_CODE
			size_t start = _sNodesBegin( d ) , end = _sNodesEnd( d ) , range = end - start;
#endif // NEW_CODE
			typename BaseConstraint::CPStencils stencils;
			F.init( d );
			F.template setStencils< false >( stencils );
#ifdef NEW_THREADS
			std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
			std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
			for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( d-1 ) );

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( const ThreadPool::ThreadNum & thread , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#else // !NEW_CODE
			for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
				if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
				{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
					FEMTreeNode* node = _sNodes.treeNodes[i];
					int start[Dim] , end[Dim];
					typename FEMTreeNode::template ConstNeighbors< OverlapSizes > neighbors;
					typename FEMTreeNode::template ConstNeighbors< OverlapSizes > pNeighbors;
					bool isInterior;
					{
						BaseFEMIntegrator::ParentOverlapBounds( UIntPack< FEMDegrees ... >() , UIntPack< CDegrees ... >() , (int)( node - node->parent->children ) , start , end );
					}
					{
						LocalDepth d ; LocalOffset off ; _localDepthAndOffset( node->parent , d , off );
						neighborKey.getNeighbors( LeftFEMCOverlapRadii() , RightFEMCOverlapRadii() , node->parent , pNeighbors );
						isInterior = BaseFEMIntegrator::IsInteriorlyOverlapped( UIntPack< FEMDegrees ... >() , UIntPack< CDegrees ... >() , d , off );
					}
					int cIdx = (int)( node - node->parent->children );
					const typename BaseConstraint::CCStencil& _stencil = stencils.data[cIdx];

					T constraint = {};

					LocalDepth d ; LocalOffset off;
					_localDepthAndOffset( node , d , off );
					int corner = (int)( node - node->parent->children );
					unsigned int size = femcLoopData.size[corner];
					const unsigned int* indices = femcLoopData.indices[corner];
					Pointer( const FEMTreeNode* ) nodes = pNeighbors.neighbors().data;
					Pointer( Point< double , CDim > ) stencilValues = _stencil.data;
					if( isInterior )
						for( unsigned int i=0 ; i<size ; i++ )
						{
							unsigned int idx = indices[i];
							if( _isValidFEM2Node( nodes[idx] ) ) constraint += _StencilDot( stencilValues[idx] , _coefficients[ nodes[idx]->nodeData.nodeIndex ] );
						}
					else
						for( unsigned int i=0 ; i<size ; i++ )
						{
							unsigned int idx = indices[i];
							if( _isValidFEM2Node( nodes[idx] ) )
							{
								LocalDepth _d ; LocalOffset _off ; _localDepthAndOffset ( nodes[idx] , _d , _off );
								constraint += _StencilDot( F.cpIntegrate( off , _off ) , _coefficients[ nodes[idx]->nodeData.nodeIndex ] );
							}
						}
					constraints[i] += constraint;
				}
			}
#ifdef NEW_THREADS
			);
#endif // NEW_THREADS
		}
		FreePointer( _coefficients );
	}
	MemoryUsage();
}

template< unsigned int Dim , class Real >
template< typename T , unsigned int ... FEMSigs , unsigned int PointD >
#ifdef NEW_THREADS
void FEMTree< Dim , Real >::addInterpolationConstraints( ThreadPool &tp , DenseNodeData< T , UIntPack< FEMSigs ... > >& constraints , LocalDepth maxDepth , const InterpolationInfo< T , PointD >& interpolationInfo ) const
#else // !NEW_THREADS
void FEMTree< Dim , Real >::addInterpolationConstraints( DenseNodeData< T , UIntPack< FEMSigs ... > >& constraints , LocalDepth maxDepth , const InterpolationInfo< T , PointD >& interpolationInfo ) const
#endif // NEW_THREADS
{
	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	typedef typename FEMIntegrator::template PointEvaluator< UIntPack< FEMSigs ... > , IsotropicUIntPack< Dim , PointD > > PointEvaluator;
	PointEvaluator evaluator( std::min< LocalDepth >( maxDepth , _maxDepth ) );

	typedef typename FEMTreeNode::template ConstNeighborKey< UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportEnd ... > , UIntPack< ( -BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportStart ) ...  > > PointSupportKey;
	maxDepth = std::min< LocalDepth >( maxDepth , _maxDepth );
	{
		typedef UIntPack< (-BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportStart ) ... >       LeftSupportRadii;
		typedef UIntPack< ( BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportEnd   ) ... >  LeftPointSupportRadii;
		typedef UIntPack< (-BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportStart ) ... > RightPointSupportRadii;
		typedef UIntPack<   BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize    ... > SupportSizes;

		for( int d=0 ; d<=maxDepth ; d++ )
		{
			std::vector< PointSupportKey > neighborKeys( omp_get_max_threads() );
			for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( maxDepth ) );
#ifdef NEW_THREADS
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Why do I have to wrap the template in a lambda?" )
#endif // SHOW_WARNINGS
			typedef PointEvaluatorState< UIntPack< FEMSigs ... > , IsotropicUIntPack< Dim , PointD > > _PointEvaluatorState;
			auto WrapperLambda = []( const _PointEvaluatorState &eState , LocalOffset off )
			{
				return eState.template dValues< Real , CumulativeDerivatives< Dim , PointD > >( off );
			};
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( const ThreadPool::ThreadNum & thread , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#else // !NEW_CODE
			for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
				if( _isValidSpaceNode( _sNodes.treeNodes[i] ) )
				{
					if( _isValidSpaceNode( _sNodes.treeNodes[i] ) )
					{
#ifdef NEW_THREADS
						_PointEvaluatorState eState;
#else // !NEW_THREADS
						PointEvaluatorState< UIntPack< FEMSigs ... > , IsotropicUIntPack< Dim , PointD > > eState;
#endif // NEW_THREADS
						FEMTreeNode* node = _sNodes.treeNodes[i];

#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
						PointSupportKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
						PointSupportKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
						PointSupportKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
						typename FEMTreeNode::template ConstNeighbors< SupportSizes > neighbors;
						neighborKey.getNeighbors( LeftPointSupportRadii() , RightPointSupportRadii() , node , neighbors );
						LocalDepth d ; LocalOffset off;
						_localDepthAndOffset( node , d , off );

						size_t begin , end;
						interpolationInfo.range( node , begin , end );
						for( size_t pIndex=begin ; pIndex<end ; pIndex++ )
						{
							const DualPointInfo< Dim , Real , T , PointD >& pData = interpolationInfo[ pIndex ];
							Point< Real , Dim > p = pData.position;
							evaluator.initEvaluationState( p , d , off , eState );

							int s[Dim];
							WindowLoop< Dim >::Run
							(
								IsotropicUIntPack< Dim , 0 >() , SupportSizes() ,
								[&]( int d , int i ){ s[d] = i; } ,
								[&]( const FEMTreeNode* _node )
							{
								if( _isValidFEM1Node( _node ) )
								{
									LocalDepth _d ; LocalOffset _off ; _localDepthAndOffset( _node , _d , _off );
#ifdef NEW_THREADS
									CumulativeDerivativeValues< Real , Dim , PointD > values = WrapperLambda( eState , _off );
#else // !NEW_THREADS
									CumulativeDerivativeValues< Real , Dim , PointD > values = eState.template dValues< Real , CumulativeDerivatives< Dim , PointD > >( _off );
#endif // NEW_THREADS
									T dot = {};
									for( int s=0 ; s<CumulativeDerivatives< Dim , PointD >::Size ; s++ ) dot += pData.dualValues[s] * values[s];
									AddAtomic( constraints[ _node->nodeData.nodeIndex ] , dot );
								}
							} ,
								neighbors.neighbors()
								);
						}
					}
				}
			}
#ifdef NEW_THREADS
			);
#endif // NEW_THREADS
		}
		MemoryUsage();
	}
}

template< unsigned int Dim , class Real >
template< typename T , typename TDotT , unsigned int ... FEMSigs1 , unsigned int ... FEMSigs2 , class Coefficients1 , class Coefficients2 , unsigned int PointD >
#ifdef NEW_THREADS
double FEMTree< Dim , Real >::_interpolationDot( UIntPack< FEMSigs1 ... > , UIntPack< FEMSigs2 ... > , ThreadPool &tp , const Coefficients1& coefficients1 , const Coefficients2& coefficients2 , TDotT Dot , const InterpolationInfo< T , PointD >* iInfo ) const
#else // !NEW_THREADS
double FEMTree< Dim , Real >::_interpolationDot( UIntPack< FEMSigs1 ... > , UIntPack< FEMSigs2 ... > , const Coefficients1& coefficients1 , const Coefficients2& coefficients2 , TDotT Dot , const InterpolationInfo< T , PointD >* iInfo ) const
#endif // NEW_THREADS
{
	typedef UIntPack< FEMSignature< FEMSigs1 >::Degree ... > FEMDegrees1;
	typedef UIntPack< FEMSignature< FEMSigs2 >::Degree ... > FEMDegrees2;
	typedef UIntPack< FEMSigs1 ... > FEMSignatures1;
	typedef UIntPack< FEMSigs2 ... > FEMSignatures2;
	double dot = 0;
	if( iInfo )
	{
#ifdef NEW_THREADS
		MultiThreadedEvaluator< FEMSignatures1 , PointD , T > mt1( tp , this , coefficients1 );
		MultiThreadedEvaluator< FEMSignatures2 , PointD , T > mt2( tp , this , coefficients2 );
#else // !NEW_THREADS
		MultiThreadedEvaluator< FEMSignatures1 , PointD , T > mt1( this , coefficients1 );
		MultiThreadedEvaluator< FEMSignatures2 , PointD , T > mt2( this , coefficients2 );
#endif // NEW_THREADS

		size_t begin , end;
		iInfo->range( _spaceRoot , begin , end );
#ifdef NEW_THREADS
		std::vector< double > dots( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
		tp.parallel_for( begin , end , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
		tp.parallel_for( begin , end , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : dot )
		for( int i=(int)begin ; i<(int)end ; i++ )
#endif // NEW_THREADS
		{
			Point< Real , Dim > p = (*iInfo)[i].position;
			Real w = (*iInfo)[i].weight;
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			CumulativeDerivativeValues< T , Dim , PointD > v1 = (*iInfo)( i , mt1.values( p , thread() ) );
			CumulativeDerivativeValues< T , Dim , PointD > v2 = mt2.values( p , thread() );
#else // !NEW_THREAD_NUM
			CumulativeDerivativeValues< T , Dim , PointD > v1 = (*iInfo)( i , mt1.values( p , thread ) );
			CumulativeDerivativeValues< T , Dim , PointD > v2 = mt2.values( p , thread );
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			CumulativeDerivativeValues< T , Dim , PointD > v1 = (*iInfo)( i , mt1.values( p , omp_get_thread_num() ) );
			CumulativeDerivativeValues< T , Dim , PointD > v2 = mt2.values( p , omp_get_thread_num() );
#endif // NEW_THREADS
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
			for( int dd=0 ; dd<CumulativeDerivatives< Dim , PointD >::Size ; dd++ ) dots[thread()] += Dot( v1[dd] , v2[dd] ) * w;
#else // !NEW_THREAD_NUM
			for( int dd=0 ; dd<CumulativeDerivatives< Dim , PointD >::Size ; dd++ ) dots[thread] += Dot( v1[dd] , v2[dd] ) * w;
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
			for( int dd=0 ; dd<CumulativeDerivatives< Dim , PointD >::Size ; dd++ ) dot += Dot( v1[dd] , v2[dd] ) * w;
#endif // NEW_THREADS
		}
#ifdef NEW_THREADS
		);
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) dot += dots[t];
#endif // NEW_THREADS
	}
	return dot;
}

template< unsigned int Dim , class Real >
template< typename T , typename TDotT , unsigned int ... FEMSigs1 , unsigned int ... FEMSigs2 , unsigned int ... Degrees1 , unsigned int ... Degrees2 , class Coefficients1 , class Coefficients2 >
#ifdef NEW_THREADS
double FEMTree< Dim , Real >::_dot( UIntPack< FEMSigs1 ... > , UIntPack< FEMSigs2 ... > , ThreadPool &tp , typename BaseFEMIntegrator::template Constraint< UIntPack< Degrees1 ... > , UIntPack< Degrees2 ... > , 1 >& F , const Coefficients1& coefficients1 , const Coefficients2& coefficients2 , TDotT Dot ) const
#else // !NEW_THREADS
double FEMTree< Dim , Real >::_dot( UIntPack< FEMSigs1 ... > , UIntPack< FEMSigs2 ... > , typename BaseFEMIntegrator::template Constraint< UIntPack< Degrees1 ... > , UIntPack< Degrees2 ... > , 1 >& F , const Coefficients1& coefficients1 , const Coefficients2& coefficients2 , TDotT Dot ) const
#endif // NEW_THREADS
{
	_setFEM1ValidityFlags( UIntPack< FEMSigs1 ... >() );
	_setFEM2ValidityFlags( UIntPack< FEMSigs2 ... >() );
	typedef typename BaseFEMIntegrator::template Constraint< UIntPack< Degrees1 ... > , UIntPack< Degrees2 ... > , 1 > BaseConstraint;
	double dot = 0;
	// Calculate the contribution from @(depth,depth)
	{
		typedef UIntPack<  BSplineOverlapSizes< Degrees1 , Degrees2 >::OverlapSize  ... >      OverlapSizes;
		typedef UIntPack< -BSplineOverlapSizes< Degrees1 , Degrees2 >::OverlapStart ... >  LeftOverlapRadii;
		typedef UIntPack<  BSplineOverlapSizes< Degrees1 , Degrees2 >::OverlapEnd   ... > RightOverlapRadii;

		for( LocalDepth d=0 ; d<=_maxDepth ; d++ )
		{
			typename BaseConstraint::CCStencil stencil;
			F.init( d );
			F.template setStencil< false >( stencil );

			std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
			for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( d ) );

#ifdef NEW_THREADS
			std::vector< double > dots( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : dot )
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#else // !NEW_CODE
			for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
				double &dot = dots[thread()];
#else // !NEW_THREAD_NUM
				double &dot = dots[thread];
#endif // NEW_THREAD_NUM
#endif // NEW_THREADS
				const FEMTreeNode* node = _sNodes.treeNodes[i];
				const T* _data1;
				if( _isValidFEM1Node( node ) && ( _data1=coefficients1(node) ) )
				{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
					typename FEMTreeNode::template ConstNeighbors< OverlapSizes > neighbors;
					neighborKey.getNeighbors( LeftOverlapRadii() , RightOverlapRadii() , node , neighbors );
					bool isInterior = _isInteriorlyOverlapped( UIntPack< Degrees1 ... >() , UIntPack< Degrees2 ... >() , node );

					LocalDepth d ; LocalOffset off;
					_localDepthAndOffset( node , d , off );
					ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
					ConstPointer( Point< double , 1 > ) stencilValues = stencil.data;
					if( isInterior )
					{
						for( int i=0 ; i<WindowSize< UIntPack< BSplineOverlapSizes< Degrees1 , Degrees2 >::OverlapSize ... > >::Size ; i++  )
						{
							const T* _data2;
							if( _isValidFEM2Node( nodes[i] ) && ( _data2=coefficients2( nodes[i] ) ) ) dot += Dot( *_data1 , *_data2 ) * stencilValues[i][0];
						}
					}
					else
					{
						for( int i=0 ; i<WindowSize< UIntPack< BSplineOverlapSizes< Degrees1 , Degrees2 >::OverlapSize ... > >::Size ; i++  )
						{
							const T* _data2;
							if( _isValidFEM2Node( nodes[i] ) && ( _data2=coefficients2( nodes[i] ) ) )
							{
								LocalDepth _d ; LocalOffset _off ; _localDepthAndOffset( nodes[i] , _d , _off );
								dot += Dot( *_data1 , *_data2 ) * F.ccIntegrate( off , _off )[0];
							}
						}
					}
				}
			}
#ifdef NEW_THREADS
			);
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) dot += dots[t];
#endif // NEW_THREADS
		}
	}
	// Calculate the contribution from @(<depth,depth)
	{
		typedef UIntPack<  BSplineOverlapSizes< Degrees2 , Degrees1 >::OverlapSize  ... >      OverlapSizes;
		typedef UIntPack< -BSplineOverlapSizes< Degrees2 , Degrees1 >::OverlapStart ... >  LeftOverlapRadii;
		typedef UIntPack<  BSplineOverlapSizes< Degrees2 , Degrees1 >::OverlapEnd   ... > RightOverlapRadii;

		DenseNodeData< T , UIntPack< FEMSigs1 ... > > cumulative1( _sNodesEnd( _maxDepth-1 ) );
		if( _maxDepth>0 ) memset( cumulative1() , 0 , sizeof(T) * _sNodesEnd( _maxDepth-1 ) );

		for( LocalDepth d=1 ; d<=_maxDepth ; d++ )
		{
			// Update the cumulative coefficients with the coefficients @(depth-1)
#pragma omp parallel for
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(d-1) ; i<_sNodesEnd(d-1) ; i++ )
#else // !NEW_CODE
			for( int i=_sNodesBegin(d-1) ; i<_sNodesEnd(d-1) ; i++ )
#endif // NEW_CODE
			{
				const T* _data1 = coefficients1( _sNodes.treeNodes[i] );
				if( _data1 ) cumulative1[i] += *_data1;
			}

			typename BaseConstraint::PCStencils stencils;
			F.init( d );
			F.template setStencils< true >( stencils );

			std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
			for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( d-1 ) );

#ifdef NEW_THREADS
			std::vector< double > dots( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : dot )
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#else // !NEW_CODE
			for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
				double &dot = dots[thread()];
#else // !NEW_THREAD_NUM
				double &dot = dots[thread];
#endif // NEW_THREAD_NUM
#endif // NEW_THREADS
				const FEMTreeNode* node = _sNodes.treeNodes[i];
				const T* _data2;
				if( _isValidFEM2Node( node ) && ( _data2=coefficients2( node ) ) )
				{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
					bool isInterior = _isInteriorlyOverlapped( UIntPack< Degrees1 ... >() , UIntPack< Degrees2 ... >() , node->parent );

					LocalDepth d ; LocalOffset off;
					_localDepthAndOffset( node , d , off );

					int cIdx = (int)( node - node->parent->children );
					typename BaseConstraint::CCStencil& _stencil = stencils.data[cIdx];
					typename FEMTreeNode::template ConstNeighbors< OverlapSizes > neighbors;
					neighborKey.getNeighbors( LeftOverlapRadii() , RightOverlapRadii() , node->parent , neighbors );

					int start[Dim] , end[Dim];
					_SetParentOverlapBounds( UIntPack< Degrees1 ... >() , UIntPack< Degrees2 ... >() , node , start , end );
					WindowLoop< Dim >::Run
					(
						start , end ,
						[&]( int , int ){;} ,
						[&]( const FEMTreeNode* node , Point< double , 1 > stencilValue )
						{
							const T* _data1;
							if( _isValidFEM1Node( node ) && ( _data1=cumulative1(node) ) )
							{
								if( isInterior ) dot += Dot( *_data1 , *_data2 ) * stencilValue[0];
								else
								{
									LocalDepth _d ; LocalOffset _off ; _localDepthAndOffset( node , _d , _off );
									dot += Dot ( *_data1 , *_data2 ) * F.pcIntegrate( _off , off )[0];
								}
							}
						} ,
						neighbors.neighbors() , _stencil()
					);
				}
			}
#ifdef NEW_THREADS
			);
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) dot += dots[t];
#endif // NEW_THREADS
			// Up sample the cumulative coefficients for the next level
#ifdef NEW_THREADS
			if( d<_maxDepth ) _upSample( UIntPack< FEMSigs1 ... >() , tp , F.tRestrictionProlongation() , d , cumulative1() );
#else // !NEW_THREADS
			if( d<_maxDepth ) _upSample( UIntPack< FEMSigs1 ... >() , F.tRestrictionProlongation() , d , cumulative1() );
#endif // NEW_THREADS
		}
	}

	// Calculate the contribution from @(>depth,depth)
	{
	typedef UIntPack<  BSplineOverlapSizes< Degrees1 , Degrees2 >::OverlapSize  ... >      OverlapSizes;
	typedef UIntPack< -BSplineOverlapSizes< Degrees1 , Degrees2 >::OverlapStart ... >  LeftOverlapRadii;
	typedef UIntPack<  BSplineOverlapSizes< Degrees1 , Degrees2 >::OverlapEnd   ... > RightOverlapRadii;

	DenseNodeData< T , UIntPack< FEMSigs2 ... > > cumulative2( _sNodesEnd( _maxDepth-1 ) );
	if( _maxDepth>0 ) memset( cumulative2() , 0 , sizeof(T) * _sNodesEnd( _maxDepth-1 ) );

		for( LocalDepth d=_maxDepth ; d>0 ; d-- )
		{
			typename BaseConstraint::CPStencils stencils;
			F.init( d );
			F.template setStencils< false >( stencils );

			std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
			for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( d-1 ) );

			// Update the cumulative constraints @(depth-1) from @(depth)
#ifdef NEW_THREADS
			std::vector< double > dots( tp.threadNum() , 0 );
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
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
#ifdef NEW_THREAD_NUM
				double &dot = dots[thread()];
#else // !NEW_THREAD_NUM
				double &dot = dots[thread];
#endif // NEW_THREAD_NUM
#endif // NEW_THREADS
				const FEMTreeNode* node = _sNodes.treeNodes[i];
				const T* _data1;
				if( _isValidFEM1Node( node ) && ( _data1=coefficients1( node ) ) )
				{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread() ];
#else // !NEW_THREAD_NUM
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
					ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
					bool isInterior = _isInteriorlyOverlapped( UIntPack< Degrees2 ... >() , UIntPack< Degrees1 ... >() , node->parent );

					LocalDepth d ; LocalOffset off;
					_localDepthAndOffset( node , d , off );

					int cIdx = (int)( node - node->parent->children );
					typename BaseConstraint::CCStencil& _stencil = stencils.data[cIdx];
					typename FEMTreeNode::template ConstNeighbors< OverlapSizes > neighbors;
					neighborKey.getNeighbors( LeftOverlapRadii() , RightOverlapRadii() , node->parent , neighbors );

					int start[Dim] , end[Dim];
					_SetParentOverlapBounds( UIntPack< Degrees1 ... >() , UIntPack< Degrees2 ... >() , node , start , end );

#ifdef __clang__
#pragma message ( "[WARNING] You've got me clang" )
					std::function< void (int,int) > updateFunction = [](int,int){};
#endif // __clang__

					WindowLoop< Dim >::Run
					(
						start , end ,
#ifdef __clang__
						updateFunction ,
#else // !__clang__
						[&]( int , int ){;} ,
#endif // __clang__
						[&]( const FEMTreeNode* node , Point< double , 1 > stencilValue )
						{
						if( _isValidFEM2Node( node ) )
							{
								T _dot;
								if( isInterior ) _dot = (*_data1) * stencilValue[0];
								else
								{
									LocalDepth _d ; LocalOffset _off ; _localDepthAndOffset( node , _d , _off );
									_dot = (*_data1) * F.cpIntegrate( off , _off )[0];
								}
								AddAtomic( cumulative2[ node->nodeData.nodeIndex ] , _dot );
							}
						} ,
						neighbors.neighbors() , _stencil()
					);
				}
			}
#ifdef NEW_THREADS
			);
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) dot += dots[t];
#endif // NEW_THREADS
			// Update the dot-product using the cumulative constraints @(depth-1)
#ifdef NEW_THREADS
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) dots[t] = 0;
#ifdef NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d-1) , _sNodesEnd(d-1) , [&]( const ThreadPool::ThreadNum &thread , size_t i )
#else // !NEW_THREAD_NUM
			tp.parallel_for( _sNodesBegin(d-1) , _sNodesEnd(d-1) , [&]( unsigned int thread , size_t i )
#endif // NEW_THREAD_NUM
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : dot )
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(d-1) ; i<_sNodesEnd(d-1) ; i++ )
#else // !NEW_CODE
			for( int i=_sNodesBegin(d-1) ; i<_sNodesEnd(d-1) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
#ifdef NEW_THREADS
#ifdef NEW_THREAD_NUM
				double &dot = dots[thread()];
#else // !NEW_THREAD_NUM
				double &dot = dots[thread];
#endif // NEW_THREAD_NUM
#endif // NEW_THREADS
				const FEMTreeNode* node = _sNodes.treeNodes[i];
				const T* _data2;
				if( _isValidFEM2Node( node ) && ( _data2=coefficients2( node ) ) ) dot += Dot( cumulative2[ node->nodeData.nodeIndex ] , *_data2 );
			}
#ifdef NEW_THREADS
			);
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) dot += dots[t];
#endif // NEW_THREADS

			// Down-sample the cumulative constraints from @(depth-1) to @(depth-2) for the next pass
#ifdef NEW_CODE
#ifdef NEW_THREADS
			if( d-1>0 ) _downSample( UIntPack< FEMSigs2 ... >() , tp , F.cRestrictionProlongation() , d-1 , GetPointer( &cumulative2[0] , cumulative2.size() ) );
#else // !NEW_THREADS
			if( d-1>0 ) _downSample( UIntPack< FEMSigs2 ... >() , F.cRestrictionProlongation() , d-1 , GetPointer( &cumulative2[0] , cumulative2.size() ) );
#endif // NEW_THREADS
#else // !NEW_CODE
			if( d-1>0 ) _downSample( UIntPack< FEMSigs2 ... >() , F.cRestrictionProlongation() , d-1 , GetPointer( &cumulative2[0] , (int)cumulative2.size() ) );
#endif // NEW_CODE
		}
	}

	return dot;
}
