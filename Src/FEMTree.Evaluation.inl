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

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , unsigned int PointD >
template< unsigned int _PointD >
CumulativeDerivativeValues< double , Dim , _PointD > FEMTree< Dim , Real >::_Evaluator< UIntPack< FEMSigs ... > , PointD >::_values( unsigned int d , const int fIdx[Dim] , const int cIdx[Dim] , const _CenterOffset off[Dim] , bool parentChild ) const
{
	double dValues[Dim][_PointD+1];
	_setDValues< _PointD >( d , fIdx , cIdx , off , parentChild , dValues );
	return Evaluate< Dim , double , _PointD >( dValues );
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , unsigned int PointD >
template< unsigned int _PointD >
CumulativeDerivativeValues< double , Dim , _PointD > FEMTree< Dim , Real >::_Evaluator< UIntPack< FEMSigs ... > , PointD >::_centerValues( unsigned int d , const int fIdx[Dim] , const int cIdx[Dim] , bool parentChild ) const
{
	_CenterOffset off[Dim];
	for( int d=0 ; d<Dim ; d++ ) off[d] = CENTER;
	return _values< _PointD >( d , fIdx , cIdx , off , parentChild );
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , unsigned int PointD >
template< unsigned int _PointD >
CumulativeDerivativeValues< double , Dim , _PointD > FEMTree< Dim , Real >::_Evaluator< UIntPack< FEMSigs ... > , PointD >::_cornerValues( unsigned int d , const int fIdx[Dim] , const int cIdx[Dim] , int corner , bool parentChild ) const
{
	_CenterOffset off[Dim];
	for( int d=0 ; d<Dim ; d++ ) off[d] = ( (corner>>d) & 1 ) ? FRONT : BACK;
	return _values< _PointD >( d , fIdx , cIdx , off , parentChild );
}

template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , unsigned int PointD >
void FEMTree< Dim , Real >::_Evaluator< UIntPack< FEMSigs ... > , PointD >::set( LocalDepth maxDepth )
{
	typedef UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > CenterSizes;
	static const unsigned int LeftCenterRadii[] = { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportEnd ... };
	typedef UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > CornerSizes;
	static const unsigned int LeftCornerRadii[] = { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportEnd ... };
	typedef UIntPack< ( BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::BCornerSize + 1 ) ... > BCornerSizes;
	static const unsigned int LeftBCornerRadii[] = { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::BCornerEnd ... };

	if( stencilData ) DeletePointer( stencilData );
	stencilData = NewPointer< StencilData >( maxDepth+1 );
	if( evaluators ) DeletePointer( evaluators );
	evaluators = NewPointer< Evaluators >( maxDepth+1 );
	if( childEvaluators ) DeletePointer( childEvaluators );
	childEvaluators = NewPointer< ChildEvaluators >( maxDepth+1 );
	_setEvaluators( maxDepth );
	for( int depth=0 ; depth<=maxDepth ; depth++ )
	{
		int center = ( 1<<depth )>>1;
		int cIdx[Dim] , fIdx[Dim];
		for( int d=0 ; d<Dim ; d++ ) cIdx[d] = center;

		// First set the stencils for the current depth
		{
			// The center stencil
			WindowLoop< Dim >::Run
			(
				ZeroUIntPack< Dim >() , CenterSizes() ,
				[&]( int d , int i ){ fIdx[d] = center + i - LeftCenterRadii[d]; } ,
				[&]( CumulativeDerivativeValues< double , Dim , PointD >& p ){ p = _centerValues( depth , fIdx , cIdx , false ); } ,
				stencilData[depth].ccCenterStencil()
			);
			// The corner stencil
			for( int c=0 ; c<(1<<Dim) ; c++ )
				WindowLoop< Dim >::Run
				(
					ZeroUIntPack< Dim >() , CornerSizes() ,
					[&]( int d , int i ){ fIdx[d] = center + i - LeftCornerRadii[d]; } ,
					[&]( CumulativeDerivativeValues< double , Dim , PointD >& p ){ p = _cornerValues( depth , fIdx , cIdx , c , false ); } ,
					stencilData[depth].ccCornerStencil[c]()
				);
			// The boundary corner stencil
			for( int c=0 ; c<(1<<Dim) ; c++ )
				WindowLoop< Dim >::Run
				(
					ZeroUIntPack< Dim >() , BCornerSizes() ,
					[&]( int d , int i ){ fIdx[d] = center + i - LeftBCornerRadii[d]; } ,
					[&]( CumulativeDerivativeValues< double , Dim , PointD >& p ){ p = _cornerValues( depth , fIdx , cIdx , c , false ); } ,
					stencilData[depth].ccBCornerStencil[c]()
				);
		}

		// Now set the stencils for the parents
		for( int c=0 ; c<(1<<Dim) ; c++ )
		{
			int cIdx[Dim] , fIdx[Dim];
			for( int d=0 ; d<Dim ; d++ ) cIdx[d] = center + ( (c>>d) & 1 );

			// The center stencil
			WindowLoop< Dim >::Run
			(
				ZeroUIntPack< Dim >() , CenterSizes() ,
				[&]( int d , int i ){ fIdx[d] = center/2 + i - LeftCenterRadii[d]; } ,
				[&]( CumulativeDerivativeValues< double , Dim , PointD >& p ){ p = _centerValues( depth , fIdx , cIdx , true ); } ,
				stencilData[depth].pcCenterStencils[c]()
			);
			// The corner stencil
			for( int cc=0 ; cc<(1<<Dim) ; cc++ )
				WindowLoop< Dim >::Run
				(
					ZeroUIntPack< Dim >() , CornerSizes() ,
					[&]( int d , int i ){ fIdx[d] = center/2 + i - LeftCornerRadii[d]; } ,
					[&]( CumulativeDerivativeValues< double , Dim , PointD >& p ){ p = _cornerValues( depth , fIdx , cIdx , cc , true ); } ,
					stencilData[depth].pcCornerStencils[c][cc]()
				);
			// The boundary corner stencil
			for( int cc=0 ; cc<(1<<Dim) ; cc++ )
				WindowLoop< Dim >::Run
				(
					ZeroUIntPack< Dim >() , BCornerSizes() ,
					[&]( int d , int i ){ fIdx[d] = center/2 + i - LeftBCornerRadii[d]; } ,
					[&]( CumulativeDerivativeValues< double , Dim , PointD >& p ){ p = _cornerValues( depth , fIdx , cIdx , cc , true ); } ,
					stencilData[depth].pcBCornerStencils[c][cc]()
				);
		}
	}
	if( _pointEvaluator ) delete _pointEvaluator;
	_pointEvaluator = new PointEvaluator< UIntPack< FEMSigs ... > , IsotropicUIntPack< Dim , PointD > >( maxDepth );
}

template< unsigned int Dim , class Real >
template< class V , unsigned int _PointD , unsigned int ... FEMSigs , unsigned int PointD >
CumulativeDerivativeValues< V , Dim , _PointD > FEMTree< Dim , Real >::_getValues( const ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey , const FEMTreeNode* node , Point< Real , Dim > p , ConstPointer( V ) solution , ConstPointer( V ) coarseSolution , const _Evaluator< UIntPack< FEMSigs ... > , PointD >& evaluator , int maxDepth ) const
{
	typedef UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > SupportSizes;

	if( IsActiveNode< Dim >( node->children ) && _localDepth( node->children )<=maxDepth ) WARN( "getValue assumes leaf node" );
	CumulativeDerivativeValues< V , Dim , _PointD > values;

	PointEvaluatorState< UIntPack< FEMSigs ... > , IsotropicUIntPack< Dim , _PointD > > state;

#ifdef SHOW_WARNINGS
#pragma message ( "[WARNING] Nudging evaluation point into the interior" )
#endif // SHOW_WARNINGS
	for( int dd=0 ; dd<Dim ; dd++ )
	{
		if     ( p[dd]==0 ) p[dd] = (Real)(0.+1e-6);
		else if( p[dd]==1 ) p[dd] = (Real)(1.-1e-6);
	}
	auto AddToValues = [&]( const typename FEMTreeNode::template ConstNeighbors< SupportSizes >& neighbors , ConstPointer( V ) coefficients )
	{
		ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
		for( unsigned int i=0 ; i<WindowSize< SupportSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) )
		{
			LocalDepth d ; LocalOffset off ; _localDepthAndOffset( nodes[i] , d , off );
			CumulativeDerivativeValues< Real , Dim , _PointD > _values = state.template dValues< Real , CumulativeDerivatives< Dim , _PointD > >( off );
			for( int d=0 ; d<CumulativeDerivatives< Dim , _PointD >::Size ; d++ ) values[d] += coefficients[ nodes[i]->nodeData.nodeIndex ] * _values[d];
		}
	};

	LocalDepth depth = _localDepth( node );
	while( GetGhostFlag< Dim >( node ) ) node = node->parent , depth--;

	{
		evaluator._pointEvaluator->initEvaluationState( p , depth , state );
		AddToValues( neighborKey.neighbors[ node->depth() ] , solution );
		if( depth>0 )
		{
			evaluator._pointEvaluator->initEvaluationState( p , depth-1 , state );
			AddToValues( neighborKey.neighbors[ node->parent->depth() ] , coarseSolution );
		}
	}
	// If there could be finer neighbors whose support overlaps the point
	if( depth<_maxDepth )
	{
		typename FEMTreeNode::template ConstNeighbors< SupportSizes > cNeighbors;
		int cIdx = 0;
		Point< Real , Dim > c ; Real w;
		_centerAndWidth( node , c , w );
		for( int d=0 ; d<Dim ; d++ ) if( p[d]>c[d] ) cIdx |= (1<<d);
		if( neighborKey.getChildNeighbors( cIdx , node->depth() , cNeighbors ) )
		{
			evaluator._pointEvaluator->initEvaluationState( p , depth+1 , state );
			AddToValues( cNeighbors , solution );
		}
	}
	return values;
}
template< unsigned int Dim , class Real >
template< class V , unsigned int _PointD , unsigned int ... FEMSigs , unsigned int PointD >
CumulativeDerivativeValues< V , Dim , _PointD > FEMTree< Dim , Real >::_getCenterValues( const ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey , const FEMTreeNode* node , ConstPointer( V ) solution , ConstPointer( V ) coarseSolution , const _Evaluator< UIntPack< FEMSigs ... > , PointD >& evaluator , int maxDepth , bool isInterior ) const
{
	if( IsActiveNode< Dim >( node->children ) && _localDepth( node->children )<=maxDepth ) ERROR_OUT( "getCenterValues assumes leaf node" );
	typedef _Evaluator< UIntPack< FEMSigs ... > , PointD > _Evaluator;
	typedef UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > SupportSizes;
	static const unsigned int supportSizes[] = { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... };

	if( IsActiveNode< Dim >( node->children ) && _localDepth( node->children )<=maxDepth ) ERROR_OUT( "getCenterValue assumes leaf node" );
	CumulativeDerivativeValues< V , Dim , _PointD > values;

	LocalDepth d ; LocalOffset cIdx;
	_localDepthAndOffset( node , d , cIdx );

	static const int corner = (1<<Dim)-1;

	static const CornerLoopData< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > loopData;
	auto AddToValuesInterior = [&]
	( 
		unsigned int size , const unsigned int* indices ,
		const typename FEMTreeNode::template ConstNeighbors< SupportSizes >& neighbors ,
		const typename _Evaluator::CornerStencil& cornerStencil ,
		ConstPointer( V ) coefficients
	)
	{
		ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
		ConstPointer( CumulativeDerivativeValues< double , Dim , PointD > ) _values = cornerStencil().data;
		for( unsigned int i=0 ; i<size ; i++ ) 
		{
			int idx = indices[i];
			if( IsActiveNode< Dim >( nodes[ idx ] ) ) for( int d=0 ; d<CumulativeDerivatives< Dim , _PointD >::Size ; d++ ) values[d] += coefficients[ nodes[ idx ]->nodeData.nodeIndex ] * (Real)_values[ idx ][d];
		}
	};
	auto AddToValuesExterior = [&]
	( 
		unsigned int size , const unsigned int* indices ,
		LocalDepth d , LocalOffset cIdx ,
		const typename FEMTreeNode::template ConstNeighbors< SupportSizes >& neighbors ,
		ConstPointer( V ) coefficients , bool parent
	)
	{
		ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
		for( unsigned int i=0 ; i<size ; i++ ) 
		{
			int idx = indices[i];
			if( IsActiveNode< Dim >( nodes[ idx ] ) )
			{
				LocalDepth _d ; LocalOffset fIdx;
				this->_localDepthAndOffset( nodes[idx] , _d , fIdx );
				CumulativeDerivativeValues< double , Dim , _PointD > _values = evaluator.template _cornerValues< _PointD >( d , fIdx , cIdx , corner , parent );
				for( int d=0 ; d<CumulativeDerivatives< Dim , _PointD >::Size ; d++ ) values[d] += coefficients[ nodes[ idx ]->nodeData.nodeIndex ] * (Real)_values[d];
			}
		}
	};

	if( isInterior )
	{
		auto AddToValues = [&]
		(
			const typename FEMTreeNode::template ConstNeighbors< SupportSizes >& neighbors ,
			const typename _Evaluator::CenterStencil& centerStencil ,
			ConstPointer( V ) coefficients
		)
		{
			ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
			ConstPointer( CumulativeDerivativeValues< double , Dim , PointD > ) _values = centerStencil.data;
			for( int i=0 ; i<WindowSize< SupportSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) )
				for( int d=0 ; d<CumulativeDerivatives< Dim , _PointD >::Size ; d++ ) values[d] += coefficients[ nodes[i]->nodeData.nodeIndex ] * (Real)_values[i][d];
		};
		AddToValues( neighborKey.neighbors[ node->depth() ] , evaluator.stencilData[d].ccCenterStencil , solution );
		if( d>0 )
		{
			int _corner = int( node - node->parent->children );
			AddToValues( neighborKey.neighbors[ node->parent->depth() ] , evaluator.stencilData[d].pcCenterStencils[_corner] , coarseSolution );
		}
	}
	else
	{
		auto AddToValues = [&]( const typename FEMTreeNode::template ConstNeighbors< SupportSizes >& neighbors , ConstPointer( V ) coefficients , bool parentChild )
		{
			ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
			for( int i=0 ; i<WindowSize< SupportSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) ) 
			{
				LocalDepth _d ; LocalOffset fIdx;
				_localDepthAndOffset( nodes[i] , _d , fIdx );
				const CumulativeDerivativeValues< double , Dim , _PointD >& _values = evaluator.template _centerValues< _PointD >( d , fIdx , cIdx , parentChild );
				for( int d=0 ; d<CumulativeDerivatives< Dim , _PointD >::Size ; d++ ) values[d] += coefficients[ nodes[i]->nodeData.nodeIndex ] * (Real)_values[d];
			}
		};

		AddToValues( neighborKey.neighbors[ node->depth() ] , solution , false );
		if( d>0 ) AddToValues( neighborKey.neighbors[ node->parent->depth() ] , coarseSolution , true );
	}
	// If there could be finer neighbors whose support overlaps the point
	if( d<_maxDepth )
	{
		typename FEMTreeNode::template ConstNeighbors< SupportSizes > cNeighbors;
		if( neighborKey.getChildNeighbors( 0 , node->depth() , cNeighbors ) )
		{
			if( isInterior ) AddToValuesInterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , cNeighbors , evaluator.stencilData[d+1].ccCornerStencil[corner] , solution );
			else
			{
				LocalDepth _d=d+1 ; LocalOffset _cIdx;
				for( int d=0 ; d<Dim ; d++ ) _cIdx[d] = cIdx[d]<<1;
				AddToValuesExterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , _d , _cIdx , cNeighbors , solution , false );
			}
		}
	}
	return values;
}

template< unsigned int Dim , class Real >
template< class V , unsigned int _PointD , unsigned int ... FEMSigs , unsigned int PointD >
CumulativeDerivativeValues< V , Dim , _PointD > FEMTree< Dim , Real >::_getCornerValues( const ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey , const FEMTreeNode* node , int corner , ConstPointer( V ) solution , ConstPointer( V ) coarseSolution , const _Evaluator< UIntPack< FEMSigs ... > , PointD >& evaluator , int maxDepth , bool isInterior ) const
{
	if( IsActiveNode< Dim >( node->children ) && _localDepth( node->children )<=maxDepth ) WARN( "getValue assumes leaf node" );
	typedef _Evaluator< UIntPack< FEMSigs ... > , PointD > _Evaluator;
	typedef UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > SupportSizes;
	static const unsigned int supportSizes[] = { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... };

	CumulativeDerivativeValues< V , Dim , _PointD > values;
	LocalDepth d ; LocalOffset cIdx;
	_localDepthAndOffset( node , d , cIdx );

	static const CornerLoopData< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::SupportSize ... > loopData;
	{
		auto AddToValuesInterior = [&]
		( 
			unsigned int size , const unsigned int* indices ,
			const typename FEMTreeNode::template ConstNeighbors< SupportSizes >& neighbors ,
			const typename _Evaluator::CornerStencil& cornerStencil ,
			ConstPointer( V ) coefficients
		)
		{
			ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
			ConstPointer( CumulativeDerivativeValues< double , Dim , PointD > ) _values = cornerStencil().data;
			for( unsigned int i=0 ; i<size ; i++ ) 
			{
				int idx = indices[i];
				if( IsActiveNode< Dim >( nodes[ idx ] ) ) for( int d=0 ; d<CumulativeDerivatives< Dim , _PointD >::Size ; d++ ) values[d] += coefficients[ nodes[ idx ]->nodeData.nodeIndex ] * (Real)_values[ idx ][d];
			}
		};
		auto AddToValuesExterior = [&]
		( 
			unsigned int size , const unsigned int* indices ,
			LocalDepth d , LocalOffset cIdx ,
			const typename FEMTreeNode::template ConstNeighbors< SupportSizes >& neighbors ,
			ConstPointer( V ) coefficients , bool parent
		)
		{
			ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
			for( unsigned int i=0 ; i<size ; i++ ) 
			{
				int idx = indices[i];
				if( IsActiveNode< Dim >( nodes[ idx ] ) )
				{
					LocalDepth _d ; LocalOffset fIdx;
					this->_localDepthAndOffset( nodes[idx] , _d , fIdx );
					CumulativeDerivativeValues< double , Dim , _PointD > _values = evaluator.template _cornerValues< _PointD >( d , fIdx , cIdx , corner , parent );
					for( int d=0 ; d<CumulativeDerivatives< Dim , _PointD >::Size ; d++ ) values[d] += coefficients[ nodes[ idx ]->nodeData.nodeIndex ] * (Real)_values[d];
				}
			}
		};
		if( isInterior ) AddToValuesInterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , neighborKey.neighbors[ node->depth() ] , evaluator.stencilData[d].ccCornerStencil[corner] , solution );
		else             AddToValuesExterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , d , cIdx , neighborKey.neighbors[ node->depth() ] , solution , false );
		if( d>0 )
		{
			int _corner = int( node - node->parent->children );
			if( isInterior ) AddToValuesInterior( loopData.pcSize[corner][_corner] , loopData.pcIndices[corner][_corner] , neighborKey.neighbors[ node->parent->depth() ] , evaluator.stencilData[d].pcCornerStencils[_corner][corner] , coarseSolution );
			else             AddToValuesExterior( loopData.pcSize[corner][_corner] , loopData.pcIndices[corner][_corner] , d , cIdx , neighborKey.neighbors[ node->parent->depth() ] , coarseSolution , true );
		}
		// If there could be finer neighbors whose support overlaps the point
		if( d<_maxDepth )
		{
			typename FEMTreeNode::template ConstNeighbors< SupportSizes > cNeighbors;
			if( neighborKey.getChildNeighbors( corner , node->depth() , cNeighbors ) )
			{
				if( isInterior ) AddToValuesInterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , cNeighbors , evaluator.stencilData[d+1].ccCornerStencil[corner] , solution );
				else
				{
					LocalDepth _d=d+1 ; LocalOffset _cIdx;
					for( int d=0 ; d<Dim ; d++ ) _cIdx[d] = (cIdx[d]<<1) | ( (corner&(1<<d)) ? 1 : 0 );
					AddToValuesExterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , _d , _cIdx , cNeighbors , solution , false );
				}
			}
		}
		return values;
	}
}
template< unsigned int Dim , class Real >
template< class V , unsigned int _PointD , unsigned int ... FEMSigs , unsigned int PointD >
CumulativeDerivativeValues< V , Dim , _PointD > FEMTree< Dim , Real >::_getCornerValues( const ConstCornerSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey , const FEMTreeNode* node , int corner , ConstPointer( V ) solution , ConstPointer( V ) coarseSolution , const _Evaluator< UIntPack< FEMSigs ... > , PointD >& evaluator , int maxDepth , bool isInterior ) const
{
	typedef _Evaluator< UIntPack< FEMSigs ... > , PointD > _Evaluator;

	typedef UIntPack< ( BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::BCornerSize + 1 ) ... > BCornerSizes;
	static const unsigned int bCornerSizes[] = { ( BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::BCornerSize + 1 ) ... };
	CumulativeDerivativeValues< V , Dim , _PointD > values;
	LocalDepth d ; LocalOffset cIdx;
	_localDepthAndOffset( node , d , cIdx );

	static const CornerLoopData< ( BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::BCornerSize + 1 ) ... > loopData;

	{
		auto AddToValuesInterior = [&]
		( 
			unsigned int size , const unsigned int* indices ,
			const typename FEMTreeNode::template ConstNeighbors< BCornerSizes >& neighbors ,
			const typename _Evaluator::BCornerStencil& cornerStencil ,
			ConstPointer( V ) coefficients
		)
		{
			ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
			ConstPointer( CumulativeDerivativeValues< double , Dim , PointD > ) _values = cornerStencil().data;
			for( unsigned int i=0 ; i<size ; i++ ) 
			{
				int idx = indices[i];
				if( IsActiveNode< Dim >( nodes[ idx ] ) ) for( int d=0 ; d<CumulativeDerivatives< Dim , _PointD >::Size ; d++ ) values[d] += coefficients[ nodes[ idx ]->nodeData.nodeIndex ] * (Real)_values[ idx ][d];
			}
		};
		auto AddToValuesExterior = [&]
		( 
			unsigned int size , const unsigned int* indices ,
			LocalDepth d , LocalOffset cIdx ,
			const typename FEMTreeNode::template ConstNeighbors< BCornerSizes >& neighbors ,
			ConstPointer( V ) coefficients , bool parent
		)
		{
			ConstPointer( FEMTreeNode * const ) nodes = neighbors.neighbors().data;
			for( unsigned int i=0 ; i<size ; i++ ) 
			{
				int idx = indices[i];
				if( IsActiveNode< Dim >( nodes[ idx ] ) )
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( nodes[idx] , _d , fIdx );
					CumulativeDerivativeValues< double , Dim , _PointD > _values = evaluator.template _cornerValues< _PointD >( d , fIdx , cIdx , corner , parent );
					for( int d=0 ; d<CumulativeDerivatives< Dim , _PointD >::Size ; d++ ) values[d] += coefficients[ nodes[idx]->nodeData.nodeIndex ] * (Real)_values[d];
				}
			}
		};
		if( isInterior ) AddToValuesInterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , neighborKey.neighbors[ node->depth() ] , evaluator.stencilData[d].ccBCornerStencil[corner] , solution );
		else             AddToValuesExterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , d , cIdx , neighborKey.neighbors[ node->depth() ] , solution , false );
		if( d>0 )
		{
			int _corner = int( node - node->parent->children );
			if( isInterior ) AddToValuesInterior( loopData.pcSize[corner][_corner] , loopData.pcIndices[corner][_corner] , neighborKey.neighbors[ node->parent->depth() ] , evaluator.stencilData[d].pcBCornerStencils[_corner][corner] , coarseSolution );
			else             AddToValuesExterior( loopData.pcSize[corner][_corner] , loopData.pcIndices[corner][_corner] , d , cIdx , neighborKey.neighbors[ node->parent->depth() ] , coarseSolution , true );
		}
		// If there could be finer neighbors whose support overlaps the point
		if( d<_maxDepth )
		{
			typename FEMTreeNode::template ConstNeighbors< BCornerSizes > cNeighbors;
			if( neighborKey.getChildNeighbors( corner , node->depth() , cNeighbors ) )
			{
				if( isInterior ) AddToValuesInterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , cNeighbors , evaluator.stencilData[d+1].ccBCornerStencil[corner] , solution );
				else
				{
					LocalDepth _d=d+1 ; LocalOffset _cIdx;
					for( int d=0 ; d<Dim ; d++ ) _cIdx[d] = (cIdx[d]<<1) | ( (corner&(1<<d)) ? 1 : 0 );
					AddToValuesExterior( loopData.ccSize[corner] , loopData.ccIndices[corner] , _d , _cIdx , cNeighbors , solution , false );
				}
			}
		}
	}
	return values;
}
////////////////////////////
// MultiThreadedEvaluator //
////////////////////////////
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , unsigned int PointD , typename T >
#ifdef USE_ALLOCATOR_POINTERS
#ifdef NEW_THREADS
FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::_MultiThreadedEvaluator( ConstPointer( FEMTree< Dim , Real > ) tree , const DenseNodeData< T , FEMSignatures >& coefficients , int threads ) : _coefficients( coefficients ) , _tree( tree )
#else // !NEW_THREADS
FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::_MultiThreadedEvaluator( ConstPointer( FEMTree< Dim , Real > ) tree , const DenseNodeData< T , FEMSignatures >& coefficients , int threads ) : _coefficients( coefficients ) , _tree( tree )
#endif // NEW_THREADS
#else // !USE_ALLOCATOR_POINTERS
#ifdef NEW_THREADS
FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::_MultiThreadedEvaluator( const FEMTree< Dim , Real >* tree , const DenseNodeData< T , FEMSignatures >& coefficients , int threads ) : _coefficients( coefficients ) , _tree( tree )
#else // !NEW_THREADS
FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::_MultiThreadedEvaluator( const FEMTree< Dim , Real >* tree , const DenseNodeData< T , FEMSignatures >& coefficients , int threads ) : _coefficients( coefficients ) , _tree( tree )
#endif // NEW_THREADS
#endif // USE_ALLOCATOR_POINTERS
{
	tree->_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	_threads = std::max< int >( 1 , threads );
	_pointNeighborKeys.resize( _threads );
	_cornerNeighborKeys.resize( _threads );
#ifdef NEW_THREADS
	_coarseCoefficients = _tree->template coarseCoefficients< T >( _coefficients );
#else // !NEW_THREADS
	_coarseCoefficients = _tree->template coarseCoefficients< T >( _coefficients );
#endif // NEW_THREADS
	_evaluator.set( _tree->_maxDepth );
	for( int t=0 ; t<_pointNeighborKeys.size() ; t++ ) _pointNeighborKeys[t].set( tree->_localToGlobal( _tree->_maxDepth ) );
	for( int t=0 ; t<_cornerNeighborKeys.size() ; t++ ) _cornerNeighborKeys[t].set( tree->_localToGlobal( _tree->_maxDepth ) );
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , unsigned int PointD , typename T >
template< unsigned int _PointD >
#ifdef USE_ALLOCATOR_POINTERS
CumulativeDerivativeValues< T , Dim , _PointD > FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::values( Point< Real , Dim > p , int thread , ConstPointer( FEMTreeNode ) node )
#else // !USE_ALLOCATOR_POINTERS
CumulativeDerivativeValues< T , Dim , _PointD > FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::values( Point< Real , Dim > p , int thread , const FEMTreeNode* node )
#endif // USE_ALLOCATOR_POINTERS
{
	if( _PointD>PointD ) ERROR_OUT( "Evaluating more derivatives than available: " , _PointD , " <= " , PointD );
	if( !node ) node = _tree->leaf( p );
	ConstPointSupportKey< FEMDegrees >& nKey = _pointNeighborKeys[thread];
	nKey.getNeighbors( node );
	return _tree->template _getValues< T , _PointD >( nKey , node , p , _coefficients() , _coarseCoefficients() , _evaluator , _tree->_maxDepth );
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , unsigned int PointD , typename T >
template< unsigned int _PointD >
#ifdef USE_ALLOCATOR_POINTERS
CumulativeDerivativeValues< T , Dim , _PointD > FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::centerValues( ConstPointer( FEMTreeNode ) node , int thread )
#else // !USE_ALLOCATOR_POINTERS
CumulativeDerivativeValues< T , Dim , _PointD > FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::centerValues( const FEMTreeNode* node , int thread )
#endif // USE_ALLOCATOR_POINTERS
{
	if( _PointD>PointD ) ERROR_OUT( "Evaluating more derivatives than available: " , _PointD, " <= " , PointD );
	ConstPointSupportKey< FEMDegrees >& nKey = _pointNeighborKeys[thread];
	nKey.getNeighbors( node );
	LocalDepth d ; LocalOffset off;
	_tree->_localDepthAndOffset( node->parent , d , off );
	return _tree->template _getCenterValues< T , _PointD >( nKey , node , _coefficients() , _coarseCoefficients() , _evaluator , _tree->_maxDepth , BaseFEMIntegrator::IsInteriorlySupported( UIntPack< FEMSigs ... >() , d , off ) );
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs , unsigned int PointD , typename T >
template< unsigned int _PointD >
#ifdef USE_ALLOCATOR_POINTERS
CumulativeDerivativeValues< T , Dim , _PointD > FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::cornerValues( ConstPointer( FEMTreeNode ) node , int corner , int thread )
#else // !USE_ALLOCATOR_POINTERS
CumulativeDerivativeValues< T , Dim , _PointD > FEMTree< Dim , Real >::_MultiThreadedEvaluator< UIntPack< FEMSigs ... > , PointD , T >::cornerValues( const FEMTreeNode* node , int corner , int thread )
#endif // USE_ALLOCATOR_POINTERS
{
	if( _PointD>PointD ) ERROR_OUT( "Evaluating more derivatives than available: " , _PointD , " <= " , PointD );
	ConstCornerSupportKey< FEMDegrees >& nKey = _cornerNeighborKeys[thread];
	nKey.getNeighbors( node );
	LocalDepth d ; LocalOffset off;
	_tree->_localDepthAndOffset( node->parent , d , off );
	return _tree->template _getCornerValues< T , _PointD >( nKey , node , corner , _coefficients() , _coarseCoefficients() , _evaluator , _tree->_maxDepth , BaseFEMIntegrator::IsInteriorlySupported( UIntPack< FEMSigs ... >() , d , off ) );
}



template< unsigned int Dim , class Real >
template< class V , class Coefficients , unsigned int D , unsigned int ... DataSigs >
V FEMTree< Dim , Real >::_evaluate( const Coefficients& coefficients , Point< Real , Dim > p , const PointEvaluator< UIntPack< DataSigs ... > , IsotropicUIntPack< Dim , D > >& pointEvaluator , const ConstPointSupportKey< UIntPack< FEMSignature< DataSigs >::Degree ... > >& dataKey ) const
{
	typedef UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... > SupportSizes;
	PointEvaluatorState< UIntPack< DataSigs ... > , ZeroUIntPack< Dim > > state;
	unsigned int derivatives[Dim];
	memset( derivatives , 0 , sizeof( derivatives ) );
	typedef PointSupportKey< UIntPack< FEMSignature< DataSigs >::Degree ... > > DataKey;
	V value = V();

	for( int d=_localToGlobal( 0 ) ; d<=dataKey.depth() ; d++ )
	{
		{
			const FEMTreeNode* node = dataKey.neighbors[d].neighbors.data[ WindowIndex< UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... > , UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportEnd ... > >::Index ];
			if( !node ) ERROR_OUT( "Point is not centered on a node" );
			pointEvaluator.initEvaluationState( p , _localDepth( node ) , state );
		}
		double scratch[Dim+1];
		scratch[0] = 1;
		ConstPointer( FEMTreeNode * const ) nodes = dataKey.neighbors[d].neighbors().data;
		for( int i=0 ; i<WindowSize< SupportSizes >::Size ; i++ ) if( _isValidFEM1Node( nodes[i] ) )
		{
			const V* v = coefficients( nodes[i] );
			if( v )
			{
				LocalDepth d ; LocalOffset off ; _localDepthAndOffset( nodes[i] , d , off );
				value += (*v) * (Real)state.value( off , derivatives );
			}
		}
	}

	return value;
}

template< unsigned int Dim , class Real >
template< bool XMajor , class V , unsigned int ... DataSigs >
#ifdef NEW_THREADS
Pointer( V ) FEMTree< Dim , Real >::regularGridEvaluate( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients , int& res , LocalDepth depth , bool primal ) const
#else // !NEW_THREADS
Pointer( V ) FEMTree< Dim , Real >::regularGridEvaluate( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients , int& res , LocalDepth depth , bool primal ) const
#endif // NEW_THREADS
{
	if( depth<=0 ) depth = _maxDepth;
#ifdef NEW_THREADS
	Pointer( V ) _coefficients = regularGridUpSample< XMajor >( coefficients , depth );
#else // !NEW_THREADS
	Pointer( V ) _coefficients = regularGridUpSample< XMajor >( coefficients , depth );
#endif // NEW_THREADS

	const int begin[] = { _BSplineBegin< DataSigs >( depth ) ... };
	const int end  [] = { _BSplineEnd< DataSigs >( depth ) ... };
	const int dim  [] = { ( _BSplineEnd< DataSigs >( depth ) - _BSplineBegin< DataSigs >( depth ) ) ... };

	res = 1<<depth;
	if( primal ) res++;
	size_t cellCount = 1;
	for( int d=0 ; d<Dim ; d++ ) cellCount *= res;
	Pointer( V ) values = NewPointer< V >( cellCount );
	memset( values , 0 , sizeof(V) * cellCount );

	if( primal )
	{
		// evaluate at the cell corners
		typedef UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::CornerSize ... > CornerSizes;
		typedef UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::CornerEnd ... > CornerEnds;

		EvaluationData::CornerEvaluator* evaluators[] = { ( new typename BSplineEvaluationData< DataSigs >::template CornerEvaluator< 0 >::Evaluator() ) ... };
		for( int d=0 ; d<Dim ; d++ ) evaluators[d]->set( depth );
		// Compute the offest from coefficient index to voxel index and the value of the stencil (if the voxel is interior)
		StaticWindow< long long , UIntPack< ( BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::CornerSize ? BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::CornerSize : 1 ) ... > > offsets;
		StaticWindow< double    , UIntPack< ( BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::CornerSize ? BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::CornerSize : 1 ) ... > > cornerValues;
		int dimMultiplier[Dim];
		if( XMajor )
		{
			dimMultiplier[0] = 1;
			for( int d=1 ; d<Dim ; d++ ) dimMultiplier[d] = dimMultiplier[d-1] * dim[d-1];
		}
		else
		{
			dimMultiplier[Dim-1] = 1;
			for( int d=Dim-2 ; d>=0 ; d-- ) dimMultiplier[d] = dimMultiplier[d+1] * dim[d+1];
		}

		{
			int center = ( 1<<depth )>>1;
			long long offset[Dim+1] ; offset[0] = 0;
			double upValue[Dim+1] ; upValue[0] = 1;
			WindowLoop< Dim >::Run
			(
				ZeroUIntPack< Dim >() , CornerSizes() ,
				[&]( int d , int i ){ offset[d+1] = offset[d] + ( i - (int)CornerEnds::Values[d] - begin[d] ) * dimMultiplier[d] ; upValue[d+1] = upValue[d] * evaluators[d]->value( center + i - (int)CornerEnds::Values[d] , center , false ); } ,
				[&]( long long& offsetValue , double& cornerValue ){ offsetValue = offset[Dim] , cornerValue = upValue[Dim]; } ,
				offsets() , cornerValues()
			);
		}
#ifdef NEW_THREADS
		ThreadPool::Parallel_for( 0 , cellCount , [&]( unsigned int , size_t c )
#else // !NEW_THREADS
#pragma omp parallel for
		for( long long c=0 ; c<(long long)cellCount ; c++ )
#endif // NEW_THREADS
		{
			V& value = values[c];
			int idx[Dim];
			{
				size_t _c = c;
				if( XMajor ) for( int d=0 ; d<Dim ; d++ ) idx[      d] = _c % res , _c /= res;
				else         for( int d=0 ; d<Dim ; d++ ) idx[Dim-1-d] = _c % res , _c /= res;
			}
			long long ii = 0;
			for( int d=0 ; d<Dim ; d++ ) ii += idx[d] * dimMultiplier[d];

			bool isInterior = true;
			for( int d=0 ; d<Dim ; d++ ) if( ( idx[d] - (int)CornerEnds::Values[d] )<begin[d] || ( idx[d] - (int)CornerEnds::Values[d] + (int)CornerSizes::Values[d] )>=end[d] ) isInterior = false;

			if( isInterior )
			{
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] This should be modified to support 0-degree elements" )
#endif // SHOW_WARNINGS
				ConstPointer( long long ) offsetValues = offsets().data;
				ConstPointer( double ) _cornerValues = cornerValues().data;
				for( int i=0 ; i<WindowSize< CornerSizes >::Size ; i++ ) value += _coefficients[ offsetValues[i]+ii ] * (Real)_cornerValues[i];
			}
			else
			{
				double upValues[Dim+1] ; upValues[0] = 1;	// Accumulates the product of the weights
				bool isValid[Dim+1] ; isValid[0] = true;
				WindowLoop< Dim >::Run
				(
					ZeroUIntPack< Dim >() , CornerSizes() ,
					[&]( int d , int i )
					{
						int ii = idx[d] + i - (int)CornerEnds::Values[d];
						if( ii>=begin[d] && ii<end[d] )
						{
							upValues[d+1] = upValues[d] * evaluators[d]->value( ii , idx[d] , false );
							isValid[d+1] = isValid[d];
						}
						else isValid[d+1] = false;
					} ,
					[&]( long long offsetValue ){ if( isValid[Dim] ) value += _coefficients[ offsetValue + ii ] * (Real)upValues[Dim]; } ,
					offsets()
				);
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
		for( int d=0 ; d<Dim ; d++ ) delete evaluators[d];
	}
	else
	{
		// evaluate at the cell centers
		typedef UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... > SupportSizes;
		typedef UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportEnd ... > SupportEnds;

		EvaluationData::CenterEvaluator* evaluators[] = { ( new typename BSplineEvaluationData< DataSigs >::template CenterEvaluator< 0 >::Evaluator() ) ... };
		for( int d=0 ; d<Dim ; d++ ) evaluators[d]->set( depth );
		// Compute the offest from coefficient index to voxel index and the value of the stencil (if the voxel is interior)
		StaticWindow< long long , UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... > > offsets;
		StaticWindow< double    , UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... > > centerValues;

		int dimMultiplier[Dim];
		if( XMajor )
		{
			dimMultiplier[0] = 1;
			for( int d=1 ; d<Dim ; d++ ) dimMultiplier[d] = dimMultiplier[d-1] * dim[d-1];
		}
		else
		{
			dimMultiplier[Dim-1] = 1;
			for( int d=Dim-2 ; d>=0 ; d-- ) dimMultiplier[d] = dimMultiplier[d+1] * dim[d+1];
		}

		{
			int center = ( 1<<depth )>>1;
			long long offset[Dim+1] ; offset[0] = 0;
			double upValue[Dim+1] ; upValue[0] = 1;
			WindowLoop< Dim >::Run
			(
				ZeroUIntPack< Dim >() , SupportSizes() ,
				[&]( int d , int i ){ offset[d+1] = offset[d] + ( i - (int)SupportEnds::Values[d] - begin[d] ) * dimMultiplier[d] ; upValue[d+1] = upValue[d] * evaluators[d]->value( center + i - (int)SupportEnds::Values[d] , center , false ); } ,
				[&]( long long& offsetValue , double& centerValue ){ offsetValue = offset[Dim] , centerValue = upValue[Dim]; } ,
				offsets() , centerValues()
			);
		}
#ifdef NEW_THREADS
		ThreadPool::Parallel_for( 0 , cellCount , [&]( unsigned int , size_t c )
#else // !NEW_THREADS
#pragma omp parallel for
		for( long long c=0 ; c<(long long)cellCount ; c++ )
#endif // NEW_THREADS
		{
			V& value = values[c];
			int idx[Dim];
			{
				size_t _c = c;
				if( XMajor ) for( int d=0 ; d<Dim ; d++ ) idx[      d] = _c % res , _c /= res;
				else         for( int d=0 ; d<Dim ; d++ ) idx[Dim-1-d] = _c % res , _c /= res;
			}
			long long ii = 0;
			for( int d=0 ; d<Dim ; d++ ) ii += idx[d] * dimMultiplier[d];

			bool isInterior = true;
			for( int d=0 ; d<Dim ; d++ ) if( ( idx[d] - (int)SupportEnds::Values[d] )<begin[d] || ( idx[d] - (int)SupportEnds::Values[d] + (int)SupportSizes::Values[d] )>=end[d] ) isInterior = false;

			if( isInterior )
			{
				ConstPointer( long long ) offsetValues = offsets().data;
				ConstPointer( double ) _centerValues = centerValues().data;
				for( int i=0 ; i<WindowSize< SupportSizes >::Size ; i++ ) value += _coefficients[ offsetValues[i] + ii ] * (Real)_centerValues[i];
			}
			else
			{
				double upValues[Dim+1] ; upValues[0] = 1;	// Accumulates the product of the weights
				bool isValid[Dim+1] ; isValid[0] = true;
				WindowLoop< Dim >::Run
				(
					ZeroUIntPack< Dim >() , SupportSizes() ,
					[&]( int d , int i )
					{
						int ii = idx[d] + i - (int)SupportEnds::Values[d];
						if( ii>=begin[d] && ii<end[d] )
						{
							upValues[d+1] = upValues[d] * evaluators[d]->value( ii , idx[d] , false );
							isValid[d+1] = isValid[d];
						}
						else isValid[d+1] = false;
					} ,
					[&]( long long offsetValue ){ if( isValid[Dim] ) value += _coefficients[ offsetValue + ii ] * (Real)upValues[Dim]; } ,
					offsets()
				);
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
		for( int d=0 ; d<Dim ; d++ ) delete evaluators[d];
	}
	MemoryUsage();
	DeletePointer( _coefficients );

	return values;
}
template< unsigned int Dim , class Real >
template< bool XMajor , class V , unsigned int ... DataSigs >
#ifdef NEW_THREADS
Pointer( V ) FEMTree< Dim , Real >::regularGridUpSample( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients , LocalDepth depth ) const
#else // !NEW_THREADS
Pointer( V ) FEMTree< Dim , Real >::regularGridUpSample( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients , LocalDepth depth ) const
#endif // NEW_THREADS
{
	if( depth<=0 ) depth = _maxDepth;
	int begin[Dim] , end[Dim];
	FEMIntegrator::BSplineBegin( UIntPack< DataSigs ... >() , depth , begin );
	FEMIntegrator::BSplineEnd  ( UIntPack< DataSigs ... >() , depth , end   );
#ifdef NEW_THREADS
	return regularGridUpSample< XMajor >( coefficients , begin , end , depth );
#else // !NEW_THREADS
	return regularGridUpSample< XMajor >( coefficients , begin , end , depth );
#endif // NEW_THREADS
}
template< unsigned int Dim , class Real >
template< bool XMajor , class V , unsigned int ... DataSigs >
#ifdef NEW_THREADS
Pointer( V ) FEMTree< Dim , Real >::regularGridUpSample( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients , const int begin[Dim] , const int end[Dim] , LocalDepth depth ) const
#else // !NEW_THREADS
Pointer( V ) FEMTree< Dim , Real >::regularGridUpSample( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients , const int begin[Dim] , const int end[Dim] , LocalDepth depth ) const
#endif // NEW_THREADS
{
	if( depth<=0 ) depth = _maxDepth;

	static const int DownSampleStart[][sizeof...(DataSigs)] = { { BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::DownSampleStart[0] ... } , { BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::DownSampleStart[1] ... } };
	static const int DownSampleEnd  [][sizeof...(DataSigs)] = { { BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::DownSampleEnd  [0] ... } , { BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::DownSampleEnd  [1] ... } };

	struct GridDimensions
	{
		int begin[Dim] , end[Dim] , dim[Dim];
		int dimMultiplier[Dim];
		GridDimensions( void ){ }
		GridDimensions( const int b[Dim] , const int e[Dim] )
		{
			memcpy( begin , b , sizeof(begin) );
			memcpy( end , e , sizeof(end) );
			for( int d=0 ; d<Dim ; d++ ) dim[d] = end[d] - begin[d];
			if( XMajor )
			{
				dimMultiplier[0] = 1;
				for( int d=1 ; d<Dim ; d++ ) dimMultiplier[d] = dimMultiplier[d-1] * dim[d-1];
			}
			else
			{
				dimMultiplier[Dim-1] = 1;
				for( int d=Dim-2 ; d>=0 ; d-- ) dimMultiplier[d] = dimMultiplier[d+1] * dim[d+1];
			}
		}
	};

	auto SetCoarseGridDimensions = []( const GridDimensions& fine , GridDimensions& coarse , int lowDepth )
	{
		int begin[Dim] , end[Dim];
		FEMIntegrator::BSplineBegin( UIntPack< DataSigs ... >() , lowDepth , begin );
		FEMIntegrator::BSplineEnd  ( UIntPack< DataSigs ... >() , lowDepth , end   );
		for( int d=0 ; d<Dim ; d++ )
		{
			coarse.begin[d] = std::max< int >( begin[d] , (fine.begin[d]>>1) + DownSampleStart[fine.begin[d]&1][d]   );
			coarse.end  [d] = std::min< int >( end  [d] , (fine.end  [d]>>1) + DownSampleEnd  [fine.end  [d]&1][d]+1 );
			coarse.dim  [d] = coarse.end[d] - coarse.begin[d];
		}
		if( XMajor )
		{
			coarse.dimMultiplier[0] = 1;
			for( int d=1 ; d<Dim ; d++ ) coarse.dimMultiplier[d] = coarse.dimMultiplier[d-1] * coarse.dim[d-1];
		}
		else
		{
			coarse.dimMultiplier[Dim-1] = 1;
			for( int d=Dim-2 ; d>=0 ; d-- ) coarse.dimMultiplier[d] = coarse.dimMultiplier[d+1] * coarse.dim[d+1];
		}
	};
	auto InBounds = []( const LocalOffset& off , const GridDimensions& gDim )
	{
		for( int d=0 ; d<Dim ; d++ ) if( off[d]<gDim.begin[d] || off[d]>=gDim.end[d] ) return false;
		return true;
	};

	std::vector< GridDimensions > gridDimensions( depth+1 );
	gridDimensions[depth] = GridDimensions( begin , end );
	for( int d=depth ; d>0 ; d-- ) SetCoarseGridDimensions( gridDimensions[d] , gridDimensions[d-1] , d-1 );

	// Initialize the coefficients at the coarsest level
	Pointer( V ) upSampledCoefficients = NullPointer( V );
	{
		LocalDepth _depth = 0;
		size_t count = 1;
		for( int dd=0 ; dd<Dim ; dd++ ) count *= gridDimensions[_depth].dim[dd];
		upSampledCoefficients = NewPointer< V >( count );
		memset( upSampledCoefficients , 0 , sizeof( V ) * count );
#ifdef NEW_THREADS
		ThreadPool::Parallel_for( _sNodesBegin(_depth) , _sNodesEnd(_depth) , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=_sNodesBegin(_depth) ; i<_sNodesEnd(_depth) ; i++ )
#else // !NEW_CODE
		for( int i=_sNodesBegin(_depth) ; i<_sNodesEnd(_depth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( !_outOfBounds( UIntPack< DataSigs ... >() , _sNodes.treeNodes[i] ) )
			{
				LocalDepth _d ; LocalOffset _off;
				_localDepthAndOffset( _sNodes.treeNodes[i] , _d , _off );
				if( InBounds( _off , gridDimensions[_depth] ) )
				{
					size_t idx = 0;
					for( int d=0 ; d<Dim ; d++ ) idx += gridDimensions[_depth].dimMultiplier[d] * ( _off[d] - gridDimensions[_depth].begin[d] );
					upSampledCoefficients[idx] = coefficients[i];
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
	// Up-sample and add in the existing coefficients
	for( LocalDepth _depth=1 ; _depth<=depth ; _depth++ )
	{
		size_t count = 1;
		for( int d=0 ; d<Dim ; d++ ) count *= gridDimensions[_depth].dim[d];
		Pointer( V ) _coefficients = NewPointer< V >( count );
		memset( _coefficients , 0 , sizeof( V ) * count );
		if( _depth<=_maxDepth )
#ifdef NEW_THREADS
			ThreadPool::Parallel_for( _sNodesBegin(_depth) , _sNodesEnd(_depth) , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
			for( node_index_type i=_sNodesBegin(_depth) ; i<_sNodesEnd(_depth) ; i++ )
#else // !NEW_CODE
			for( int i=_sNodesBegin(_depth) ; i<_sNodesEnd(_depth) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
				if( !_outOfBounds( UIntPack< DataSigs ... >() , _sNodes.treeNodes[i] ) )
				{
					LocalDepth _d ; LocalOffset _off;
					_localDepthAndOffset( _sNodes.treeNodes[i] , _d , _off );
					if( InBounds( _off , gridDimensions[_depth] ) )
					{
						size_t idx = 0;
						for( int d=0 ; d<Dim ; d++ ) idx += gridDimensions[_depth].dimMultiplier[d] * ( _off[d] - gridDimensions[_depth].begin[d] );
						_coefficients[idx] = coefficients[i];
					}
				}
			}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
#ifdef NEW_THREADS
		_RegularGridUpSample< XMajor >( UIntPack< DataSigs ... >() , gridDimensions[_depth-1].begin , gridDimensions[_depth-1].end , gridDimensions[_depth].begin , gridDimensions[_depth].end , _depth , ( ConstPointer(V) )upSampledCoefficients , _coefficients );
#else // !NEW_THREADS
		_RegularGridUpSample< XMajor >( UIntPack< DataSigs ... >() , gridDimensions[_depth-1].begin , gridDimensions[_depth-1].end , gridDimensions[_depth].begin , gridDimensions[_depth].end , _depth , ( ConstPointer(V) )upSampledCoefficients , _coefficients );
#endif // NEW_THREADS
		DeletePointer( upSampledCoefficients );
		upSampledCoefficients = _coefficients;
	}
	return upSampledCoefficients;
}
template< unsigned int Dim , class Real >
template< class V , unsigned int ... DataSigs >
#ifdef NEW_THREADS
V FEMTree< Dim , Real >::average( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients ) const
#else // !NEW_THREADS
V FEMTree< Dim , Real >::average( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients ) const
#endif // NEW_THREADS
{
	Real begin[Dim] , end[Dim];
	for( int d=0 ; d<Dim ; d++ ) begin[d] = (Real)0. , end[d] = (Real)1.;
#ifdef NEW_THREADS
	return average( coefficients , begin , end );
#else // !NEW_THREADS
	return average( coefficients , begin , end );
#endif // NEW_THREADS
}
template< unsigned int Dim , class Real >
template< class V , unsigned int ... DataSigs >
#ifdef NEW_THREADS
V FEMTree< Dim , Real >::average( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients , const Real begin[Dim] , const Real end[Dim] ) const
#else // !NEW_THREADS
V FEMTree< Dim , Real >::average( const DenseNodeData< V , UIntPack< DataSigs ... > >& coefficients , const Real begin[Dim] , const Real end[Dim] ) const
#endif // NEW_THREADS
{
	_setFEM1ValidityFlags( UIntPack< DataSigs ... >() );
#ifdef NEW_THREADS
	std::vector< V > avgs( ThreadPool::NumThreads() );
#else // !NEW_THREADS
	std::vector< V > avgs( omp_get_max_threads() );
#endif // NEW_THREADS
	for( int i=0 ; i<avgs.size() ; i++ ) avgs[i] = {};
	double _begin[Dim] , _end[Dim];
	for( int d=0 ; d<Dim ; d++ ) _begin[d] = begin[d] , _end[d] = end[d];
	for( int d=0 ; d<=_maxDepth ; d++ )
	{
		int center = ( 1<<d )>>1;
		int off[Dim];
		double __begin[Dim] , __end[Dim];
		for( int dd=0 ; dd<Dim ; dd++ ) off[dd] = center , __begin[dd] = 0 , __end[dd] = 1;
		double integral = FEMIntegrator::Integral( UIntPack< DataSigs ... >() , d , off , __begin , __end );
#ifdef NEW_THREADS
		ThreadPool::Parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int thread , size_t i )

#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#else //!NEW_CODE
		for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( _isValidFEM1Node( _sNodes.treeNodes[i] ) )
			{
				int d , off[Dim];
				_localDepthAndOffset( _sNodes.treeNodes[i] , d , off );
#ifdef NEW_THREADS
				if( BaseFEMIntegrator::IsInteriorlySupported( UIntPack< FEMSignature< DataSigs >::Degree ... >() , d , off , _begin , _end ) ) avgs[ thread ] += (V)( coefficients[i] * (Real)integral );
				else                                                                                                                           avgs[ thread ] += (V)( coefficients[i] * (Real)FEMIntegrator::Integral( UIntPack< DataSigs ... >() , d , off , _begin , _end ) );
#else // !NEW_THREADS
				if( BaseFEMIntegrator::IsInteriorlySupported( UIntPack< FEMSignature< DataSigs >::Degree ... >() , d , off , _begin , _end ) ) avgs[ omp_get_thread_num() ] += (V)( coefficients[i] * (Real)integral );
				else                                                                                                                           avgs[ omp_get_thread_num() ] += (V)( coefficients[i] * (Real)FEMIntegrator::Integral( UIntPack< DataSigs ... >() , d , off , _begin , _end ) );
#endif // NEW_THREADS
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
	V avg = {};
	for( int i=0 ; i<avgs.size() ; i++ ) avg += avgs[i];
	Real scale = (Real)1.;
	for( int d=0 ; d<Dim ; d++ ) scale *= end[d] - begin[d];
	return avg / scale;
}

template< unsigned int Dim , class Real >
template< unsigned int PointD , unsigned int ... FEMSigs >
#ifdef NEW_THREADS
SparseNodeData< CumulativeDerivativeValues< Real , Dim , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::leafValues( const DenseNodeData< Real , UIntPack< FEMSigs ... > >& coefficients , int maxDepth ) const
#else // !NEW_THREADS
SparseNodeData< CumulativeDerivativeValues< Real , Dim , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::leafValues( const DenseNodeData< Real , UIntPack< FEMSigs ... > >& coefficients , int maxDepth ) const
#endif // NEW_THREADS
{
	if( maxDepth<0 ) maxDepth = _maxDepth;
	_setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
	SparseNodeData< CumulativeDerivativeValues< Real , Dim , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > values;
#ifdef NEW_THREADS
	DenseNodeData< Real , UIntPack< FEMSigs ... > > _coefficients = coarseCoefficients< Real >( coefficients );
#else // !NEW_THREADS
	DenseNodeData< Real , UIntPack< FEMSigs ... > > _coefficients = coarseCoefficients< Real >( coefficients );
#endif // NEW_THREADS
	_Evaluator< UIntPack< FEMSigs ... > , PointD > evaluator;
	evaluator.set( maxDepth );
	for( LocalDepth d=maxDepth ; d>=0 ; d-- )
	{
#ifdef NEW_THREADS
		std::vector< ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > > > neighborKeys( ThreadPool::NumThreads() );
#else // !NEW_THREADS
		std::vector< ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > > > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( d ) );
#ifdef NEW_THREADS
		ThreadPool::Parallel_for( _sNodesBegin(d) , _sNodesEnd(d) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#ifdef NEW_CODE
		for( node_index_type i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#else // !NEW_CODE
		for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( _isValidSpaceNode( _sNodes.treeNodes[i] ) )
			{
#ifdef NEW_THREADS
				ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey = neighborKeys[ thread ];
#else // !NEW_THREADS
				ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
				FEMTreeNode* node = _sNodes.treeNodes[i];
				if( !IsActiveNode< Dim >( node->children ) || d==maxDepth )
				{
					neighborKey.getNeighbors( node );
					bool isInterior = _isInteriorlySupported( UIntPack< FEMSignature< FEMSigs >::Degree ... >() , node->parent );
					values[ node ] = _getCenterValues< Real , PointD >( neighborKey , node , coefficients() , _coefficients() , evaluator , maxDepth , isInterior );
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
	return values;
}
