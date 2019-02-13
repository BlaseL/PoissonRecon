
template< class T , class const_iterator > size_t SparseMatrixInterface< T , const_iterator >::entries( void ) const
{
	size_t entries = 0;
	for( size_t i=0 ; i<rows() ; i++ ) entries += rowSize( i );
	return entries;
}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::squareNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) n += iter->Value * iter->Value;
	}
	return n;

}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::squareASymmetricNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter1 = begin( i ) ; iter1!=e ; iter1++ )
		{
#ifdef NEW_CODE_SPARSE_MATRIX
			size_t j = iter1->N;
#else // !NEW_CODE_SPARSE_MATRIX
			int j = iter1->N;
#endif // NEW_CODE_SPARSE_MATRIX
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
#ifdef NEW_CODE_SPARSE_MATRIX
				size_t k = iter2->N;
#else // !NEW_CODE_SPARSE_MATRIX
				int k = iter2->N;
#endif // NEW_CODE_SPARSE_MATRIX
				if( k==i ) value += iter2->Value;
			}
			n += (iter1->Value-value) * (iter1->Value-value);
		}
	}
	return n;
}
#ifdef NEW_CODE_SPARSE_MATRIX
template< class T , class const_iterator >
template< typename Index >
double SparseMatrixInterface< T , const_iterator >::squareASymmetricNorm( Index &idx1 , Index &idx2 ) const
#else // !NEW_CODE_SPARSE_MATRIX
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::squareASymmetricNorm( int& idx1 , int& idx2 ) const
#endif // NEW_CODE_SPARSE_MATRIX
{
	double n=0;
	double max=0;
	for( size_t i=0 ; i<rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ )
		{
#ifdef NEW_CODE_SPARSE_MATRIX
			size_t j = iter->N;
#else // !NEW_CODE_SPARSE_MATRIX
			int j = iter->N;
#endif // NEW_CODE_SPARSE_MATRIX
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
#ifdef NEW_CODE_SPARSE_MATRIX
				size_t k = iter2->N;
#else // !NEW_CODE_SPARSE_MATRIX
				int k = iter2->N;
#endif // NEW_CODE_SPARSE_MATRIX
				if( k==i ) value += iter2->Value;
			}
			double temp = (iter->Value-value) * (iter->Value-value);
			n += temp;
#ifdef NEW_CODE_SPARSE_MATRIX
			if( temp>=max ) idx1 = (Index)i , idx2 = (Index)j , max = temp;
#else // !NEW_CODE_SPARSE_MATRIX
			if( temp>=max ) idx1 = i , idx2 = j , max=temp;
#endif // NEW_CODE_SPARSE_MATRIX
		}
	}
	return n;
}
template< class T , class const_iterator >
template< class T2 >
#ifdef NEW_CODE_SPARSE_MATRIX
void SparseMatrixInterface< T , const_iterator >::multiply( ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag ) const
#else // !NEW_CODE_SPARSE_MATRIX
void SparseMatrixInterface< T , const_iterator >::multiply( ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag ) const
#endif // NEW_CODE_SPARSE_MATRIX
{
	ConstPointer( T2 ) in = In;
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)rows() ; i++ )
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<rows() ; i++ )
#endif // NEW_CODE_SPARSE_MATRIX
	{
		T2 temp;
		memset( &temp , 0 , sizeof(T2) );
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += (T2)( _in[ iter->N ] * iter->Value );
		if( multiplyFlag & MULTIPLY_NEGATE ) temp = -temp;
		if( multiplyFlag & MULTIPLY_ADD ) Out[i] += temp;
		else                              Out[i]  = temp;
	}
}
template< class T , class const_iterator >
template< class T2 >
#ifdef NEW_CODE_SPARSE_MATRIX
void SparseMatrixInterface< T , const_iterator >::multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag ) const
#else // !NEW_CODE_SPARSE_MATRIX
void SparseMatrixInterface< T , const_iterator >::multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag ) const
#endif // NEW_CODE_SPARSE_MATRIX
{
	ConstPointer( T2 ) in = In;
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)rows() ; i++ )
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<rows() ; i++ )
#endif // NEW_CODE_SPARSE_MATRIX
	{
		T2 temp;
		memset( &temp , 0 , sizeof(T2) );
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * iter->Value;
		temp *= scale;
		if( multiplyFlag & MULTIPLY_NEGATE ) temp = -temp;
		if( multiplyFlag & MULTIPLY_ADD ) Out[i] += temp;
		else                              Out[i]  = temp;
	}
}

template< class T , class const_iterator >
void SparseMatrixInterface< T , const_iterator >::setDiagonal( Pointer( T ) diagonal ) const
{
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)rows() ; i++ )
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<rows() ; i++ )
#endif // NEW_CODE_SPARSE_MATRIX
	{
		diagonal[i] = (T)0;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) if( iter->N==i ) diagonal[i] += iter->Value;
	}
}

template< class T , class const_iterator >
void SparseMatrixInterface< T , const_iterator >::setDiagonalR( Pointer( T ) diagonal ) const
{
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)rows() ; i++ )
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<rows() ; i++ )
#endif // NEW_CODE_SPARSE_MATRIX
	{
		diagonal[i] = (T)0;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) if( iter->N==i ) diagonal[i] += iter->Value;
		if( diagonal[i] ) diagonal[i] = (T)( 1./diagonal[i] );
	}
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::jacobiIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , ConstPointer( T2 ) in , Pointer( T2 ) out , bool dReciprocal ) const
{
	multiply( in , out );
	if( dReciprocal )
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) * diagonal[i];
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) * diagonal[i];
#endif // NEW_CODE_SPARSE_MATRIX
	else
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) / diagonal[i];
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) / diagonal[i];
#endif // NEW_CODE_SPARSE_MATRIX
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const
{
	if( dReciprocal )
	{
#define ITERATE( j )                                                                                \
	{                                                                                               \
		T2 _b = b[j];                                                                               \
		const_iterator e = end( j );                                                                \
		for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
		x[j] += _b * diagonal[j];                                                                   \
	}
#ifdef NEW_CODE_SPARSE_MATRIX
		if( forward ) for( long long j=0 ; j<(long long)rows() ; j++ ){ ITERATE( j ); }
		else          for( long long j=(long long)rows()-1 ; j>=0 ; j-- ){ ITERATE( j ); }
#else // !NEW_CODE_SPARSE_MATRIX
		if( forward ) for( int j=0 ; j<int( rows() ) ; j++ ){ ITERATE( j ); }
		else          for( int j=int( rows() )-1 ; j>=0 ; j-- ){ ITERATE( j ); }
#endif // NEW_CODE_SPARSE_MATRIX
#undef ITERATE
	}
	else
	{
#define ITERATE( j )                                                                                \
	{                                                                                               \
		T2 _b = b[j];                                                                               \
		const_iterator e = end( j );                                                                \
		for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
		x[j] += _b / diagonal[j];                                                                   \
	}

#ifdef NEW_CODE_SPARSE_MATRIX
		if( forward ) for( long long j=0 ; j<(long long)rows() ; j++ ){ ITERATE( j ); }
		else          for( long long j=(long long)rows()-1 ; j>=0 ; j-- ){ ITERATE( j ); }
#else // !NEW_CODE_SPARSE_MATRIX
		if( forward ) for( int j=0 ; j<int( rows() ) ; j++ ){ ITERATE( j ); }
		else          for( int j=int( rows() )-1 ; j>=0 ; j-- ){ ITERATE( j ); }
#endif // NEW_CODE_SPARSE_MATRIX
#undef ITERATE
	}
}
template< class T , class const_iterator >
#ifdef NEW_CODE_SPARSE_MATRIX
template< typename Index , class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< Index >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool dReciprocal ) const
#else // !NEW_CODE_SPARSE_MATRIX
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< int >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool dReciprocal ) const
#endif // NEW_CODE_SPARSE_MATRIX
{
	if( dReciprocal )
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long j=0 ; j<(long long)multiColorIndices.size() ; j++ )
#else // !NEW_CODE_SPARSE_MATRIX
		for( int j=0 ; j<(int)multiColorIndices.size() ; j++ )
#endif // NEW_CODE_SPARSE_MATRIX
		{
#ifdef NEW_CODE_SPARSE_MATRIX
			Index jj = multiColorIndices[j];
#else // !NEW_CODE_SPARSE_MATRIX
			int jj = multiColorIndices[j];
#endif // NEW_CODE_SPARSE_MATRIX
			T2 _b = b[jj];
			const_iterator e = end( jj );
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;
			x[jj] += _b * diagonal[jj];
		}
	else
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long j=0 ; j<(long long)multiColorIndices.size() ; j++ )
#else // !NEW_CODE_SPARSE_MATRIX
		for( int j=0 ; j<(int)multiColorIndices.size() ; j++ )
#endif // NEW_CODE_SPARSE_MATRIX
		{
#ifdef NEW_CODE_SPARSE_MATRIX
			Index jj = multiColorIndices[j];
#else // !NEW_CODE_SPARSE_MATRIX
			int jj = multiColorIndices[j];
#endif // NEW_CODE_SPARSE_MATRIX
			T2 _b = b[jj];
			const_iterator e = end( jj );
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;
			x[jj] += _b / diagonal[jj];
		}
}

template< class T , class const_iterator >
#ifdef NEW_CODE_SPARSE_MATRIX
template< typename Index , class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< std::vector< Index > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const
#else // !NEW_CODE_SPARSE_MATRIX
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const
#endif // NEW_CODE_SPARSE_MATRIX
{
#ifdef _WIN32
#define SetOMPParallel __pragma( omp parallel for )
#else // !_WIN32
#define SetOMPParallel _Pragma( "omp parallel for" )
#endif // _WIN32

#ifdef NEW_CODE_SPARSE_MATRIX
	if( dReciprocal )
	{
#define ITERATE( indices )                                                                               \
	{                                                                                                    \
SetOMPParallel                                                                                           \
		for( long long k=0 ; k<(long long)indices.size() ; k++ )                                       \
		{                                                                                                \
			long long jj = indices[k];                                                                   \
			T2 _b = b[jj];                                                                               \
			const_iterator e = end( jj );                                                                \
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[jj] += _b * diagonal[jj];                                                                  \
		}                                                                                                \
	}
		if( forward ) for( long long j=0 ; j<(long long)multiColorIndices.size() ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( long long j=(long long)multiColorIndices.size()-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
	else
	{
#define ITERATE( indices )                                                                               \
	{                                                                                                    \
SetOMPParallel                                                                                           \
		for( long long k=0 ; k<(long long)indices.size() ; k++ )                                       \
		{                                                                                                \
			long long jj = indices[k];                                                                   \
			T2 _b = b[jj];                                                                               \
			const_iterator e = end( jj );                                                                \
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[jj] += _b / diagonal[jj];                                                                  \
		}                                                                                                \
	}
		if( forward ) for( long long j=0 ; j<(long long)multiColorIndices.size()  ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( long long j=(long long)multiColorIndices.size()-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
#undef SetOMPParallel
#else // !NEW_CODE_SPARSE_MATRIX
	if( dReciprocal )
	{
#define ITERATE( indices )                                                                               \
	{                                                                                                    \
SetOMPParallel                                                                                           \
		for( int k=0 ; k<int( indices.size() ) ; k++ )                                                   \
		{                                                                                                \
			int jj = indices[k];                                                                         \
			T2 _b = b[jj];                                                                               \
			const_iterator e = end( jj );                                                                \
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[jj] += _b * diagonal[jj];                                                                  \
		}                                                                                                \
	}
		if( forward ) for( int j=0 ; j<multiColorIndices.size()  ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( int j=int( multiColorIndices.size() )-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
	else
	{
#define ITERATE( indices )                                                                               \
	{                                                                                                    \
SetOMPParallel                                                                                           \
		for( int k=0 ; k<int( indices.size() ) ; k++ )                                                   \
		{                                                                                                \
			int jj = indices[k];                                                                         \
			T2 _b = b[jj];                                                                               \
			const_iterator e = end( jj );                                                                \
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[jj] += _b / diagonal[jj];                                                                  \
		}                                                                                                \
	}
		if( forward ) for( int j=0 ; j<multiColorIndices.size()  ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( int j=int( multiColorIndices.size() )-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
#undef SetOMPParallel
#endif // NEW_CODE_SPARSE_MATRIX
}
#ifdef NEW_CODE_SPARSE_MATRIX
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , size_t dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot )
#else // !NEW_CODE_SPARSE_MATRIX
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot )
#endif // NEW_CODE_SPARSE_MATRIX
{
	eps *= eps;
	Pointer( T ) r = AllocPointer< T >( dim );
	Pointer( T ) d = AllocPointer< T >( dim );
	Pointer( T ) q = AllocPointer< T >( dim );

	Real delta_new = 0 , delta_0;
	M( ( ConstPointer( T ) )x , r );
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)dim ; i++ ) d[i] = r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] );
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<dim ; i++ ) d[i] = r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] );
#endif // NEW_CODE_SPARSE_MATRIX

	delta_0 = delta_new;
	if( delta_new<=eps )
	{
		FreePointer( r );
		FreePointer( d );
		FreePointer( q );
		return 0;
	}
	int ii;
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		M( ( ConstPointer( T ) )d , q );
		Real dDotQ = 0;
#pragma omp parallel for reduction( + : dDotQ )
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
#endif // NEW_CODE_SPARSE_MATRIX
		if( !dDotQ ) break;

		Real alpha = delta_new / dDotQ;
		Real delta_old = delta_new;
		delta_new = 0;
		if( (ii%50)==(50-1) )
		{
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
			M( ( ConstPointer( T ) )x , r );
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
		}
		else
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX

		Real beta = delta_new / delta_old;
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
#endif // NEW_CODE_SPARSE_MATRIX
	}
	FreePointer( r );
	FreePointer( d );
	FreePointer( q );
	return ii;
}
#ifdef NEW_CODE_SPARSE_MATRIX
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , const Preconditioner& P , size_t dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot  )
#else // !NEW_CODE_SPARSE_MATRIX
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , const Preconditioner& P , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot  )
#endif // NEW_CODE_SPARSE_MATRIX
{
	eps *= eps;
	Pointer( T ) r = AllocPointer< T >( dim );
	Pointer( T ) d = AllocPointer< T >( dim );
	Pointer( T ) q = AllocPointer< T >( dim );
	Pointer( T ) Pb = AllocPointer< T >( dim );
	Pointer( T ) temp = AllocPointer< T >( dim );

	auto PM = [&] ( ConstPointer(T) x , Pointer(T) y )
	{
		M( x , temp );
		P( ( ConstPointer(T) )temp , y );
	};

	Real delta_new = 0 , delta_0;
	P( b , Pb );
	PM( ( ConstPointer( T ) )x , r );
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)dim ; i++ ) d[i] = r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] );
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<dim ; i++ ) d[i] = r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] );
#endif // NEW_CODE_SPARSE_MATRIX

	delta_0 = delta_new;
	if( delta_new<=eps )
	{
		FreePointer( Pb );
		FreePointer( r );
		FreePointer( d );
		FreePointer( q );
		FreePointer( temp );
		return 0;
	}
	int ii;
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		PM( ( ConstPointer( T ) )d , q );
		Real dDotQ = 0;
#pragma omp parallel for reduction( + : dDotQ )
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
#endif //NEW_CODE_SPARSE_MATRIX
		if( !dDotQ ) break;

		Real alpha = delta_new / dDotQ;
		Real delta_old = delta_new;
		delta_new = 0;
		if( (ii%50)==(50-1) )
		{
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
			PM( ( ConstPointer( T ) )x , r );
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
		}
		else
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX

		Real beta = delta_new / delta_old;
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
#endif // NEW_CODE_SPARSE_MATRIX
	}
	FreePointer( Pb );
	FreePointer( r );
	FreePointer( d );
	FreePointer( q );
	FreePointer( temp );
	return ii;
}
