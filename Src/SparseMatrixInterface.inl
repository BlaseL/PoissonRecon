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
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::squareASymmetricNorm( size_t &idx1 , size_t &idx2 ) const
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
			if( temp>=max ) idx1 = i , idx2 = j , max=temp;
		}
	}
	return n;
}
template< class T , class const_iterator >
template< class T2 >
#ifdef NEW_CODE_SPARSE_MATRIX
#ifdef NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::multiply( ThreadPool &tp , ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag ) const
#else // !NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::multiply( ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag ) const
#endif // NEW_THREADS
#else // !NEW_CODE_SPARSE_MATRIX
void SparseMatrixInterface< T , const_iterator >::multiply( ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag ) const
#endif // NEW_CODE_SPARSE_MATRIX
{
	ConstPointer( T2 ) in = In;
#ifdef NEW_THREADS
	tp.parallel_for( 0 , rows() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)rows() ; i++ )
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<rows() ; i++ )
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
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
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}
template< class T , class const_iterator >
template< class T2 >
#ifdef NEW_CODE_SPARSE_MATRIX
#ifdef NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::multiplyScaled( ThreadPool &tp , T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag ) const
#else // !NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag ) const
#endif // NEW_THREADS
#else // !NEW_CODE_SPARSE_MATRIX
void SparseMatrixInterface< T , const_iterator >::multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag ) const
#endif // NEW_CODE_SPARSE_MATRIX
{
	ConstPointer( T2 ) in = In;
#ifdef NEW_THREADS
	tp.parallel_for( 0 , rows() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)rows() ; i++ )
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<rows() ; i++ )
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
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
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< class T , class const_iterator >
#ifdef NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::setDiagonal( ThreadPool &tp , Pointer( T ) diagonal ) const
#else // !NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::setDiagonal( Pointer( T ) diagonal ) const
#endif // NEW_THREADS
{
#ifdef NEW_THREADS
	tp.parallel_for( 0 , rows() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)rows() ; i++ )
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<rows() ; i++ )
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
	{
		diagonal[i] = (T)0;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) if( iter->N==i ) diagonal[i] += iter->Value;
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< class T , class const_iterator >
#ifdef NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::setDiagonalR( ThreadPool &tp , Pointer( T ) diagonal ) const
#else // !NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::setDiagonalR( Pointer( T ) diagonal ) const
#endif // NEW_THREADS
{
#ifdef NEW_THREADS
	tp.parallel_for( 0 , rows() , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)rows() ; i++ )
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<rows() ; i++ )
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
	{
		diagonal[i] = (T)0;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) if( iter->N==i ) diagonal[i] += iter->Value;
		if( diagonal[i] ) diagonal[i] = (T)( 1./diagonal[i] );
	}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< class T , class const_iterator >
template< class T2 >
#ifdef NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::jacobiIteration( ThreadPool &tp , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , ConstPointer( T2 ) in , Pointer( T2 ) out , bool dReciprocal ) const
#else // !NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::jacobiIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , ConstPointer( T2 ) in , Pointer( T2 ) out , bool dReciprocal ) const
#endif // NEW_THREADS
{
#ifdef NEW_THREADS
	multiply( tp , in , out );
#else // !NEW_THREADS
	multiply( in , out );
#endif // NEW_THREADS
	if( dReciprocal )
#ifdef NEW_THREADS
		tp.parallel_for( 0 , rows() , [&]( unsigned int , size_t i ){ out[i] = in[i] + ( b[i] - out[i] ) * diagonal[i]; } );
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) * diagonal[i];
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) * diagonal[i];
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
	else
#ifdef NEW_THREADS
		tp.parallel_for( 0 , rows() , [&]( unsigned int , size_t i ){ out[i] = in[i] + ( b[i] - out[i] ) / diagonal[i]; } );
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) / diagonal[i];
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) / diagonal[i];
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
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
template< class T2 >
#ifdef NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::gsIteration( ThreadPool &tp , const std::vector< size_t >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool dReciprocal ) const
#else // !NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< size_t >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool dReciprocal ) const
#endif // NEW_THREADS
#else // !NEW_CODE_SPARSE_MATRIX
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< int >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool dReciprocal ) const
#endif // NEW_CODE_SPARSE_MATRIX
{
	if( dReciprocal )
#ifdef NEW_THREADS
		tp.parallel_for( 0 , multiColorIndices.size() , [&]( unsigned int , size_t j )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long j=0 ; j<(long long)multiColorIndices.size() ; j++ )
#else // !NEW_CODE_SPARSE_MATRIX
		for( int j=0 ; j<(int)multiColorIndices.size() ; j++ )
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
		{
#ifdef NEW_CODE_SPARSE_MATRIX
			size_t jj = multiColorIndices[j];
#else // !NEW_CODE_SPARSE_MATRIX
			int jj = multiColorIndices[j];
#endif // NEW_CODE_SPARSE_MATRIX
			T2 _b = b[jj];
			const_iterator e = end( jj );
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;
			x[jj] += _b * diagonal[jj];
		}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
	else
#ifdef NEW_THREADS
		tp.parallel_for( 0 , multiColorIndices.size() , [&]( unsigned int , size_t j )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long j=0 ; j<(long long)multiColorIndices.size() ; j++ )
#else // !NEW_CODE_SPARSE_MATRIX
		for( int j=0 ; j<(int)multiColorIndices.size() ; j++ )
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
		{
#ifdef NEW_CODE_SPARSE_MATRIX
			size_t jj = multiColorIndices[j];
#else // !NEW_CODE_SPARSE_MATRIX
			int jj = multiColorIndices[j];
#endif // NEW_CODE_SPARSE_MATRIX
			T2 _b = b[jj];
			const_iterator e = end( jj );
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;
			x[jj] += _b / diagonal[jj];
		}
#ifdef NEW_THREADS
	);
#endif // NEW_THREADS
}

template< class T , class const_iterator >
#ifdef NEW_CODE_SPARSE_MATRIX
template< class T2 >
#ifdef NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::gsIteration( ThreadPool &tp , const std::vector< std::vector< size_t > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const
#else // !NEW_THREADS
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< std::vector< size_t > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const
#endif // NEW_THREADS
#else // !NEW_CODE_SPARSE_MATRIX
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const
#endif // NEW_CODE_SPARSE_MATRIX
{
#ifdef NEW_THREADS
#ifdef FIXED_BLOCK_SIZE
#else // !FIXED_BLOCK_SIZE
	size_t blockSize = tp.blockSize();
	tp.setBlockSize( 400 );
#endif // FIXED_BLOCK_SIZE
	if( dReciprocal )
	{
#define ITERATE( indices )                                                                               \
	{                                                                                                    \
		tp.parallel_for( 0 , indices.size() , [&]( unsigned int , size_t k )                             \
		{                                                                                                \
			size_t jj = indices[k];                                                                      \
			T2 _b = b[jj];                                                                               \
			const_iterator e = end( jj );                                                                \
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[jj] += _b * diagonal[jj];                                                                  \
		}                                                                                                \
		);                                                                                               \
	}
		if( forward ) for( size_t j=0 ; j<multiColorIndices.size() ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( size_t j=multiColorIndices.size()-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
	else
	{
#define ITERATE( indices )                                                                               \
	{                                                                                                    \
		tp.parallel_for( 0 , indices.size() , [&]( unsigned int , size_t k )                             \
		{                                                                                                \
			size_t jj = indices[k];                                                                      \
			T2 _b = b[jj];                                                                               \
			const_iterator e = end( jj );                                                                \
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[jj] += _b / diagonal[jj];                                                                  \
		}                                                                                                \
		);                                                                                               \
	}
		if( forward ) for( size_t j=0 ; j<multiColorIndices.size()  ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( size_t j=multiColorIndices.size()-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
#ifdef FIXED_BLOCK_SIZE
#else // !FIXED_BLOCK_SIZE
	tp.setBlockSize( blockSize );
#endif // FIXED_BLOCK_SIZE
#else // !NEW_THREADS
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
		for( long long k=0 ; k<(long long)indices.size() ; k++ )                                         \
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
		for( long long k=0 ; k<(long long)indices.size() ; k++ )                                         \
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
#endif // NEW_THREADS

}
#ifdef NEW_CODE_SPARSE_MATRIX
#ifdef NEW_THREADS
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > size_t SolveCG( ThreadPool &tp , const SPDFunctor& M , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot )
#else // !NEW_THREADS
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > size_t SolveCG( const SPDFunctor& M , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot )
#endif // NEW_THREADS
#else // !NEW_CODE_SPARSE_MATRIX
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot )
#endif // NEW_CODE_SPARSE_MATRIX
{
#ifdef NEW_THREADS
	std::vector< Real > scratch( tp.threadNum() , 0 );
#endif // NEW_THREADS
	eps *= eps;
	Pointer( T ) r = AllocPointer< T >( dim );
	Pointer( T ) d = AllocPointer< T >( dim );
	Pointer( T ) q = AllocPointer< T >( dim );

	Real delta_new = 0 , delta_0;
	M( ( ConstPointer( T ) )x , r );
#ifdef NEW_THREADS
	tp.parallel_for( 0 , dim , [&]( unsigned int thread , size_t i ){ d[i] = r[i] = b[i] - r[i] , scratch[thread] += Dot( r[i] , r[i] ); } );
	for( unsigned int t=0 ; t<tp.threadNum() ; t++ ){ delta_new += scratch[t] ; scratch[t] = 0; }
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)dim ; i++ ) d[i] = r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] );
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<dim ; i++ ) d[i] = r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS

	delta_0 = delta_new;
	if( delta_new<=eps )
	{
		FreePointer( r );
		FreePointer( d );
		FreePointer( q );
		return 0;
	}
#ifdef NEW_CODE
	size_t ii;
#else // !NEW_CODE
	int ii;
#endif // NEW_CODE
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		M( ( ConstPointer( T ) )d , q );
		Real dDotQ = 0;
#ifdef NEW_THREADS
		tp.parallel_for( 0 , dim , [&]( unsigned int thread , size_t i ){ scratch[thread] += Dot( d[i] , q[i] ); } );
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ){ dDotQ += scratch[t] ; scratch[t] = 0; }
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : dDotQ )
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
		if( !dDotQ ) break;

		Real alpha = delta_new / dDotQ;
		Real delta_old = delta_new;
		delta_new = 0;
		if( (ii%50)==(50-1) )
		{
#ifdef NEW_THREADS
			tp.parallel_for( 0 , dim , [&]( unsigned int , size_t i ){ x[i] += (T)( d[i] * alpha ); } );
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
			M( ( ConstPointer( T ) )x , r );
#ifdef NEW_THREADS
			tp.parallel_for( 0 , dim , [&]( unsigned int thread , size_t i ){ r[i] = b[i] - r[i] , scratch[thread] += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha ); } );
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ){ delta_new += scratch[t] ; scratch[t] = 0; }
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
		}
		else
#ifdef NEW_THREADS
		{
			tp.parallel_for( 0 , dim , [&]( unsigned int thread , size_t i ){ r[i] -=(T)( q[i] * alpha ) , scratch[thread] += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha ); } );
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ){ delta_new += scratch[t] ; scratch[t] = 0; }
		}
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS

		Real beta = delta_new / delta_old;
#ifdef NEW_THREADS
		tp.parallel_for( 0 , dim , [&]( unsigned int , size_t i ){ d[i] = r[i] + (T)( d[i] * beta ); } );
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
	}
	FreePointer( r );
	FreePointer( d );
	FreePointer( q );
	return ii;
}
#ifdef NEW_CODE_SPARSE_MATRIX
#ifdef NEW_THREADS
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > size_t SolveCG( ThreadPool &tp , const SPDFunctor& M , const Preconditioner& P , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot  )
#else // !NEW_THREADS
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > size_t SolveCG( const SPDFunctor& M , const Preconditioner& P , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot  )
#endif // NEW_THREADS
#else // !NEW_CODE_SPARSE_MATRIX
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , const Preconditioner& P , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot  )
#endif // NEW_CODE_SPARSE_MATRIX
{
#ifdef NEW_THREADS
	std::vector< Real > scratch( tp.threadNum() , 0 );
#endif // NEW_THREADS
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
#ifdef NEW_THREADS
	tp.parallel_for( 0 , dim , [&]( unsigned int thread , size_t i ){ d[i] = r[i] = Pb[i] - r[i] , scratch[thread] += Dot( r[i] , r[i] ); } );
	for( unsigned int t=0 ; t<tp.threadNum() ; t++ ){ delta_new += scratch[t] ; scratch[t] = 0; }
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
	for( long long i=0 ; i<(long long)dim ; i++ ) d[i] = r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] );
#else // !NEW_CODE_SPARSE_MATRIX
	for( int i=0 ; i<dim ; i++ ) d[i] = r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS

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
#ifdef NEW_CODE
	size_t ii;
#else // !NEW_CODE
	int ii;
#endif // NEW_CODE
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		PM( ( ConstPointer( T ) )d , q );
		Real dDotQ = 0;
#ifdef NEW_THREADS
		tp.parallel_for( 0 , dim , [&]( unsigned int thread , size_t i ){ scratch[thread] += Dot( d[i] , q[i] ); } );
		for( unsigned int t=0 ; t<tp.threadNum() ; t++ ){ dDotQ += scratch[t] ; scratch[t] = 0; }
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : dDotQ )
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
#endif //NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
		if( !dDotQ ) break;

		Real alpha = delta_new / dDotQ;
		Real delta_old = delta_new;
		delta_new = 0;
		if( (ii%50)==(50-1) )
		{
#ifdef NEW_THREADS
			tp.parallel_for( 0 , dim , [&]( unsigned int , size_t i ){ x[i] += (T)( d[i] * alpha ); } );
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
			PM( ( ConstPointer( T ) )x , r );
#ifdef NEW_THREADS
			tp.parallel_for( 0 , dim , [&]( unsigned int thread , size_t i ){ r[i] = Pb[i] - r[i] , scratch[thread] += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha ); } );
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) delta_new += scratch[t];
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
		}
		else
		{
#ifdef NEW_THREADS
			tp.parallel_for( 0 , dim , [&]( unsigned int thread , size_t i ){ r[i] -=(T)( q[i] * alpha ) , scratch[thread] += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha ); } );
			for( unsigned int t=0 ; t<tp.threadNum() ; t++ ) delta_new += scratch[t];
#else // !NEW_THREADS
#pragma omp parallel for reduction( + : delta_new )
#ifdef NEW_CODE_SPARSE_MATRIX
			for( long long i=0 ; i<(long long)dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );
#else // !NEW_CODE_SPARSE_MATRIX
			for( int i=0 ; i<dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
		}

		Real beta = delta_new / delta_old;
#ifdef NEW_THREADS
		tp.parallel_for( 0 , dim , [&]( unsigned int , size_t i ){ d[i] = r[i] + (T)( d[i] * beta ); } );
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE_SPARSE_MATRIX
		for( long long i=0 ; i<(long long)dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
#else // !NEW_CODE_SPARSE_MATRIX
		for( int i=0 ; i<dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
	}
	FreePointer( Pb );
	FreePointer( r );
	FreePointer( d );
	FreePointer( q );
	FreePointer( temp );
	return ii;
}
