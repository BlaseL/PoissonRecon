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

#ifndef SPARSE_MATRIX_INTERFACE_INCLUDED
#define SPARSE_MATRIX_INTERFACE_INCLUDED

#define NEW_CODE_SPARSE_MATRIX

#define FORCE_TWO_BYTE_ALIGNMENT 1
#ifdef NEW_THREADS
#include "MyMiscellany.h"
#endif // NEW_THREADS
#include "Array.h"

#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(push)
#pragma pack(2)
#endif // FORCE_TWO_BYTE_ALIGNMENT
#ifdef NEW_CODE_SPARSE_MATRIX
template< class T , class IndexType >
#else // !NEW_CODE_SPARSE_MATRIX
#pragma message( "[WARNING] change me not to use default int" )
template< class T , class IndexType=int >
#endif // NEW_CODE_SPARSE_MATRIX
struct MatrixEntry
{
	MatrixEntry( void )             { N =-1 , Value = 0; }
	MatrixEntry( IndexType i )      { N = i , Value = 0; }
	MatrixEntry( IndexType n , T v ){ N = n , Value = v; }
	IndexType N;
	T Value;
};

#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(pop)
#endif // FORCE_TWO_BYTE_ALIGNMENT

enum
{
	MULTIPLY_ADD = 1 ,
	MULTIPLY_NEGATE = 2
};

//#pragma message( "[WARNING] make me templated off of IndexType as well" )
template< class T , class const_iterator > class SparseMatrixInterface
{
public:
	virtual const_iterator begin( size_t row ) const = 0;
	virtual const_iterator end  ( size_t row ) const = 0;
	virtual size_t rows   ( void )             const = 0;
	virtual size_t rowSize( size_t idx )       const = 0;

	size_t entries( void ) const;

	double squareNorm( void ) const;
	double squareASymmetricNorm( void ) const;
#ifdef NEW_CODE_SPARSE_MATRIX
	double squareASymmetricNorm( size_t &idx1 , size_t &idx2 ) const;
#else // !NEW_CODE_SPARSE_MATRIX
	double squareASymmetricNorm( int& idx1 , int& idx2 ) const;
#endif // NEW_CODE_SPARSE_MATRIX

#ifdef NEW_CODE_SPARSE_MATRIX
#ifdef NEW_THREADS
#if !defined( FORCE_PARALLEL )
	template< class T2 > void multiply      (           ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const
	{
		ThreadPool tp;
		return multiply( tp , In , Out , multiplyFlag );
	}
	template< class T2 > void multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const
	{
		ThreadPool tp;
		return multiplyScaled( tp , scale , In , Out , multiplyFlag );
	}
	template< class T2 > void multiply      (                Pointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const
	{
		ThreadPool tp;
		multiply( tp , In , Out , multiplyFlag );
	}
	template< class T2 > void multiplyScaled( T scale ,      Pointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const
	{
		ThreadPool tp;
		multiply( tp , scale , In , Out , multiplyFlag );
	}
#endif // !FORCE_PARALLEL

	template< class T2 > void multiply      ( ThreadPool &tp ,           ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const;
	template< class T2 > void multiplyScaled( ThreadPool &tp , T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const;
	template< class T2 > void multiply      (  ThreadPool &tp ,               Pointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const { multiply      ( tp,          ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
	template< class T2 > void multiplyScaled( ThreadPool &tp , T scale ,      Pointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const { multiplyScaled( tp,  scale , ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
#else // !NEW_THREADS
	template< class T2 > void multiply      (           ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const;
	template< class T2 > void multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const;
	template< class T2 > void multiply      (                Pointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const { multiply      (         ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
	template< class T2 > void multiplyScaled( T scale ,      Pointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const { multiplyScaled( scale , ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
#endif // NEW_THREADS
#else // !NEW_CODE_SPARSE_MATRIX
	template< class T2 > void multiply      (           ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const;
	template< class T2 > void multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const;
	template< class T2 > void multiply      (                Pointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const { multiply      (         ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
	template< class T2 > void multiplyScaled( T scale ,      Pointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const { multiplyScaled( scale , ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
#endif // NEW_CODE_SPARSE_MATRIX

#ifdef NEW_THREADS
#if !defined( FORCE_PARALLEL )
	void setDiagonal( Pointer( T ) diagonal ) const
	{
		ThreadPool tp;
		return setDiagonal( tp , diagonal );
	}
	void setDiagonalR( Pointer( T ) diagonal ) const
	{
		ThreadPool tp;
		return setDiagonalR( tp , diagonal );
	}
	template< class T2 > void jacobiIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , ConstPointer( T2 ) in , Pointer( T2 ) out , bool dReciprocal ) const
	{
		ThreadPool tp;
		return jacobiIteration( tp , diagonal , b , in , out , dReciprocal );
	}
	template< class T2 > void gsIteration( const              std::vector< size_t >  & multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x ,                bool dReciprocal ) const
	{
		ThreadPool tp;
		return gsIteration( tp , multiColorIndices , diagonal , b , x , dReciprocal );
	}
	template< class T2 > void gsIteration( const std::vector< std::vector< size_t > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const
	{
		ThreadPool tp;
		return gsIteration( tp , multiColorIndices , diagonal , b , x , forward , dReciprocal );
	}
#endif // !FORCE_PARALLEL

	void setDiagonal( ThreadPool &tp , Pointer( T ) diagonal ) const;
	void setDiagonalR( ThreadPool &tp , Pointer( T ) diagonal ) const;
	template< class T2 > void jacobiIteration( ThreadPool &tp , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , ConstPointer( T2 ) in , Pointer( T2 ) out , bool dReciprocal ) const;
	template< class T2 > void gsIteration( ThreadPool &tp , const              std::vector< size_t >  & multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x ,                bool dReciprocal ) const;
	template< class T2 > void gsIteration( ThreadPool &tp , const std::vector< std::vector< size_t > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const;
#else // !NEW_THREADS
	void setDiagonal( Pointer( T ) diagonal ) const;
	void setDiagonalR( Pointer( T ) diagonal ) const;
	template< class T2 > void jacobiIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , ConstPointer( T2 ) in , Pointer( T2 ) out , bool dReciprocal ) const;
#ifdef NEW_CODE_SPARSE_MATRIX
	template< class T2 > void gsIteration( const              std::vector< size_t >  & multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x ,                bool dReciprocal ) const;
	template< class T2 > void gsIteration( const std::vector< std::vector< size_t > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const;
#else // !NEW_CODE_SPARSE_MATRIX
	template< class T2 > void gsIteration( const              std::vector< int >  & multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x ,                bool dReciprocal ) const;
	template< class T2 > void gsIteration( const std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const;
#endif // NEW_CODE_SPARSE_MATRIX
#endif // NEW_THREADS
	template< class T2 > void gsIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const;
};

// Assuming that the SPDOperator class defines:
//		auto SPDOperator::()( ConstPointer( T ) , Pointer( T ) ) const
#ifdef NEW_CODE_SPARSE_MATRIX
#ifdef NEW_THREADS
#if !defined( FORCE_PARALLEL )
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > size_t SolveCG( const SPDFunctor& M , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot )
{
	ThreadPool tp;
	return SolveCG( tp , M , dim , b , iters ,  x , eps , Dot );
}
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > size_t SolveCG( const SPDFunctor& M , const Preconditioner& P , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot )
{
	ThreadPool tp;
	return SolveCG( tp , M , P , dim , b , iters , x , eps , Dot );
}
#endif // !FORCE_PARALLEL

template< class SPDFunctor , class T , typename Real , class TDotTFunctor > size_t SolveCG( ThreadPool &tp , const SPDFunctor& M , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > size_t SolveCG( ThreadPool &tp , const SPDFunctor& M , const Preconditioner& P , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
#else // !NEW_THREADS
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > size_t SolveCG( const SPDFunctor& M , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > size_t SolveCG( const SPDFunctor& M , const Preconditioner& P , size_t dim , ConstPointer( T ) b , size_t iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
#endif // NEW_THREADS
#else // !NEW_CODE_SPARSE_MATRIX
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , const Preconditioner& P , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
#endif // NEW_CODE_SPARSE_MATRIX

#include "SparseMatrixInterface.inl"
#endif // SPARSE_MATRIX_INTERFACE_INCLUDED
