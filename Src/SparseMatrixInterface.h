#ifndef SPARSE_MATRIX_INTERFACE_INCLUDED
#define SPARSE_MATRIX_INTERFACE_INCLUDED

#define NEW_CODE_SPARSE_MATRIX

#define FORCE_TWO_BYTE_ALIGNMENT 1
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
	template< typename Index >
	double squareASymmetricNorm( Index &idx1 , Index &idx2 ) const;
#else // !NEW_CODE_SPARSE_MATRIX
	double squareASymmetricNorm( int& idx1 , int& idx2 ) const;
#endif // NEW_CODE_SPARSE_MATRIX

#ifdef NEW_CODE_SPARSE_MATRIX
	template< class T2 > void multiply      (           ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const;
	template< class T2 > void multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const;
	template< class T2 > void multiply      (                Pointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const { multiply      (         ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
	template< class T2 > void multiplyScaled( T scale ,      Pointer( T2 ) In , Pointer( T2 ) Out , char multiplyFlag=0 ) const { multiplyScaled( scale , ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
#else // !NEW_CODE_SPARSE_MATRIX
	template< class T2 > void multiply      (           ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const;
	template< class T2 > void multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const;
	template< class T2 > void multiply      (                Pointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const { multiply      (         ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
	template< class T2 > void multiplyScaled( T scale ,      Pointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const { multiplyScaled( scale , ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
#endif // NEW_CODE_SPARSE_MATRIX

	void setDiagonal( Pointer( T ) diagonal ) const;
	void setDiagonalR( Pointer( T ) diagonal ) const;
	template< class T2 > void jacobiIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , ConstPointer( T2 ) in , Pointer( T2 ) out , bool dReciprocal ) const;
	template< class T2 > void gsIteration(                                                              ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const;
#ifdef NEW_CODE_SPARSE_MATRIX
	template< typename Index , class T2 > void gsIteration( const              std::vector< Index >  & multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x ,                bool dReciprocal ) const;
	template< typename Index , class T2 > void gsIteration( const std::vector< std::vector< Index > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const;
#else // !NEW_CODE_SPARSE_MATRIX
	template< class T2 > void gsIteration( const              std::vector< int >  & multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x ,                bool dReciprocal ) const;
	template< class T2 > void gsIteration( const std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const;
#endif // NEW_CODE_SPARSE_MATRIX
};

// Assuming that the SPDOperator class defines:
//		auto SPDOperator::()( ConstPointer( T ) , Pointer( T ) ) const
#ifdef NEW_CODE_SPARSE_MATRIX
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , size_t dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , const Preconditioner& P , size_t dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
#else // !NEW_CODE_SPARSE_MATRIX
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , const Preconditioner& P , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
#endif // NEW_CODE_SPARSE_MATRIX

#include "SparseMatrixInterface.inl"
#endif // SPARSE_MATRIX_INTERFACE_INCLUDED
