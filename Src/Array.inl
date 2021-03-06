/*
Copyright (c) 2011, Michael Kazhdan and Ming Chuang
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
#define ARRAY_NEW_CODE
#undef ARRAY_NEW_CODE_2

#include <string.h>
#include <stdio.h>
#include <emmintrin.h>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#endif // _WIN32
#include <stddef.h>
#ifdef ARRAY_NEW_CODE
#include <type_traits>
#include <cstddef>
#endif // ARRAY_NEW_CODE


#ifdef ARRAY_NEW_CODE
#else // !ARRAY_NEW_CODE

inline bool isfinitef( float fp ){ float f=fp; return ((*(unsigned *)&f)&0x7f800000)!=0x7f800000; }


template< class C >        bool IsValid( const C& c );
#if _DEBUG
template< >         inline bool IsValid< float >( const float& f ) { return isfinitef( f ) &&  ( f==0.f || abs(f)>1e-31f ); }
#else // !_DEBUG
template< >         inline bool IsValid< float >( const float& f ) { return isfinitef( f ); }
#endif // _DEBUG
template< >         inline bool IsValid< __m128 >( const __m128& m )
{
	const __m128* addr = &m;
	if( size_t(addr) & 15 ) return false;
	else                    return true;
}
template< class C > inline bool IsValid( const C& c ){ return true; }
#endif // ARRAY_NEW_CODE


template< class C >
class Array
{
#ifdef ARRAY_NEW_CODE_2
	template< class _C > friend class ConstArray;
#endif // ARRAY_NEW_CODE_2
	template< class D > friend class Array;
#ifdef ARRAY_NEW_CODE
	void _assertBounds( std::ptrdiff_t idx ) const
#else // !ARRAY_NEW_CODE
	void _assertBounds( long long idx ) const
#endif // ARRAY_NEW_CODE
	{
#ifdef ARRAY_NEW_CODE_2
		if( idx<min || idx>=max )
		{
#ifdef NEW_CODE
			StackTracer::Trace();
#endif // NEW_CODE
			ERROR_OUT( "Array index out-of-bounds: " , min , " <= " , idx , " < " , max , "\nFile: " , _fileName , "; Line: " , _line , "; Function: " , _functionName );
		}
#else // !ARRAY_NEW_CODE_2
		if( idx<min || idx>=max )
		{
#ifdef NEW_CODE
			StackTracer::Trace();
#endif // NEW_CODE
			ERROR_OUT( "Array index out-of-bounds: " , min , " <= " , idx , " < " , max );
		}
#endif // ARRAY_NEW_CODE_2
	}
protected:
	C *data , *_data;
#ifdef ARRAY_NEW_CODE
	std::ptrdiff_t min , max;
#else // !ARRAY_NEW_CODE
	long long min , max;
#endif // ARRAY_NEW_CODE
#ifdef ARRAY_NEW_CODE_2
	std::string _fileName , _functionName;
	std::string _line;
#endif // ARRAY_NEW_CODE_2

public:
#ifdef ARRAY_NEW_CODE
	std::ptrdiff_t minimum( void ) const { return min; }
	std::ptrdiff_t maximum( void ) const { return max; }
#else // !ARRAY_NEW_CODE
	long long minimum( void ) const { return min; }
	long long maximum( void ) const { return max; }
#endif // ARRAY_NEW_CODE

#ifdef ARRAY_NEW_CODE_2
	static inline Array New( size_t size , std::string fileName , int line , std::string functionName )
#else // !ARRAY_NEW_CODE_2
	static inline Array New( size_t size , const char* name=NULL )
#endif // ARRAY_NEW_CODE_2
	{
		Array a;
		a._data = a.data = new C[size];
		a.min = 0;
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Casting unsigned to signed" )
#endif // SHOW_WARNINGS
#ifdef ARRAY_NEW_CODE
		a.max = (std::ptrdiff_t)size;
#else // !ARRAY_NEW_CODE
		a.max = ( long long ) size;
#endif // ARRAY_NEW_CODE
#ifdef ARRAY_NEW_CODE_2
		a._fileName = fileName , a._line = line , a._functionName = functionName;
#endif // ARRAY_NEW_CODE_2
		return a;
	}
#ifdef ARRAY_NEW_CODE_2
	static inline Array Alloc( size_t size , bool clear , std::string fileName , int line , std::string functionName )
#else // !ARRAY_NEW_CODE_2
	static inline Array Alloc( size_t size , bool clear , const char* name=NULL )
#endif // ARRAY_NEW_CODE_2
	{
		Array a;
		a._data = a.data = ( C* ) malloc( size * sizeof( C ) );
		if( clear ) memset( a.data ,  0 , size * sizeof( C ) );
		a.min = 0;
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Casting unsigned to signed" )
#endif // SHOW_WARNINGS
#ifdef ARRAY_NEW_CODE
		a.max = (std::ptrdiff_t)size;
#else // !ARRAY_NEW_CODE
		a.max = ( long long ) size;
#endif // ARRAY_NEW_CODE
#ifdef ARRAY_NEW_CODE_2
		a._fileName = fileName , a._line = line , a._functionName = functionName;
#endif // ARRAY_NEW_CODE_2
		return a;
	}

#ifdef ARRAY_NEW_CODE_2
	static inline Array AlignedAlloc( size_t size , size_t alignment , bool clear , std::string fileName , int line , std::string functionName )
#else // !ARRAY_NEW_CODE_2
	static inline Array AlignedAlloc( size_t size , size_t alignment , bool clear , const char* name=NULL )
#endif // ARRAY_NEW_CODE_2
	{
		Array a;
		a.data = ( C* ) aligned_malloc( sizeof(C) * size , alignment );
		a._data = ( C* )( ( ( void** )a.data )[-1] );
		if( clear ) memset( a.data ,  0 , size * sizeof( C ) );
		a.min = 0;
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Casting unsigned to signed" )
#endif // SHOW_WARNINGS
#ifdef ARRAY_NEW_CODE
		a.max = (std::ptrdiff_t)size;
#else // !ARRAY_NEW_CODE
		a.max = ( long long ) size;
#endif // ARRAY_NEW_CODE
#ifdef ARRAY_NEW_CODE_2
		a._fileName = fileName , a._line = line , a._functionName = functionName;
#endif // ARRAY_NEW_CODE_2
		return a;
	}

#ifdef ARRAY_NEW_CODE_2
	static inline Array ReAlloc( Array& a , size_t size , bool clear , std::string fileName , int line , std::string functionName )
#else // !ARRAY_NEW_CODE_2
	static inline Array ReAlloc( Array& a , size_t size , bool clear , const char* name=NULL )
#endif // ARRAY_NEW_CODE_2
	{
		Array _a;
		_a._data = _a.data = ( C* ) realloc( a.data , size * sizeof( C ) );
		if( clear ) memset( _a.data ,  0 , size * sizeof( C ) );
		a._data = NULL;
		_a.min = 0;
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Casting unsigned to signed" )
#endif // SHOW_WARNINGS
#ifdef ARRAY_NEW_CODE
		_a.max = (std::ptrdiff_t)size;
#else // !ARRAY_NEW_CODE
		_a.max = ( long long ) size;
#endif // ARRAY_NEW_CODE
#ifdef ARRAY_NEW_CODE_2
		a._fileName = fileName , a._line = line , a._functionName = functionName;
#endif // ARRAY_NEW_CODE_2
		return _a;
	}

	Array( void )
	{
		data = _data = NULL;
		min = max = 0;
	}
	template< class D >
	Array( Array< D >& a )
	{
		_data = NULL;
		if( !a )
		{
			data =  NULL;
			min = max = 0;
		}
		else
		{
			// [WARNING] Changing szC and szD to size_t causes some really strange behavior.
#ifdef ARRAY_NEW_CODE
			std::ptrdiff_t szC = (std::ptrdiff_t)sizeof( C );
			std::ptrdiff_t szD = (std::ptrdiff_t)sizeof( D );
#else // !ARRAY_NEW_CODE
			long long szC = sizeof( C );
			long long szD = sizeof( D );
#endif // ARRAY_NEW_CODE
			data = (C*)a.data;
			min = ( a.minimum() * szD ) / szC;
			max = ( a.maximum() * szD ) / szC;
			if( min*szC!=a.minimum()*szD || max*szC!=a.maximum()*szD ) ERROR_OUT( "Could not convert array [ " , a.minimum() , " , " , a.maximum() , " ] * " , szD , " => [ " , min , " , " , max , " ] * " , szC );
#ifdef ARRAY_NEW_CODE_2
			_fileName = a._fileName;
			_line = a._line;
			_functionName = a._functionName;
#endif // ARRAY_NEW_CODE_2
		}
	}
#ifdef ARRAY_NEW_CODE
	static Array FromPointer( C* data , std::ptrdiff_t max )
#else // !ARRAY_NEW_CODE
	static Array FromPointer( C* data , long long max )
#endif // ARRAY_NEW_CODE
	{
		Array a;
		a._data = NULL;
		a.data = data;
		a.min = 0;
		a.max = max;
		return a;
	}
#ifdef ARRAY_NEW_CODE
	static Array FromPointer( C* data , std::ptrdiff_t min , std::ptrdiff_t max )
#else // !ARRAY_NEW_CODE
	static Array FromPointer( C* data , long long min , long long max )
#endif // ARRAY_NEW_CODE
	{
		Array a;
		a._data = NULL;
		a.data = data;
		a.min = min;
		a.max = max;
		return a;
	}
	inline bool operator == ( const Array< C >& a ) const { return data==a.data; }
	inline bool operator != ( const Array< C >& a ) const { return data!=a.data; }
	inline bool operator == ( const C* c ) const { return data==c; }
	inline bool operator != ( const C* c ) const { return data!=c; }
	inline C* operator -> ( void )
	{
		_assertBounds( 0 );
		return data;
	}
	inline const C* operator -> ( void ) const
	{
		_assertBounds( 0 );
		return data;
	}
#ifdef NEW_CODE
	inline C &operator * ( void )
	{
		_assertBounds( 0 );
		return data[0];
	}
	inline const C &operator * ( void ) const
	{
		_assertBounds( 0 );
		return data[0];
	}
#endif // NEW_CODE
#ifdef ARRAY_NEW_CODE
	template< typename Int >
	inline C& operator[]( Int idx )
	{
		static_assert( std::is_integral< Int >::value , "Integral type expected" );
		_assertBounds( idx );
		return data[idx];
	}

	template< typename Int >
	inline const C& operator[]( Int idx ) const
	{
		static_assert( std::is_integral< Int >::value , "Integral type expected" );
		_assertBounds( idx );
		return data[idx];
	}

	template< typename Int >
	inline Array operator + ( Int idx ) const
	{
		static_assert( std::is_integral< Int >::value , "Integral type expected" );
		Array a;
		a._data = _data;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	template< typename Int >
	inline Array& operator += ( Int idx  )
	{
		static_assert( std::is_integral< Int >::value , "Integral type expected" );
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	template< typename Int > Array  operator -  ( Int idx ) const { return (*this) +  (-idx); }
	template< typename Int > Array& operator -= ( Int idx )       { return (*this) += (-idx); }

#else // !ARRAY_NEW_CODE
	inline C& operator[]( long long idx )
	{
		_assertBounds( idx );
		return data[idx];
	}
	inline const C& operator[]( long long idx ) const
	{
		_assertBounds( idx );
		return data[idx];
	}
	inline Array operator + ( int idx ) const
	{
		Array a;
		a._data = _data;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline Array operator + ( long long idx ) const
	{
		Array a;
		a._data = _data;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline Array operator + ( unsigned int idx ) const
	{
		Array a;
		a._data = _data;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline Array operator + ( unsigned long long idx ) const
	{
		Array a;
		a._data = _data;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline Array& operator += ( int idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	inline Array& operator += ( long long idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	inline Array& operator += ( unsigned int idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	inline Array& operator += ( unsigned long long idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	Array  operator -  ( int idx ) const { return (*this) +  (-idx); }
	Array  operator -  ( long long idx ) const { return (*this) +  (-idx); }
	Array  operator -  ( unsigned int idx ) const { return (*this) +  (-idx); }
	Array  operator -  ( unsigned long long idx ) const { return (*this) +  (-idx); }
	Array& operator -= ( int idx )    { return (*this) += (-idx); }
	Array& operator -= ( long long idx )    { return (*this) += (-idx); }
	Array& operator -= ( unsigned int idx )    { return (*this) += (-idx); }
	Array& operator -= ( unsigned long long idx )    { return (*this) += (-idx); }
	Array& operator -- ( void ) { return (*this) -= 1; }
#endif // ARRAY_NEW_CODE
	inline Array& operator ++ ( void  ) { return (*this) += 1; }
	inline Array operator++( int ){ Array< C > temp = (*this) ; (*this) +=1 ; return temp; }
	inline Array operator--( int ){ Array< C > temp = (*this) ; (*this) -=1 ; return temp; }
#ifdef ARRAY_NEW_CODE
	std::ptrdiff_t operator - ( const Array& a ) const { return data - a.data; }
#else // !ARRAY_NEW_CODE
	long long operator - ( const Array& a ) const { return ( long long )( data - a.data ); }
#endif // ARRAY_NEW_CODE

	void Free( void )
	{
		if( _data )
		{
			free( _data );
		}
		(*this) = Array( );
	}
	void Delete( void )
	{
		if( _data )
		{
			delete[] _data;
		}
		(*this) = Array( );
	}
	C* pointer( void ){ return data; }
	const C* pointer( void ) const { return data; }
	bool operator !( void ) const { return data==NULL; }
	operator bool( ) const { return data!=NULL; }
};

template< class C >
class ConstArray
{
	template< class D > friend class ConstArray;
#ifdef ARRAY_NEW_CODE
	void _assertBounds( std::ptrdiff_t idx ) const
#else // !ARRAY_NEW_CODE
	void _assertBounds( long long idx ) const
#endif // ARRAY_NEW_CODE
	{
#ifdef ARRAY_NEW_CODE_2
		if( idx<min || idx>=max ) ERROR_OUT( "ConstArray index out-of-bounds: " , min , " <= " , idx , " < " , max , "\nFile: " , _fileName , "; Line: " , _line , "; Function: " , _functionName );
#else // !ARRAY_NEW_CODE_2
		if( idx<min || idx>=max ) ERROR_OUT( "ConstArray index out-of-bounds: " , min , " <= " , idx , " < " , max );
#endif // ARRAY_NEW_CODE_2
	}
protected:
	const C *data;
#ifdef ARRAY_NEW_CODE
	std::ptrdiff_t min , max;
#else // !ARRAY_NEW_CODE
	long long min , max;
#endif // ARRAY_NEW_CODE
#ifdef ARRAY_NEW_CODE_2
	std::string _fileName , _functionName;
	std::string _line;
#endif // ARRAY_NEW_CODE_2
public:
#ifdef ARRAY_NEW_CODE
	std::ptrdiff_t minimum( void ) const { return min; }
	std::ptrdiff_t maximum( void ) const { return max; }
#else // !ARRAY_NEW_CODE
	long long minimum( void ) const { return min; }
	long long maximum( void ) const { return max; }
#endif // ARRAY_NEW_CODE

	inline ConstArray( void )
	{
		data = NULL;
		min = max = 0;
	}
	inline ConstArray( const Array< C >& a )
	{
		// [WARNING] Changing szC and szD to size_t causes some really strange behavior.
		data = ( const C* )a.pointer( );
		min = a.minimum();
		max = a.maximum();
#ifdef ARRAY_NEW_CODE_2
		_fileName = a._fileName;
		_line = a._line;
		_functionName = a._functionName;
#endif // ARRAY_NEW_CODE_2
	}
	template< class D >
	inline ConstArray( const Array< D >& a )
	{
		// [WARNING] Changing szC and szD to size_t causes some really strange behavior.
#ifdef ARRAY_NEW_CODE
		std::ptrdiff_t szC = (std::ptrdiff_t)sizeof( C );
		std::ptrdiff_t szD = (std::ptrdiff_t)sizeof( D );
#else // !ARRAY_NEW_CODE
		long long szC = ( long long ) sizeof( C );
		long long szD = ( long long ) sizeof( D );
#endif // ARRAY_NEW_CODE
		data = ( const C* )a.pointer( );
		min = ( a.minimum() * szD ) / szC;
		max = ( a.maximum() * szD ) / szC;
		if( min*szC!=a.minimum()*szD || max*szC!=a.maximum()*szD ) ERROR_OUT( "Could not convert const array [ " , a.minimum() , " , " , a.maximum() , " ] * " , szD , " => [ " , min , " , " , max , " ] * " , szC );
#ifdef ARRAY_NEW_CODE_2
		_fileName = a._fileName;
		_line = a._line;
		_functionName = a._functionName;
#endif // ARRAY_NEW_CODE_2
	}
	template< class D >
	inline ConstArray( const ConstArray< D >& a )
	{
		// [WARNING] Chaning szC and szD to size_t causes some really strange behavior.
#ifdef ARRAY_NEW_CODE
		std::ptrdiff_t szC = (std::ptrdiff_t)sizeof( C );
		std::ptrdiff_t szD = (std::ptrdiff_t)sizeof( D );
#else // !ARRAY_NEW_CODE
		long long szC = sizeof( C );
		long long szD = sizeof( D );
#endif // ARRAY_NEW_CODE
		data = ( const C*)a.pointer( );
		min = ( a.minimum() * szD ) / szC;
		max = ( a.maximum() * szD ) / szC;
		if( min*szC!=a.minimum()*szD || max*szC!=a.maximum()*szD ) ERROR_OUT( "Could not convert array [ " , a.minimum() , " , " , a.maximum() , " ] * " , szD , " => [ " , min , " , " , max , " ] * " , szC );
#ifdef ARRAY_NEW_CODE_2
		_fileName = a._fileName;
		_line = a._line;
		_functionName = a._functionName;
#endif // ARRAY_NEW_CODE_2
	}
#ifdef ARRAY_NEW_CODE
	static ConstArray FromPointer( const C* data , std::ptrdiff_t max )
#else // !ARRAY_NEW_CODE
	static ConstArray FromPointer( const C* data , long long max )
#endif // ARRAY_NEW_CODE
	{
		ConstArray a;
		a.data = data;
		a.min = 0;
		a.max = max;
		return a;
	}
#ifdef ARRAY_NEW_CODE
	static ConstArray FromPointer( const C* data , std::ptrdiff_t min , std::ptrdiff_t max )
#else // !ARRAY_NEW_CODE
	static ConstArray FromPointer( const C* data , long long min , long long max )
#endif // ARRAY_NEW_CODE
	{
		ConstArray a;
		a.data = data;
		a.min = min;
		a.max = max;
		return a;
	}

	inline bool operator == ( const ConstArray< C >& a ) const { return data==a.data; }
	inline bool operator != ( const ConstArray< C >& a ) const { return data!=a.data; }
	inline bool operator == ( const C* c ) const { return data==c; }
	inline bool operator != ( const C* c ) const { return data!=c; }
	inline const C* operator -> ( void )
	{
		_assertBounds( 0 );
		return data;
	}
#ifdef NEW_CODE
	inline const C &operator * ( void )
	{
		_assertBounds( 0 );
		return data[0];
	}
#endif // NEW_CODE
#ifdef ARRAY_NEW_CODE
	template< typename Int >
	inline const C& operator[]( Int idx ) const
	{
		static_assert( std::is_integral< Int >::value , "Integral type expected" );
		_assertBounds( idx );
		return data[idx];
	}
	template< typename Int >
	inline ConstArray operator + ( Int idx ) const
	{
		static_assert( std::is_integral< Int >::value , "Integral type expected" );
		ConstArray a;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	template< typename  Int >
	inline ConstArray& operator += ( Int idx  )
	{
		static_assert( std::is_integral< Int >::value , "Integral type expected" );
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	template< typename Int > ConstArray  operator -  ( Int idx ) const { return (*this) +  (-idx); }
	template< typename Int > ConstArray& operator -= ( Int idx )       { return (*this) += (-idx); }
#else // !ARRAY_NEW_CODE
	inline const C& operator[]( long long idx ) const
	{
		_assertBounds( idx );
		return data[idx];
	}
	inline ConstArray operator + ( int idx ) const
	{
		ConstArray a;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline ConstArray operator + ( long long idx ) const
	{
		ConstArray a;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline ConstArray operator + ( unsigned int idx ) const
	{
		ConstArray a;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline ConstArray operator + ( unsigned long long idx ) const
	{
		ConstArray a;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline ConstArray& operator += ( int idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	inline ConstArray& operator += ( long long idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	inline ConstArray& operator += ( unsigned int idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	inline ConstArray& operator += ( unsigned long long idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	ConstArray  operator -  ( int idx ) const { return (*this) +  (-idx); }
	ConstArray  operator -  ( long long idx ) const { return (*this) +  (-idx); }
	ConstArray  operator -  ( unsigned int idx ) const { return (*this) +  (-idx); }
	ConstArray  operator -  ( unsigned long long idx ) const { return (*this) +  (-idx); }
	ConstArray& operator -= ( int idx )    { return (*this) += (-idx); }
	ConstArray& operator -= ( long long idx )    { return (*this) += (-idx); }
	ConstArray& operator -= ( unsigned int idx )    { return (*this) += (-idx); }
	ConstArray& operator -= ( unsigned long long idx )    { return (*this) += (-idx); }
	ConstArray& operator -- ( void ) { return (*this) -= 1; }
#endif // ARRAY_NEW_CODE
	inline ConstArray& operator ++ ( void ) { return (*this) += 1; }
	inline ConstArray operator++( int ){ ConstArray< C > temp = (*this) ; (*this) +=1 ; return temp; }
	inline ConstArray operator--( int ){ ConstArray< C > temp = (*this) ; (*this) -=1 ; return temp; }
#ifdef ARRAY_NEW_CODE
	std::ptrdiff_t operator - ( const ConstArray& a ) const { return data - a.data; }
	std::ptrdiff_t operator - ( const Array< C >& a ) const { return data - a.pointer(); }
#else // !ARRAY_NEW_CODE
	long long operator - ( const ConstArray& a ) const { return ( long long )( data - a.data ); }
	long long operator - ( const Array< C >& a ) const { return ( long long )( data - a.pointer() ); }
#endif // ARRAY_NEW_CODE

	const C* pointer( void ) const { return data; }
	bool operator !( void ) { return data==NULL; }
	operator bool( ) { return data!=NULL; }
};

template< class C >
Array< C > memcpy( Array< C > destination , const void* source , size_t size )
{
	if( size>destination.maximum()*sizeof(C) ) ERROR_OUT( "Size of copy exceeds destination maximum: " , size , " > " , destination.maximum()*sizeof( C ) );
	if( size ) memcpy( &destination[0] , source , size );
	return destination;
}
template< class C , class D >
Array< C > memcpy( Array< C > destination , Array< D > source , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) ) ERROR_OUT( "Size of copy exceeds destination maximum: " , size , " > " , destination.maximum()*sizeof( C ) );
	if( size>source.maximum()*sizeof( D ) ) ERROR_OUT( "Size of copy exceeds source maximum: " , size , " > " , source.maximum()*sizeof( D ) );
	if( size ) memcpy( &destination[0] , &source[0] , size );
	return destination;
}
template< class C , class D >
Array< C > memcpy( Array< C > destination , ConstArray< D > source , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) ) ERROR_OUT( "Size of copy exceeds destination maximum: " , size , " > " , destination.maximum()*sizeof( C ) );
	if( size>source.maximum()*sizeof( D ) ) ERROR_OUT( "Size of copy exceeds source maximum: " , size , " > " , source.maximum()*sizeof( D ) );
	if( size ) memcpy( &destination[0] , &source[0] , size );
	return destination;
}
template< class D >
void* memcpy( void* destination , Array< D > source , size_t size )
{
	if( size>source.maximum()*sizeof( D ) ) ERROR_OUT( "Size of copy exceeds source maximum: " , size , " > " , source.maximum()*sizeof( D ) );
	if( size ) memcpy( destination , &source[0] , size );
	return destination;
}
template< class D >
void* memcpy( void* destination , ConstArray< D > source , size_t size )
{
	if( size>source.maximum()*sizeof( D ) ) ERROR_OUT( "Size of copy exceeds source maximum: " , size , " > " , source.maximum()*sizeof( D ) );
	if( size ) memcpy( destination , &source[0] , size );
	return destination;
}
template< class C >
Array< C > memset( Array< C > destination , int value , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) ) ERROR_OUT( "Size of set exceeds destination maximum: " , size , " > " , destination.maximum()*sizeof( C ) );
	if( size ) memset( &destination[0] , value , size );
	return destination;
}

template< class C >
size_t fread( Array< C > destination , size_t eSize , size_t count , FILE* fp )
{
	if( count*eSize>destination.maximum()*sizeof( C ) ) ERROR_OUT( "Size of read exceeds source maximum: " , count*eSize , " > " , destination.maximum()*sizeof( C ) );
	return fread( &destination[0] , eSize , count , fp );
}
template< class C >
size_t fwrite( Array< C > source , size_t eSize , size_t count , FILE* fp )
{
	if( count*eSize>source.maximum()*sizeof( C ) ) ERROR_OUT( "Size of write exceeds source maximum: " , count*eSize , " > " , source.maximum()*sizeof( C ) );
	return fwrite( &source[0] , eSize , count , fp );
}
template< class C >
size_t fwrite( ConstArray< C > source , size_t eSize , size_t count , FILE* fp )
{
	if( count*eSize>source.maximum()*sizeof( C ) ) ERROR_OUT( "Size of write exceeds source maximum: " , count*eSize , " > " , source.maximum()*sizeof( C ) );
	return fwrite( &source[0] , eSize , count , fp );
}
template< class C >
void qsort( Array< C > base , size_t numElements , size_t elementSize , int (*compareFunction)( const void* , const void* ) )
{
	if( sizeof(C)!=elementSize ) ERROR_OUT( "Element sizes differ: " , sizeof(C) , " != " , elementSize );
	if( base.minimum()>0 || base.maximum()<numElements ) ERROR_OUT( "Array access out of bounds: " , base.minimum() , " <= 0 <= " , base.maximum() , " <= " , numElements );
	qsort( base.pointer() , numElements , elementSize , compareFunction );
}
