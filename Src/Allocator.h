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

#ifndef ALLOCATOR_INCLUDED
#define ALLOCATOR_INCLUDED
#include <vector>

#ifdef NEW_CODE
struct AllocatorState
{
	size_t index , remains;
	AllocatorState( void ) : index(0) , remains(0) {}
};
#else // !NEW_CODE
struct AllocatorState{ int index , remains; };
#endif // NEW_CODE
/** This templated class assists in memory allocation and is well suited for instances
  * when it is known that the sequence of memory allocations is performed in a stack-based
  * manner, so that memory allocated last is released first. It also preallocates memory
  * in chunks so that multiple requests for small chunks of memory do not require separate
  * system calls to the memory manager.
  * The allocator is templated off of the class of objects that we would like it to allocate,
  * ensuring that appropriate constructors and destructors are called as necessary.
  */
template< class T >
class SingleThreadedAllocator
{
#ifdef NEW_CODE
	size_t _blockSize;
	AllocatorState _state;
	std::vector< T* > _memory;
#else // !NEW_CODE
	int blockSize;
	int index , remains;
	std::vector< T* > memory;
#endif // NEW_CODE
public:
#ifdef NEW_CODE
	SingleThreadedAllocator( void ) : _blockSize(0) {}
#else // !NEW_CODE
	SingleThreadedAllocator( void ){ blockSize = index = remains = 0; }
#endif // NEW_CODE
	~SingleThreadedAllocator( void ){ reset(); }

	/** This method is the allocators destructor. It frees up any of the memory that
	  * it has allocated. */
	void reset( void )
	{
#ifdef NEW_CODE
		for( size_t i=0 ; i<_memory.size() ; i++ ) delete[] _memory[i];
		_memory.clear();
		_blockSize = 0;
		_state = AllocatorState();
#else // !NEW_CODE
		for( size_t i=0 ; i<memory.size() ; i++ ) delete[] memory[i];
		memory.clear();
		blockSize=index=remains=0;
#endif // NEW_CODE
	}
	/** This method returns the memory state of the allocator. */
	AllocatorState getState( void ) const
	{
#ifdef NEW_CODE
		return _state;
#else // !NEW_CODE
		AllocatorState s;
		s.index=index;
		s.remains=remains;
		return s;
#endif // NEW_CODE
	}


	/** This method rolls back the allocator so that it makes all of the memory previously
	  * allocated available for re-allocation. Note that it does it not call the constructor
	  * again, so after this method has been called, assumptions about the state of the values
	  * in memory are no longer valid. */
	void rollBack( void )
	{
#ifdef NEW_CODE
		if( _memory.size() )
		{
			for( size_t i=0 ; i<_memory.size() ; i++ ) for( size_t j=0 ; j<_blockSize ; j++ )
			{
				_memory[i][j].~T();
				new( &_memory[i][j] ) T();
			}
			_state = AllocatorState();
		}
#else // !NEW_CODE
		if( memory.size() )
		{
			for( size_t i=0 ; i<memory.size() ; i++ )
			{
				for( int j=0 ; j<blockSize ; j++ )
				{
					memory[i][j].~T();
					new(&memory[i][j]) T();
				}
			}
			index=0;
			remains=blockSize;
		}
#endif // NEW_CODE
	}
	/** This method rolls back the allocator to the previous memory state and makes all of the memory previously
	  * allocated available for re-allocation. Note that it does it not call the constructor
	  * again, so after this method has been called, assumptions about the state of the values
	  * in memory are no longer valid. */
	void rollBack( const AllocatorState& state )
	{
#ifdef NEW_CODE
		if( state.index<_state.index || ( state.index==_state.index && state.remains<_state.remains ) )
		{
			if( state.index<_state.index )
			{
				for( size_t j=state.remains ; j<_blockSize ; j++ )
				{
					_memory[ state.index ][j].~T();
					new( &_memory[ state.index ][j] ) T();
				}
				for( size_t i=state.index+1 ; i<_state.index-1 ; i++ ) for( size_t j=0 ; j<_blockSize ; j++ )
				{
					_memory[i][j].~T();
					new( &_memory[i][j] ) T();
				}
				for( size_t j=0 ; j<_state.remains ; j++ )
				{
					_memory[ _state.index ][j].~T();
					new( &_memory[ _state.index ][j] ) T();
				}
				_state = state;
			}
			else
			{
				for( size_t j=0 ; j<state.remains ; j++ )
				{
					_memory[ _state.index ][j].~T();
					new( &_memory[ _state.index ][j] ) T();
				}
				_state.remains = state.remains;
			}
		}
#else // !NEW_CODE
		if(state.index<index || (state.index==index && state.remains<remains)){
			if(state.index<index){
				for(int j=state.remains;j<blockSize;j++){
					memory[state.index][j].~T();
					new(&memory[state.index][j]) T();
				}
				for(int i=state.index+1;i<index-1;i++){
					for(int j=0;j<blockSize;j++){
						memory[i][j].~T();
						new(&memory[i][j]) T();
					}
				}
				for(int j=0;j<remains;j++){
					memory[index][j].~T();
					new(&memory[index][j]) T();
				}
				index=state.index;
				remains=state.remains;
			}
			else{
				for(int j=0;j<state.remains;j++){
					memory[index][j].~T();
					new(&memory[index][j]) T();
				}
				remains=state.remains;
			}
		}
#endif // NEW_CODE
	}

	/** This method initiallizes the constructor and the blockSize variable specifies the
	  * the number of objects that should be pre-allocated at a time. */
#ifdef NEW_CODE
	void set( size_t blockSize )
#else // !NEW_CODE
	void set( int blockSize )
#endif // NEW_CODE
	{
		reset();
#ifdef NEW_CODE
		_blockSize = blockSize;
		_state.index = -1;
		_state.remains = 0;
#else // !NEW_CODE
		this->blockSize = blockSize;
		index=-1;
		remains=0;
#endif // NEW_CODE
	}

	/** This method returns a pointer to an array of elements objects. If there is left over pre-allocated
	  * memory, this method simply returns a pointer to the next free piece of memory, otherwise it pre-allocates
	  * more memory. Note that if the number of objects requested is larger than the value blockSize with which
	  * the allocator was initialized, the request for memory will fail.
	  */
#ifdef NEW_CODE
	T* newElements( size_t elements=1 )
	{
		T* mem;
		if( !elements ) return NULL;
		if( elements>_blockSize ) ERROR_OUT( "elements bigger than block-size: " , elements , " > " , _blockSize );
		if( _state.remains<elements )
		{
			if( _state.index==_memory.size()-1 )
			{
				mem = new T[ _blockSize ];
				if( !mem ) ERROR_OUT( "Failed to allocate memory" );
				_memory.push_back( mem );
			}
			_state.index++;
			_state.remains = _blockSize;
		}
		mem = &( _memory[ _state.index ][ _blockSize-_state.remains ] );
		_state.remains -= elements;
		return mem;
	}
#else // !NEW_CODE
	T* newElements( int elements=1 )
	{
		T* mem;
		if( !elements ) return NULL;
		if( elements>blockSize ) ERROR_OUT( "elements bigger than block-size: " , elements , " > " , blockSize );
		if( remains<elements )
		{
			if( index==memory.size()-1 )
			{
				mem = new T[blockSize];
				if( !mem ) ERROR_OUT( "Failed to allocate memory" );
				memory.push_back( mem );
			}
			index++;
			remains=blockSize;
		}
		mem = &(memory[index][blockSize-remains]);
		remains -= elements;
		return mem;
	}
#endif // NEW_CODE
};
#ifdef NEW_CODE
#include <thread>
#endif // NEW_CODE
template< class T >
class Allocator
{
#ifdef TEST_ALLOCATOR_LOCK
public:
	bool enableThreadIDTest;
	std::thread::id rootThreadID;
protected:
#endif // TEST_ALLOCATOR_LOCK
	SingleThreadedAllocator< T >* _allocators;
#ifdef NEW_CODE
	size_t _maxThreads;
#else // !NEW_CODE
	int _maxThreads;
#endif // NEW_CODE
public:
	Allocator( void )
	{
#ifdef TEST_ALLOCATOR_LOCK
		enableThreadIDTest = true;
		rootThreadID = std::this_thread::get_id();
#endif // TEST_ALLOCATOR_LOCK
		_maxThreads = omp_get_max_threads();
		_allocators = new SingleThreadedAllocator< T >[_maxThreads];
	}
	~Allocator( void ){ delete[] _allocators; }

#ifdef NEW_CODE
	void set( size_t blockSize ){ for( size_t t=0 ; t<_maxThreads ; t++ ) _allocators[t].set( blockSize ); }
	T* newElements( size_t elements=1 )
	{
#ifdef TEST_ALLOCATOR_LOCK
		std::thread::id id = std::this_thread::get_id();
		if( id!=rootThreadID && enableThreadIDTest ) ERROR_OUT( "thread IDs don't match:  "  , rootThreadID , " <-> " , id );
#endif // TEST_ALLOCATOR_LOCK
		return _allocators[ omp_get_thread_num() ].newElements( elements );
	}
#else // !NEW_CODE
	void set( int blockSize ){ for( int t=0 ; t<_maxThreads ; t++ ) _allocators[t].set( blockSize ); }
	T* newElements( int elements=1 ){ return _allocators[ omp_get_thread_num() ].newElements( elements ); }
#endif // NEW_CODE
};

#endif // ALLOCATOR_INCLUDE
