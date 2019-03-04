/*
Copyright (c) 2013, Michael Kazhdan
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
*/
#define NEW_CODE

#include "PreProcessor.h"

#if defined( _WIN32 ) || defined( _WIN64 )
#define WINDOWS
#else
#undef WINDOWS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#ifdef WINDOWS
#include <windows.h>
#endif // WINDOWS
#include <atomic>
#include <condition_variable>
#include <chrono>
#include <functional>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "Array.h"

enum
{
	ATOMIC_NONE ,
	ATOMIC_ATOMIC ,
	ATOMIC_ATOMIC_FLOAT ,
	ATOMIC_ATOMIC_DOUBLE ,
	ATOMIC_MUTEX ,
	ATOMIC_REDUCTION ,
	ATOMIC_COUNT
};

const char * AtomicNames[] = { "none" , "atomic" , "atomic float" , "atomic double" , "mutex" , "reduction" };

cmdLineParameter< int > InnerIterations( "inIters" ) , OuterIterations( "outIters" , 1 );
cmdLineParameter< int > AtomicType( "atomic" , ATOMIC_NONE );
cmdLineParameter< int >
	ParallelType( "parallel" , (int)ThreadPool::OPEN_MP ) ,
	ScheduleType( "schedule" , (int)ThreadPool::DefaultSchedule ) ,
	ThreadChunkSize( "tChunkSize" , (int)ThreadPool::DefaultChunkSize ) ,
	Threads( "threads" , (int)std::thread::hardware_concurrency() );


cmdLineReadable* params[] = { &InnerIterations , &OuterIterations , &Threads , &ParallelType , &ScheduleType , &ThreadChunkSize ,  &AtomicType , NULL };

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <inner iterations count>\n" , InnerIterations.name );
	printf( "\t[--%s <outer iterations count>=%d]\n" , OuterIterations.name , OuterIterations.value );
	printf( "\t[--%s <threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <num threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <parallel type>=%d]\n" , ParallelType.name , ParallelType.value );
	for( size_t i=0 ; i<ThreadPool::ParallelNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ParallelNames[i].c_str() );
	printf( "\t[--%s <schedue type>=%d]\n" , ScheduleType.name , ScheduleType.value );
	for( size_t i=0 ; i<ThreadPool::ScheduleNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ScheduleNames[i].c_str() );
	printf( "\t[--%s <thread chunk size>=%d]\n" , ThreadChunkSize.name , ThreadChunkSize.value );
	printf( "\t[--%s <synchronization type>=%d]\n" , AtomicType.name , AtomicType.value );
	for( int i=0 ; i<ATOMIC_COUNT ;  i++ ) printf( "\t\t%d] %s\n" , i , AtomicNames[i] );
}

void AddAtomic( float &dest , float value )
{
	float current = dest;
	float sum = current+value;
#ifdef WINDOWS
	long &_current = *(long *)&current;
	long &_sum = *(long *)&sum;
	while( InterlockedCompareExchange( (long*)&dest , _sum , _current )!=_current ) current = dest , sum = dest+value;
#else // !WINDOWS
	uint32_t &_current = *(uint32_t *)&current;
	uint32_t &_sum = *(uint32_t *)&sum;
	while( __sync_val_compare_and_swap( (uint32_t *)&dest , _current , _sum )!=_current ) current = dest , sum = dest+value;
#endif // WINDOWS
};

void AddAtomic( double &dest , double value )
{
	double current = dest;
	double sum = current+value;
#ifdef WINDOWS
	__int64 &_current = *(__int64 *)&current;
	__int64 &_sum = *(__int64 *)&sum;
	while( InterlockedCompareExchange64( (__int64*)&dest , _sum , _current )!=_current ) current = dest , sum = dest+value;
#else // !WINDOWS
	uint64_t &_current = *(uint64_t *)&current;
	uint64_t &_sum = *(uint64_t *)&sum;
	while( __sync_val_compare_and_swap( (uint64_t*)&dest , _current , _sum )!=_current ) current = dest , sum = dest+value;
#endif // WINDOWS
};

double RunMutex( size_t outIters , size_t inIters )
{
	auto Task = []( size_t i , size_t &sum )
	{
		for( size_t j=0 ; j<i ; j++ )
		{
			static std::mutex m;
			std::lock_guard< std::mutex > lock(m);
			sum += j;
		}
	};

	double t = Time();
	{
		size_t sum = 0;
		for( size_t i=0 ; i<outIters ; i++ ) ThreadPool::Parallel_for( 0 , inIters , [&]( unsigned int , size_t ii ){ Task( ii , sum ); } );
		printf( "\t\tMutex:             %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	return Time()-t;
}

double RunAtomicDouble( size_t outIters , size_t inIters )
{
	auto Task = []( size_t i , double &sum )
	{
		for( size_t j=0 ; j<i ; j++ ) AddAtomic( sum , (double)j );
	};

	double t = Time();
	{
		double sum = 0;
		for( size_t i=0 ; i<outIters ; i++ ) ThreadPool::Parallel_for( 0 , inIters , [&]( unsigned int , size_t ii ){ Task( ii , sum ); } );
		printf( "\t\tAtomic add double: %.2f(s)\t%g\n" , Time()-t , sum );
	}
	return Time()-t;
}

double RunAtomicFloat( size_t outIters , size_t inIters )
{
	auto Task = []( size_t i , float &sum )
	{
		for( size_t j=0 ; j<i ; j++ ) AddAtomic( sum , (float)j );
	};

	double t = Time();
	{
		float sum = 0;
		for( size_t i=0 ; i<outIters ; i++ ) ThreadPool::Parallel_for( 0 , inIters , [&Task,&sum]( unsigned int , size_t ii ){ Task( ii , sum ); } );
		printf( "\t\tAtomic add float:  %.2f(s)\t%g\n" , Time()-t , sum );
	}
	return Time()-t;
}

double RunAtomic( size_t outIters , size_t inIters )
{
	auto Task = []( size_t i , std::atomic< size_t > &sum )
	{
		for( size_t j=0 ; j<i ; j++ ) sum += j;
	};

	double t = Time();
	{
		double t = Time();
		std::atomic< size_t > sum;
		sum.store(0);
		for( size_t i=0 ; i<outIters ; i++ ) ThreadPool::Parallel_for( 0 , inIters , [&]( unsigned int , size_t ii ){ Task( ii , sum ); } );
		printf( "\t\tAtomic:            %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	return Time()-t;
}

double Run( size_t outIters , size_t inIters )
{
	auto Task = []( size_t i , size_t &sum ){ for( size_t j=0 ; j<i ; j++ ) sum += j; };

	double t = Time();
	{
		size_t sum = 0;
		for( size_t i=0 ; i<outIters ; i++ ) ThreadPool::Parallel_for( 0 , inIters , [&]( unsigned int , size_t ii ){ Task( ii , sum ); } );
		printf( "\t\tNone:              %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	return Time()-t;
}

void foo( void )
{
	Pointer( int ) foo = NewPointer< int >( 10 );
	printf( "%d\n" , foo[10] );
}
int main( int argc , char* argv[] )
{
foo();
	cmdLineParse( argc-1 , &argv[1] , params );

	if( !InnerIterations.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	size_t outIters = OuterIterations.value;
	size_t inIters = InnerIterations.value;

	printf( "Iterations: %llu x %llu -> %llu\n" , (unsigned long long)outIters , (unsigned long long)inIters , (unsigned long long)outIters*inIters );
	printf( "Threads: %d\n" , ThreadPool::NumThreads() );


	ThreadPool::DefaultChunkSize = ThreadChunkSize.value;
	for( int i=0 ; i<ThreadPool::ParallelNames.size() ; i++ )
	{
		ThreadPool::Init( (ThreadPool::ParallelType)i , Threads.value );
		for( int j=0 ; j<ThreadPool::ScheduleNames.size() ; j++ )
		{
			std::cout << "Parallel / Schedule: " << ThreadPool::ParallelNames[i] << " / " << ThreadPool::ScheduleNames[j] << std::endl;
			ThreadPool::DefaultSchedule = (ThreadPool::ScheduleType)j;

			if( AtomicType.value<0 )
			{
				Run            ( outIters , inIters );
				RunAtomic      ( outIters , inIters );
				RunAtomicFloat ( outIters , inIters );
				RunAtomicDouble( outIters , inIters );
				RunMutex       ( outIters , inIters );
			}
			else
			{
				switch( AtomicType.value )
				{
				case ATOMIC_NONE:          Run            ( outIters , inIters ) ; break;
				case ATOMIC_ATOMIC:        RunAtomic      ( outIters , inIters ) ; break;
				case ATOMIC_ATOMIC_FLOAT:  RunAtomicFloat ( outIters , inIters ) ; break;
				case ATOMIC_ATOMIC_DOUBLE: RunAtomicDouble( outIters , inIters ) ; break;
				case ATOMIC_MUTEX:         RunMutex       ( outIters , inIters ) ; break;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
