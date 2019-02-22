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
#define NEW_THREADS

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

enum
{
	ATOMIC_NONE ,
	ATOMIC_VOLATILE ,
	ATOMIC_ATOMIC ,
	ATOMIC_ATOMIC_FLOAT ,
	ATOMIC_ATOMIC_DOUBLE ,
	ATOMIC_MUTEX ,
	ATOMIC_REDUCTION ,
	ATOMIC_COUNT
};

const char * AtomicNames[] = { "none" , "volatile" , "atomic" , "atomic float" , "atomic double" , "mutex" , "reduction" };

cmdLineParameter< int > Count( "count" ) , Threads( "threads" , 1 );
cmdLineParameter< int > AtomicType( "atomic" , ATOMIC_NONE );

cmdLineReadable* params[] = { &Count , &Threads , &AtomicType , NULL };

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <itereation count (in thousands)>\n" , Count.name );
	printf( "\t[--%s <threads>=%d]\n" , Threads.name , Threads.value );
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

void RunReduction( size_t iterations )
{
	printf( "Reduction:\n" );
	auto Task = []( size_t i , size_t &sum ){ for( size_t j=0 ; j<i ; j++ ) sum += j; };

	{
		double t = Time();
		size_t sum = 0;
#pragma omp parallel for num_threads( MKThread::Threads ) reduction( + : sum )
		for( long long i=0 ; i<(long long)iterations ; i++ ) Task( i , sum );
		printf( "\tOpenMP:       %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		double t = Time();
		size_t sum = 0;
		std::vector< size_t > sums( MKThread::Threads , 0 );
		MKThread::parallel_for( 0 , iterations , [&]( unsigned int t , size_t i ){ Task( i , sums[t] ); } );
		for( unsigned int t=0 ; t<MKThread::Threads ; t++ ) sum += sums[t];
		printf( "\tMKThread:     %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		ThreadPool tq( MKThread::Threads );
		double t = Time();
		size_t sum = 0;
		std::vector< size_t > sums( MKThread::Threads , 0 );
		tq.parallel_for( 0 , iterations , [&]( unsigned int t , size_t i ){ Task( i , sums[t] ); } );
		for( unsigned int t=0 ; t<MKThread::Threads ; t++ ) sum += sums[t];
		printf( "\tThread Queue: %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
}

void RunMutex( size_t iterations )
{
	printf( "Mutex:\n" );
	auto Task1 = []( size_t i , size_t &sum )
	{
		for( size_t j=0 ; j<i ; j++ )
#pragma omp critical
			sum += j;
	};

	auto Task2 = []( size_t i , size_t &sum )
	{
		for( size_t j=0 ; j<i ; j++ )
		{
			static std::mutex m;
			std::lock_guard< std::mutex > lock(m);
			sum += j;
		}
	};

	{
		double t = Time();
		size_t sum = 0;
#pragma omp parallel for num_threads( MKThread::Threads )
		for( long long i=0 ; i<(long long)iterations ; i++ ) Task1( i , sum );
		printf( "\tOpenMP:       %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		double t = Time();
		size_t sum = 0;
		MKThread::parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task2( i , sum ); } );
		printf( "\tMKThread:     %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		ThreadPool tq( MKThread::Threads );
		double t = Time();
		size_t sum = 0;
		tq.parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task2( i , sum ); } );
		printf( "\tThread Queue: %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
}

void RunAtomicDouble( size_t iterations )
{
	printf( "Atomic double:\n" );
	auto Task1 = []( size_t i , double &sum )
	{
		double s = 0;
		for( size_t j=0 ; j<i ; j++ ) s += (double)j;
#pragma omp atomic
		sum += s;
	};

	auto Task2 = []( size_t i , double &sum )
	{
		double s = 0;
		for( size_t j=0 ; j<i ; j++ ) s += (double)j;
		AddAtomic( sum , s );
	};

	{
		double t = Time();
		double sum = 0;
#pragma omp parallel for num_threads( MKThread::Threads )
		for( long long i=0 ; i<(long long)iterations ; i++ ) Task1( i , sum );
		printf( "\tOpenMP:       %.2f(s)\t%g\n" , Time()-t , sum );
	}
	{
		double t = Time();
		double sum = 0;
		MKThread::parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task2( i , sum ); } );
		printf( "\tMKThread:     %.2f(s)\t%g\n" , Time()-t , sum );
	}
	{
		ThreadPool tq( MKThread::Threads );
		double t = Time();
		double sum = 0;
		tq.parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task2( i , sum ); } );
		printf( "\tThread Queue: %.2f(s)\t%g\n" , Time()-t , sum );
	}
}

void RunAtomicFloat( size_t iterations )
{
	printf( "Atomic float:\n" );
	auto Task1 = []( size_t i , float &sum )
	{
		float s = 0;
		for( size_t j=0 ; j<i ; j++ ) s += (float)j;
#pragma omp atomic
		sum += s;
	};

	auto Task2 = []( size_t i , float &sum )
	{
		float s = 0;
		for( size_t j=0 ; j<i ; j++ ) s += (float)j;
		AddAtomic( sum , s );
	};

	{
		double t = Time();
		float sum = 0;
#pragma omp parallel for num_threads( MKThread::Threads )
		for( long long i=0 ; i<(long long)iterations ; i++ ) Task1( i , sum );
		printf( "\tOpenMP:       %.2f(s)\t%g\n" , Time()-t , sum );
	}
	{
		double t = Time();
		float sum = 0;
		MKThread::parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task2( i , sum ); } );
		printf( "\tMKThread:     %.2f(s)\t%g\n" , Time()-t , sum );
	}
	{
		ThreadPool tq( MKThread::Threads );
		double t = Time();
		float sum = 0;
		tq.parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task2( i , sum ); } );
		printf( "\tThread Queue: %.2f(s)\t%g\n" , Time()-t , sum );
	}
}

void RunAtomic( size_t iterations )
{
	printf( "Atomic:\n" );
	auto Task1 = []( size_t i , size_t &sum )
	{
		size_t s = 0;
		for( size_t j=0 ; j<i ; j++ ) s += j;
#pragma omp atomic
		sum += s;
	};

	auto Task2 = []( size_t i , std::atomic< size_t > &sum )
	{
		size_t s = 0;
		for( size_t j=0 ; j<i ; j++ ) s += j;
		sum += s;
	};

	{
		double t = Time();
		size_t sum = 0;
#pragma omp parallel for num_threads( MKThread::Threads )
		for( long long i=0 ; i<(long long)iterations ; i++ ) Task1( i , sum );
		printf( "\tOpenMP:       %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		double t = Time();
		std::atomic< size_t > sum;
		sum.store(0);
		MKThread::parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task2( i , sum ); } );
		printf( "\tMKThread:     %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		ThreadPool tq( MKThread::Threads );
		double t = Time();
		std::atomic< size_t > sum;
		sum.store(0);
		tq.parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task2( i , sum ); } );
		printf( "\tThread Queue: %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
}

void RunVolatile( size_t iterations )
{
	printf( "Volatile:\n" );
	auto Task = []( size_t i , volatile size_t &sum ){ for( size_t j=0 ; j<i ; j++ ) sum += j; };

	{
		double t = Time();
		volatile size_t sum = 0;
#pragma omp parallel for num_threads( MKThread::Threads )
		for( long long i=0 ; i<(long long)iterations ; i++ ) Task( i , sum );
		printf( "\tOpenMP:       %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		double t = Time();
		volatile size_t sum = 0;
		MKThread::parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task( i , sum ); } );
		printf( "\tMKThread:     %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		ThreadPool tq( MKThread::Threads );
		double t = Time();
		volatile size_t sum = 0;
		tq.parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task( i , sum ); } );
		printf( "\tThread Queue: %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
}

void Run( size_t iterations )
{
	printf( "None:\n" );
	auto Task = []( size_t i , size_t &sum ){ for( size_t j=0 ; j<i ; j++ ) sum += j; };

	{
		double t = Time();
		size_t sum = 0;
#pragma omp parallel for num_threads( MKThread::Threads )
		for( long long i=0 ; i<(long long)iterations ; i++ ) Task( i , sum );
		printf( "\tOpenMP:   %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		double t = Time();
		size_t sum = 0;
		MKThread::parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task( i , sum ); } );
		printf( "\tMKThread: %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
	{
		ThreadPool tq( MKThread::Threads );
		double t = Time();
		size_t sum = 0;
		tq.parallel_for( 0 , iterations , [&]( unsigned int , size_t i ){ Task( i , sum ); } );
		printf( "\tMKThread: %.2f(s)\t%llu\n" , Time()-t , (unsigned long long)sum );
	}
}

int main( int argc , char* argv[] )
{

	cmdLineParse( argc-1 , &argv[1] , params );

	if( !Count.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	MKThread::Threads = Threads.value;
	size_t iterations = ( (size_t)Count.value )<<10;

	printf( "Iterations: %llu\n" , (unsigned long long)iterations );
	printf( "Threads: %d\n" , MKThread::Threads );

	if( AtomicType.value<0 )
	{
		Run            ( iterations );
		RunVolatile    ( iterations );
		RunAtomic      ( iterations );
		RunAtomicFloat ( iterations );
		RunAtomicDouble( iterations );
		RunMutex       ( iterations );
		RunReduction   ( iterations );
	}
	else
	{
		switch( AtomicType.value )
		{
		case ATOMIC_NONE:          Run            ( iterations ) ; break;
		case ATOMIC_VOLATILE:      RunVolatile    ( iterations ) ; break;
		case ATOMIC_ATOMIC:        RunAtomic      ( iterations ) ; break;
		case ATOMIC_ATOMIC_FLOAT:  RunAtomicFloat ( iterations ) ; break;
		case ATOMIC_ATOMIC_DOUBLE: RunAtomicDouble( iterations ) ; break;
		case ATOMIC_MUTEX:         RunMutex       ( iterations ) ; break;
		case ATOMIC_REDUCTION:     RunReduction   ( iterations ) ; break;
		}
	}
	return EXIT_SUCCESS;
}
