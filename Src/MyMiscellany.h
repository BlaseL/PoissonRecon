/*
Copyright (c) 2017, Michael Kazhdan
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
#ifndef MY_MISCELLANY_INCLUDED
#define MY_MISCELLANY_INCLUDED

#include "PreProcessor.h"

//////////////////
// OpenMP Stuff //
//////////////////
#ifdef _OPENMP
#include <omp.h>
#else // !_OPENMP
#ifdef NEW_THREADS
#else // !NEW_THREADS
inline int omp_get_num_procs  ( void ){ return 1; }
inline int omp_get_max_threads( void ){ return 1; }
inline int omp_get_thread_num ( void ){ return 0; }
inline void omp_set_num_threads( int ){}
inline void omp_set_nested( int ){}
struct omp_lock_t{};
inline void omp_init_lock( omp_lock_t* ){}
inline void omp_set_lock( omp_lock_t* ){}
inline void omp_unset_lock( omp_lock_t* ){}
inline void omp_destroy_lock( omp_lock_t* ){}
#endif // NEW_THREADS
#endif // _OPENMP

////////////////
// Time Stuff //
////////////////
#include <string.h>
#include <sys/timeb.h>
#ifndef WIN32
#include <sys/time.h>
#endif // WIN32

inline double Time( void )
{
#ifdef WIN32
	struct _timeb t;
	_ftime( &t );
	return double( t.time ) + double( t.millitm ) / 1000.0;
#else // WIN32
	struct timeval t;
	gettimeofday( &t , NULL );
	return t.tv_sec + double( t.tv_usec ) / 1000000;
#endif // WIN32
}

#include <cstdio>
#include <ctime>
#include <chrono>
struct Timer
{
	Timer( void ){ _startCPUClock = std::clock() , _startWallClock = std::chrono::system_clock::now(); }
	double cpuTime( void ) const{ return (std::clock() - _startCPUClock) / (double)CLOCKS_PER_SEC; };
	double wallTime( void ) const{  std::chrono::duration<double> diff = (std::chrono::system_clock::now() - _startWallClock) ; return diff.count(); }
protected:
	std::clock_t _startCPUClock;
	std::chrono::time_point< std::chrono::system_clock > _startWallClock;
};

///////////////
// I/O Stuff //
///////////////
#if defined( _WIN32 ) || defined( _WIN64 )
const char FileSeparator = '\\';
#else // !_WIN
const char FileSeparator = '/';
#endif // _WIN

#ifndef SetTempDirectory
#if defined( _WIN32 ) || defined( _WIN64 )
#define SetTempDirectory( tempDir , sz ) GetTempPath( (sz) , (tempDir) )
#else // !_WIN32 && !_WIN64
#define SetTempDirectory( tempDir , sz ) if( std::getenv( "TMPDIR" ) ) strcpy( tempDir , std::getenv( "TMPDIR" ) );
#endif // _WIN32 || _WIN64
#endif // !SetTempDirectory

#include <stdarg.h>
#include <vector>
#include <string>
struct MessageWriter
{
	char* outputFile;
	bool echoSTDOUT;
	MessageWriter( void ){ outputFile = NULL , echoSTDOUT = true; }
	void operator() ( const char* format , ... )
	{
		if( outputFile )
		{
			FILE* fp = fopen( outputFile , "a" );
			va_list args;
			va_start( args , format );
			vfprintf( fp , format , args );
			fclose( fp );
			va_end( args );
		}
		if( echoSTDOUT )
		{
			va_list args;
			va_start( args , format );
			vprintf( format , args );
			va_end( args );
		}
	}
	void operator() ( std::vector< char* >& messages  , const char* format , ... )
	{
		if( outputFile )
		{
			FILE* fp = fopen( outputFile , "a" );
			va_list args;
			va_start( args , format );
			vfprintf( fp , format , args );
			fclose( fp );
			va_end( args );
		}
		if( echoSTDOUT )
		{
			va_list args;
			va_start( args , format );
			vprintf( format , args );
			va_end( args );
		}
		// [WARNING] We are not checking the string is small enough to fit in 1024 characters
		messages.push_back( new char[1024] );
		char* str = messages.back();
		va_list args;
		va_start( args , format );
		vsprintf( str , format , args );
		va_end( args );
		if( str[strlen(str)-1]=='\n' ) str[strlen(str)-1] = 0;
	}
	void operator() ( std::vector< std::string >& messages  , const char* format , ... )
	{
		if( outputFile )
		{
			FILE* fp = fopen( outputFile , "a" );
			va_list args;
			va_start( args , format );
			vfprintf( fp , format , args );
			fclose( fp );
			va_end( args );
		}
		if( echoSTDOUT )
		{
			va_list args;
			va_start( args , format );
			vprintf( format , args );
			va_end( args );
		}
		// [WARNING] We are not checking the string is small enough to fit in 1024 characters
		char message[1024];
		va_list args;
		va_start( args , format );
		vsprintf( message , format , args );
		va_end( args );
		if( message[strlen(message)-1]=='\n' ) message[strlen(message)-1] = 0;
		messages.push_back( std::string( message ) );
	}
};

/////////////////////////////////////
// Exception, Warnings, and Errors //
/////////////////////////////////////
#include <exception>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
namespace MKExceptions
{
	template< typename ... Arguments > void _AddToMessageStream( std::stringstream &stream , Arguments ... arguments );
	inline void _AddToMessageStream( std::stringstream &stream ){ return; }
	template< typename Argument , typename ... Arguments > void _AddToMessageStream( std::stringstream &stream , Argument argument , Arguments ... arguments )
	{
		stream << argument;
		_AddToMessageStream( stream , arguments ... );
	}

	template< typename ... Arguments >
	std::string MakeMessageString( std::string header , std::string fileName , int line , std::string functionName , Arguments ... arguments )
	{
		size_t headerSize = header.size();
		std::stringstream stream;

		// The first line is the header, the file name , and the line number
		stream << header << " " << fileName << " (Line " << line << ")" << std::endl;

		// Inset the second line by the size of the header and write the function name
		for( size_t i=0 ; i<=headerSize ; i++ ) stream << " ";
		stream << functionName << std::endl;

		// Inset the third line by the size of the header and write the rest
		for( size_t i=0 ; i<=headerSize ; i++ ) stream << " ";
		_AddToMessageStream( stream , arguments ... );

		return stream.str();
	}
	struct Exception : public std::exception
	{
		const char *what( void ) const noexcept { return _message.c_str(); }
		template< typename ... Args >
		Exception( const char *fileName , int line , const char *functionName , const char *format , Args ... args )
		{
			_message = MakeMessageString( "[EXCEPTION]" , fileName , line , functionName , format , args ... );
		}
	private:
		std::string _message;
	};

	template< typename ... Args > void Throw( const char *fileName , int line , const char *functionName , const char *format , Args ... args ){ throw Exception( fileName , line , functionName , format , args ... ); }
	template< typename ... Args >
	void Warn( const char *fileName , int line , const char *functionName , const char *format , Args ... args )
	{
		std::cerr << MakeMessageString( "[WARNING]" , fileName , line , functionName , format , args ... ) << std::endl;
	}
	template< typename ... Args >
	void ErrorOut( const char *fileName , int line , const char *functionName , const char *format , Args ... args )
	{
		std::cerr << MakeMessageString( "[ERROR]" , fileName , line , functionName , format , args ... ) << std::endl;
		exit( 0 );
	}
}
#ifndef WARN
#define WARN( ... ) MKExceptions::Warn( __FILE__ , __LINE__ , __FUNCTION__ , __VA_ARGS__ )
#endif // WARN
#ifndef WARN_ONCE
#define WARN_ONCE( ... ) { static bool firstTime = true ; if( firstTime ) MKExceptions::Warn( __FILE__ , __LINE__ , __FUNCTION__ , __VA_ARGS__ ) ; firstTime = false; }
#endif // WARN_ONCE
#ifndef THROW
#define THROW( ... ) MKExceptions::Throw( __FILE__ , __LINE__ , __FUNCTION__ , __VA_ARGS__ )
#endif // THROW
#ifndef ERROR_OUT
#define ERROR_OUT( ... ) MKExceptions::ErrorOut( __FILE__ , __LINE__ , __FUNCTION__ , __VA_ARGS__ )
#endif // ERROR_OUT

#ifdef NEW_CODE

#include <signal.h>
#if defined(_WIN32) || defined( _WIN64 )
inline void StackTrace( void )
{
}
#else // !WINDOWS
#include <execinfo.h>
#include <unistd.h>
#include <cxxabi.h>
#include <mutex>

inline void StackTrace( void )
{
	static std::mutex mutex;
	std::lock_guard< std::mutex > lock(mutex);

	// Code borrowed from:
	// https://stackoverflow.com/questions/77005/how-to-automatically-generate-a-stacktrace-when-my-program-crashes
	// and
	// https://stackoverflow.com/questions/15129089/is-there-a-way-to-dump-stack-trace-with-line-number-from-a-linux-release-binary/15130037
	void * trace[128];
	int size = backtrace( trace , 128 );

	char ** messages = backtrace_symbols( trace , size );
	for( int i=1 ; i< size && messages!=NULL ; ++i )
	{
		char *mangled_name=0 , *offset_begin=0 , *offset_end=0;

		char syscom[1024];
		sprintf( syscom , "addr2line %p -e PoissonRecon/Bin/Linux/PoissonRecon" , trace[i] ); //last parameter is the name of this app
		if( !system( syscom ) ) ;

		// find parantheses and +address offset surrounding mangled name
		for( char *p=messages[i] ; *p ; ++p )
		{
			if     ( *p=='(' ) mangled_name = p; 
			else if( *p=='+' ) offset_begin = p;
			else if( *p==')' )
			{
				offset_end = p;
				break;
			}
		}

		// if the line could be processed, attempt to demangle the symbol
		if( mangled_name && offset_begin && offset_end && mangled_name<offset_begin )
		{
			*mangled_name++ = '\0';
			*offset_begin++ = '\0';
			*offset_end++ = '\0';

			int status;
			char * real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);

			// if demangling is successful, output the demangled function name
			if( !status )
			{
				std::cerr << "\t(" << i << ") " << messages[i] << " : " << real_name << "+" << offset_begin << offset_end  << std::endl;
				std::cout << "\t(" << i << ") " << messages[i] << " : " << real_name << "+" << offset_begin << offset_end  << std::endl;
			}
			// otherwise, output the mangled function name
			else
			{
				std::cerr << "\t(" << i << ") " << messages[i] << " : " << mangled_name << "+" << offset_begin << offset_end << std::endl;
				std::cout << "\t(" << i << ") " << messages[i] << " : " << mangled_name << "+" << offset_begin << offset_end << std::endl;
			}
			free( real_name );
		}
		// otherwise, print the whole line
		else
		{
			std::cerr << "\t(" << i << ") " << messages[i] << std::endl;
			std::cout << "\t(" << i << ") " << messages[i] << std::endl;
		}
	}

	free( messages );
}
#endif // WINDOWS
inline void SignalHandler( int signal )
{
	printf( "Signal: %d\n" , signal );
	StackTrace();
	exit( 0 );
}


template< typename Value > bool SetAtomic( volatile Value *value , Value newValue , Value oldValue );
template< typename Data > void AddAtomic( Data& a , Data b );

#ifdef NEW_THREADS
////////////////////
// MKThread Stuff //
////////////////////
#include <thread>
#include <mutex>
#include <vector>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <chrono>
#include <future>

#include <memory>
namespace hh
{
	inline int get_max_threads( void ) { return std::max< int >( (int)std::thread::hardware_concurrency() , 1 ); }

	class ThreadPoolIndexedTask
	{
	public:
		using Task = std::function< void (int) >;
		ThreadPoolIndexedTask( void )
		{
#ifdef USE_THREAD_MINUS_ONE
			const int num_threads = get_max_threads()-1;
#else // !USE_THREAD_MINUS_ONE
			const int num_threads = get_max_threads();
#endif // USE_THREAD_MINUS_ONE
			_threads.reserve( num_threads );
			for( int i=0 ; i<num_threads ; i++ ) _threads.emplace_back( &ThreadPoolIndexedTask::worker_main , this );
		}
		~ThreadPoolIndexedTask( void )
		{
printf( "ThreadPoolIndexedTask dtor called\n" );
#ifdef USE_THREAD_MINUS_ONE
			if( _threads.size() )
#endif // USE_THREAD_MINUS_ONE
			{
				std::unique_lock< std::mutex > lock( _mutex );
				if( !_running ) ERROR_OUT( "not running" );
				if( _num_remaining_tasks ) ERROR_OUT( "tasks remaining" );
				if( _task_index!=_num_tasks ) ERROR_OUT( "task index and num tasks don't match" );
				_running = false;
				_task_index = 0;
				_num_tasks = 1;
				_condition_variable_worker.notify_all();
			}
			for( int i=0 ; i<_threads.size() ; i++ ) _threads[i].join();
		}
#ifdef USE_THREAD_MINUS_ONE
		int num_threads( void ) const { return (int)_threads.size()+1; }
#else // !USE_THREAD_MINUS_ONE
		int num_threads( void ) const { return (int)_threads.size(); }
#endif // USE_THREAD_MINUS_ONE
		bool already_active( void ) const { return _num_remaining_tasks!=0; }  // detect nested execution
		void execute( int num_tasks , const Task& task_function )
		{
			if( already_active() )
			{
				WARN( "Nested execution of ThreadPoolIndexedTask is run serially" );
				for( int i=0 ; i<num_tasks ; i++ ) task_function( i );
			}
			else
			{
				std::unique_lock< std::mutex > lock( _mutex );
				_task_function = task_function;
				_num_tasks = num_tasks;
				_num_remaining_tasks = num_tasks;
				_task_index = 0;
				_condition_variable_worker.notify_all();
				_condition_variable_master.wait( lock , [this]{ return !_num_remaining_tasks; } );
			}
		}
		static ThreadPoolIndexedTask& default_threadpool( void )
		{
			static std::unique_ptr< ThreadPoolIndexedTask > thread_pool;
			// This is safe because thread_pool is nullptr only in the main thread before any other thread is launched.
			if( !thread_pool ) thread_pool = std::unique_ptr< ThreadPoolIndexedTask >( new ThreadPoolIndexedTask() );
			return *thread_pool;
		}

	private:
		std::mutex _mutex;
		bool _running = true;
		std::vector< std::thread > _threads;
		Task _task_function;
		int _num_tasks = 0;
		int _num_remaining_tasks = 0;
		int _task_index = 0;
		std::condition_variable _condition_variable_worker;
		std::condition_variable _condition_variable_master;

		void worker_main( void )
		{
			std::unique_lock< std::mutex > lock( _mutex );
			// Consider: https://stackoverflow.com/questions/233127/how-can-i-propagate-exceptions-between-threads
			// However, rethrowing the exception in the main thread loses the stack state, so not useful for debugging.
			for (;;)
			{
				_condition_variable_worker.wait( lock , [this]{ return _task_index < _num_tasks; } );
				if( !_running ) break;
				while( _task_index<_num_tasks )
				{
					int i = _task_index++;
					lock.unlock();
					_task_function(i);
					lock.lock();
					if( _num_remaining_tasks<=0 ) ERROR_OUT( "num remaining tasks not greater than zero" );
					if (!--_num_remaining_tasks) _condition_variable_master.notify_all();
				}
			}
		}
	};
}

struct ThreadPool
{
	enum ParallelType
	{
#ifdef _OPENMP
		OPEN_MP ,
#endif // _OPENMP
		THREAD_POOL ,
		THREAD_POOL_HH ,
		ASYNC ,
		NONE
	};
	static const std::vector< std::string > ParallelNames;

	enum ScheduleType
	{
		STATIC ,
		DYNAMIC
	};
	static const std::vector< std::string > ScheduleNames;

	static size_t DefaultChunkSize;
	static ScheduleType DefaultSchedule;

	template< typename ... Functions >
	static void ParallelSections( const Functions & ... functions )
	{
		std::vector< std::future< void > > futures( sizeof...(Functions) );
		_ParallelSections( &futures[0] , functions ... );
		for( size_t t=0 ; t<futures.size() ; t++ ) futures[t].get();
	}

	static void Parallel_for( size_t begin , size_t end , const std::function< void ( unsigned int , size_t ) > &iterationFunction , ScheduleType schedule=DefaultSchedule , size_t chunkSize=DefaultChunkSize )
	{
		if( begin>=end ) return;
		size_t range = end - begin;
		size_t chunks = ( range + chunkSize - 1 ) / chunkSize;
#ifdef USE_FEWER_THREADS
		unsigned int threads = (unsigned int)std::min< size_t >( NumThreads() , chunks );
#else // !USE_FEWER_THREADS
		unsigned int threads = (unsigned int)NumThreads();
#endif // USE_FEWER_THREADS
		std::atomic< size_t > index;
		index.store( 0 );


#ifdef USE_FEWER_THREADS
		if( _ParallelType==NONE || threads==1 )
#else // !USE_FEWER_THREADS
		if( range<chunkSize || _ParallelType==NONE || threads==1 )
#endif // USE_FEWER_THREADS
		{
			for( size_t i=begin ; i<end ; i++ ) iterationFunction( 0 , i );
			return;
		}

		auto _ChunkFunction = [ &iterationFunction , begin , end , chunkSize ]( unsigned int thread , size_t chunk )
		{
			const size_t _begin = begin + chunkSize*chunk;
			const size_t _end = std::min< size_t >( end , _begin+chunkSize );
			for( size_t i=_begin ; i<_end ; i++ ) iterationFunction( thread , i );
		};
		auto _StaticThreadFunction = [ &_ChunkFunction , chunks , threads ]( unsigned int thread )
		{
			for( size_t chunk=thread ; chunk<chunks ; chunk+=threads ) _ChunkFunction( thread , chunk );
		};
		auto _DynamicThreadFunction = [ &_ChunkFunction , chunks , &index ]( unsigned int thread )
		{
			size_t chunk;
			while( ( chunk=index.fetch_add(1) )<chunks ) _ChunkFunction( thread , chunk );
		};

		if     ( schedule==STATIC  ) _ThreadFunction = _StaticThreadFunction;
		else if( schedule==DYNAMIC ) _ThreadFunction = _DynamicThreadFunction;

		if( false ){}
#ifdef _OPENMP
		else if( _ParallelType==OPEN_MP )
		{
			if( schedule==STATIC )
#pragma omp parallel for num_threads( threads ) schedule( static , 1 )
				for( int c=0 ; c<chunks ; c++ ) _ChunkFunction( omp_get_thread_num() , c );
			else if( schedule==DYNAMIC )
#pragma omp parallel for num_threads( threads ) schedule( dynamic , 1 )
				for( int c=0 ; c<chunks ; c++ ) _ChunkFunction( omp_get_thread_num() , c );
		}
#endif // _OPENMP
		else if( _ParallelType=THREAD_POOL_HH )
		{
#ifdef USE_FEWER_THREADS
#else // !USE_FEWER_THREADS
			const int num_threads = (int)std::min< size_t >( threads , range );
#endif // USE_FEWER_THREADS
			hh::ThreadPoolIndexedTask* const thread_pool = &hh::ThreadPoolIndexedTask::default_threadpool();
			if( !thread_pool || thread_pool->already_active() )
			{
				// Traverse the range elements sequentially.
				for( size_t index=begin ; index<end ; ++index ) iterationFunction( 0 , index );
			}
#ifdef USE_FEWER_THREADS
#ifdef USE_THREAD_MINUS_ONE
			else
			{
				thread_pool->execute( threads-1 , _ThreadFunction );
				_ThreadFunction( threads-1 );
			}
#else // !USE_THREAD_MINUS_ONE
			else thread_pool->execute( threads , _ThreadFunction );
#endif // USE_THREAD_MINUS_ONE
#else // !USE_FEWER_THREADS
			else thread_pool->execute( num_threads , _ThreadFunction );
#endif // USE_FEWER_THREADS
		}
		else if( _ParallelType==ASYNC )
		{
			static std::vector< std::future< void > > futures;
#ifdef USE_THREAD_MINUS_ONE
			futures.resize( threads-1 );
			for( unsigned int t=1 ; t<threads ; t++ ) futures[t-1] = std::async( std::launch::async , _ThreadFunction , t );
			_ThreadFunction( 0 );
			for( unsigned int t=1 ; t<threads ; t++ ) futures[t-1].get();
#else // !USE_THREAD_MINUS_ONE
			futures.resize( threads );
			for( unsigned int t=0 ; t<threads ; t++ ) futures[t] = std::async( std::launch::async , _ThreadFunction , t );
			for( unsigned int t=0 ; t<threads ; t++ ) futures[t].get();
#endif // USE_THREAD_MINUS_ONE
		}
		else if( _ParallelType==THREAD_POOL )
		{
#ifdef USE_FEWER_THREADS
#if 1
#ifdef USE_THREAD_MINUS_ONE
			unsigned int targetTasks = 0;
			if( !SetAtomic( &_RemainingTasks , threads-1 , targetTasks ) )
#else // !USE_THREAD_MINUS_ONE
			if( !SetAtomic( &_RemainingTasks , threads , 0 ) )
#endif // USE_THREAD_MINUS_ONE
#else
			if( _RemainingTasks )
#endif
#else // !USE_FEWER_THREADS
			if( _WorkToBeDone )
#endif // USE_FEWER_THREADS
			{
				WARN( "nested for loop, reverting to serial" );
				for( size_t i=begin ; i<end ; i++ ) iterationFunction( 0 , i );
			}
			else
			{
#ifdef USE_FEWER_THREADS
				std::unique_lock< std::mutex > lock( _Mutex );
#ifdef USE_THREAD_MINUS_ONE
#if 1
				_TotalTasks = threads-1;
#else
				_TotalTasks = _RemainingTasks = threads-1;
#endif
				_NextTask = 0;
				for( unsigned int t=0 ; t<threads-1 ; t++ ) _WaitingForWorkOrClose.notify_one();
				lock.unlock();
				_ThreadFunction( threads-1 );
				lock.lock();
#else // !USE_THREAD_MINUS_ONE
#if 1
				_TotalTasks = threads;
#else
				_TotalTasks = _RemainingTasks = threads;
#endif
				_NextTask = 0;
				for( unsigned int t=0 ; t<threads ; t++ ) _WaitingForWorkOrClose.notify_one();
#endif // USE_THREAD_MINUS_ONE
				_DoneWithWork.wait( lock , [&]( void ){ return _RemainingTasks==0; } );
#else // !USE_FEWER_THREADS
				_WorkToBeDone = (unsigned int)_threads.size();
				_WaitingForWorkOrClose.notify_all();
				{
					std::unique_lock< std::mutex > lock( _mutex );
					_DoneWithWork.wait( lock , [&]( void ){ return _WorkToBeDone==0; } );
				}
#endif // USE_FEWER_THREADS
			}
		}
	}

#ifdef USE_THREAD_MINUS_ONE
	static unsigned int NumThreads( void ){ return (unsigned int)_Threads.size()+1; }
#else // !USE_THREAD_MINUS_ONE
	static unsigned int NumThreads( void ){ return (unsigned int)_Threads.size(); }
#endif // USE_THREAD_MINUS_ONE

	static void Init( ParallelType parallelType , unsigned int numThreads=std::thread::hardware_concurrency() )
	{
		_ParallelType = parallelType;
		if( _Threads.size() && !_Close )
		{
			_Close = true;
			_WaitingForWorkOrClose.notify_all();
			for( unsigned int t=0 ; t<_Threads.size() ; t++ ) _Threads[t].join();
		}
		_Close = true;
#ifdef USE_THREAD_MINUS_ONE
		numThreads--;
#endif // USE_THREAD_MINUS_ONE
		_Threads.resize( numThreads );
		if( _ParallelType==THREAD_POOL )
		{
#ifdef USE_FEWER_THREADS
			_RemainingTasks = 0;
#else // !USE_FEWER_THREADS
			_WorkToBeDone = 0;
#endif // USE_FEWER_THREADS
			_Close = false;
#ifdef USE_FEWER_THREADS
			for( unsigned int t=0 ; t<numThreads ; t++ ) _Threads[t] = std::thread( _ThreadInitFunction );
#else // !USE_FEWER_THREADS
			for( unsigned int t=0 ; t<numThreads ; t++ ) _Threads[t] = std::thread( _ThreadInitFunction , t );
#endif // USE_FEWER_THREADS
		}
	}
	static void Terminate( void )
	{
		if( _Threads.size() && !_Close )
		{
			_Close = true;
			_WaitingForWorkOrClose.notify_all();
			for( unsigned int t=0 ; t<_Threads.size() ; t++ ) _Threads[t].join();
			_Threads.resize( 0 );
		}
	}
private:
	ThreadPool( const ThreadPool & ){}
	ThreadPool &operator = ( const ThreadPool & ){}

	template< typename Function >
	static void _ParallelSections( std::future< void > *futures , const Function &function ){ *futures = std::async( std::launch::async , function ); }
	template< typename Function , typename ... Functions >
	static void _ParallelSections( std::future< void > *futures , const Function &function , const Functions& ... functions )
	{
		*futures = std::async( std::launch::async , function );
		_ParallelSections( futures+1 , functions ... );
	}
#ifdef USE_FEWER_THREADS
	static void _ThreadInitFunction( void )
#else // !USE_FEWER_THREADS
	static void _ThreadInitFunction( unsigned int thread )
#endif // USE_FEWER_THREADS
	{
		// Wait for the first job to come in
		std::unique_lock< std::mutex > lock( _Mutex );
		_WaitingForWorkOrClose.wait( lock );
		while( !_Close )
		{
#ifdef USE_FEWER_THREADS
#if 1
			if( _NextTask<_TotalTasks )
#else
			while( _NextTask<_TotalTasks )
#endif
			{
				unsigned int currentTask = _NextTask++;
				lock.unlock();
				// do the job
				_ThreadFunction( currentTask );
				lock.lock();
				_RemainingTasks--;
			}
			if( !_RemainingTasks ) _DoneWithWork.notify_all();
			_WaitingForWorkOrClose.wait( lock );
#else // !USE_FEWER_THREADS
			lock.unlock();
			// do the job
			_ThreadFunction( thread );

			// Notify and wait for the next job
			lock.lock();
			_WorkToBeDone--;
			if( !_WorkToBeDone ) _DoneWithWork.notify_all();
			_WaitingForWorkOrClose.wait( lock );
#endif // USE_FEWER_THREADS
		}
	}

	static bool _Close;
#ifdef USE_FEWER_THREADS
	static unsigned int _TotalTasks;
	static volatile unsigned int _RemainingTasks;
	static unsigned int _NextTask;
#else // !USE_FEWER_THREADS
	static unsigned int _WorkToBeDone;
#endif // USE_FEWER_THREADS
	static std::mutex _Mutex;
	static std::condition_variable _WaitingForWorkOrClose , _DoneWithWork;
	static std::vector< std::thread > _Threads;
	static std::function< void ( unsigned int ) > _ThreadFunction;
	static ParallelType _ParallelType;
};

size_t ThreadPool::DefaultChunkSize = 128;
ThreadPool::ScheduleType ThreadPool::DefaultSchedule = ThreadPool::DYNAMIC;
bool ThreadPool::_Close;
#ifdef USE_FEWER_THREADS
unsigned int ThreadPool::_TotalTasks;
volatile unsigned int ThreadPool::_RemainingTasks;
unsigned int ThreadPool::_NextTask;
#else // !USE_FEWER_THREADS
unsigned int ThreadPool::_WorkToBeDone;
#endif // USE_FEWER_THREADS
std::mutex ThreadPool::_Mutex;
std::condition_variable ThreadPool::_WaitingForWorkOrClose;
std::condition_variable ThreadPool::_DoneWithWork;
std::vector< std::thread > ThreadPool::_Threads;
std::function< void ( unsigned int ) > ThreadPool::_ThreadFunction;
ThreadPool::ParallelType ThreadPool::_ParallelType;

const std::vector< std::string >ThreadPool::ParallelNames =
{
#ifdef _OPENMP
	"open mp" ,
#endif // _OPENMP
	"thread pool" ,
	"thread pool (hh)" ,
	"async" ,
	"none"
};
const std::vector< std::string >ThreadPool::ScheduleNames = { "static" , "dynamic" };


#endif // NEW_THREADS

#ifdef NEW_CODE
#include <mutex>

#ifdef NEW_THREADS
#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#endif // _WIN32 || _WIN64
#endif // NEW_THREADS
template< typename Value >
bool SetAtomic32( volatile Value *value , Value newValue , Value oldValue )
{
#if defined( _WIN32 ) || defined( _WIN64 )
	long &_oldValue = *(long *)&oldValue;
	long &_newValue = *(long *)&newValue;
	return InterlockedCompareExchange( (long*)value , _newValue , _oldValue )==_oldValue;
#else // !_WIN32 && !_WIN64
	uint32_t &_oldValue = *(uint32_t *)&oldValue;
	uint32_t &_newValue = *(uint32_t *)&newValue;
//	return __sync_bool_compare_and_swap( (uint32_t *)value , _oldValue , _newValue );
	return __atomic_compare_exchange_n( (uint32_t *)value , (uint32_t *)&oldValue , _newValue , false , __ATOMIC_SEQ_CST , __ATOMIC_SEQ_CST );
#endif // _WIN32 || _WIN64
}
template< typename Value >
bool SetAtomic64( volatile Value *value , Value newValue , Value oldValue )
{
#if defined( _WIN32 ) || defined( _WIN64 )
	__int64 &_oldValue = *(__int64 *)&oldValue;
	__int64 &_newValue = *(__int64 *)&newValue;
	return InterlockedCompareExchange64( (__int64*)value , _newValue , _oldValue )==_oldValue;
#else // !_WIN32 && !_WIN64
	uint64_t &_oldValue = *(uint64_t *)&oldValue;
	uint64_t &_newValue = *(uint64_t *)&newValue;
//	return __sync_bool_compare_and_swap ( (uint64_t *)&value , _oldValue , _newValue );
	return __atomic_compare_exchange_n( (uint64_t *)value , (uint64_t *)&oldValue , _newValue , false , __ATOMIC_SEQ_CST , __ATOMIC_SEQ_CST );
#endif // _WIN32 || _WIN64
}

template< typename Number >
void AddAtomic32( Number &a , Number b )
{
#ifdef NEW_THREADS
#if 0
	Number current = a;
	Number sum = current+b;
	while( !SetAtomic32( &a , sum , current ) ) current = a , sum = a+b;
#else
#if defined( _WIN32 ) || defined( _WIN64 )
	Number current = a;
	Number sum = current+b;
	long &_current = *(long *)&current;
	long &_sum = *(long *)&sum;
	while( InterlockedCompareExchange( (long*)&a , _sum , _current )!=_current ) current = a , sum = a+b;
#else // !_WIN32 && !_WIN64
	Number current = a;
	Number sum = current+b;
	uint32_t &_current = *(uint32_t *)&current;
	uint32_t &_sum = *(uint32_t *)&sum;
	while( __sync_val_compare_and_swap( (uint32_t *)&a , _current , _sum )!=_current ) current = a , sum = a+b;
#endif // _WIN32 || _WIN64
#endif
#else // !NEW_THREADS
#pragma omp atomic
	a += b;
#endif // NEW_THREADS
}

template< typename Number >
void AddAtomic64( Number &a , Number b )
{
#ifdef NEW_THREADS
#if 1
	Number current = a;
	Number sum = current+b;
	while( !SetAtomic64( &a , sum , current ) ) current = a , sum = a+b;
#else
#if defined( _WIN32 ) || defined( _WIN64 )
	Number current = a;
	Number sum = current+b;
	__int64 &_current = *(__int64 *)&current;
	__int64 &_sum = *(__int64 *)&sum;
	while( InterlockedCompareExchange64( (__int64*)&a , _sum , _current )!=_current ) current = a , sum = a+b;
#else // !_WIN32 && !_WIN64
	Number current = a;
	Number sum = current+b;
	uint64_t &_current = *(uint64_t *)&current;
	uint64_t &_sum = *(uint64_t *)&sum;
	while( __sync_val_compare_and_swap( (uint64_t *)&a , _current , _sum )!=_current ) current = a , sum = a+b;
#endif // _WIN32 || _WIN64
#endif
#else // !NEW_THREADS
#pragma omp atomic
	a += b;
#endif // NEW_THREADS
}

template< typename Value >
bool SetAtomic( volatile Value *value , Value newValue , Value oldValue )
{
	switch( sizeof(Value) )
	{
	case 4: return SetAtomic32( value , newValue , oldValue );
	case 8: return SetAtomic64( value , newValue , oldValue );
	default:
		WARN_ONCE( "should not use this function: " , sizeof(Value) );
		static std::mutex setAtomicMutex;
		std::lock_guard< std::mutex > lock( setAtomicMutex );
		if( *value==oldValue ){ *value = newValue ; return true; }
		else return false;
	}
}

template< typename Data >
void AddAtomic( Data& a , Data b )
{
	switch( sizeof(Data) )
	{
	case 4: return AddAtomic32( a , b );
	case 8: return AddAtomic64( a , b );
	default:
		WARN_ONCE( "should not use this function: " , sizeof(Data) );
#ifdef NEW_THREADS
		static std::mutex addAtomicMutex;
		std::lock_guard< std::mutex > lock( addAtomicMutex );
		a += b;
#else // !NEW_THREADS
#pragma omp critical
		a += b;
#endif // NEW_THREADS
	}
}
#endif // NEW_CODE
/////////////////////////
// NumberWrapper Stuff //
/////////////////////////
#include <vector>
struct EmptyNumberWrapperClass{};

template< typename Number , typename Type=EmptyNumberWrapperClass , size_t I=0 >
struct NumberWrapper
{
	typedef Number type;

	Number n;

	NumberWrapper( Number _n=0 ) : n(_n){}
	NumberWrapper operator + ( NumberWrapper _n ) const { return NumberWrapper( n + _n.n ); }
	NumberWrapper operator - ( NumberWrapper _n ) const { return NumberWrapper( n - _n.n ); }
	NumberWrapper operator * ( NumberWrapper _n ) const { return NumberWrapper( n * _n.n ); }
	NumberWrapper operator / ( NumberWrapper _n ) const { return NumberWrapper( n / _n.n ); }
	NumberWrapper &operator += ( NumberWrapper _n ){ n += _n.n ; return *this; }
	NumberWrapper &operator -= ( NumberWrapper _n ){ n -= _n.n ; return *this; }
	NumberWrapper &operator *= ( NumberWrapper _n ){ n *= _n.n ; return *this; }
	NumberWrapper &operator /= ( NumberWrapper _n ){ n /= _n.n ; return *this; }
	bool operator == ( NumberWrapper _n ) const { return n==_n.n; }
	bool operator != ( NumberWrapper _n ) const { return n!=_n.n; }
	bool operator <  ( NumberWrapper _n ) const { return n<_n.n; }
	bool operator >  ( NumberWrapper _n ) const { return n>_n.n; }
	bool operator <= ( NumberWrapper _n ) const { return n<=_n.n; }
	bool operator >= ( NumberWrapper _n ) const { return n>=_n.n; }
	NumberWrapper operator ++ ( int ) { NumberWrapper _n(n) ; n++ ; return _n; }
	NumberWrapper operator -- ( int ) { NumberWrapper _n(n) ; n-- ; return _n; }
	NumberWrapper &operator ++ ( void ) { n++ ; return *this; }
	NumberWrapper &operator -- ( void ) { n-- ; return *this; }
	explicit operator Number () const { return n; }
};

#if 0
template< typename Number , typename Type , size_t I >
struct std::atomic< NumberWrapper< Number , Type , I > >
{
	typedef Number type;

	std::atomic< Number > n;

	atomic( Number _n=0 ) : n(_n){}
	atomic( const std::atomic< Number > &_n ) : n(_n){}
	atomic( NumberWrapper< Number , Type , I > _n ) : n(_n.n){}
	atomic &operator = ( Number _n ){ n = _n ; return *this; }
//	atomic &operator = ( const atomic &a ){ n = a.n ; return *this; }
//	atomic &operator = ( const NumberWrapper< Number , Type , I > &_n ){ n = _n.n ; return *this; }
	atomic operator + ( atomic _n ) const { return atomic( n + _n.n ); }
	atomic operator - ( atomic _n ) const { return atomic( n * _n.n ); }
	atomic operator * ( atomic _n ) const { return atomic( n * _n.n ); }
	atomic operator / ( atomic _n ) const { return atomic( n / _n.n ); }
	atomic &operator += ( atomic _n ){ n += _n.n ; return *this; }
	atomic &operator -= ( atomic _n ){ n -= _n.n ; return *this; }
	atomic &operator *= ( atomic _n ){ n *= _n.n ; return *this; }
	atomic &operator /= ( atomic _n ){ n /= _n.n ; return *this; }
	bool operator == ( atomic _n ) const { return n==_n.n; }
	bool operator != ( atomic _n ) const { return n!=_n.n; }
	bool operator <  ( atomic _n ) const { return n<_n.n; }
	bool operator >  ( atomic _n ) const { return n>_n.n; }
	bool operator <= ( atomic _n ) const { return n<=_n.n; }
	bool operator >= ( atomic _n ) const { return n>=_n.n; }
	atomic operator ++ ( int ) { atomic _n(n) ; n++ ; return _n; }
	atomic operator -- ( int ) { atomic _n(n) ; n-- ; return _n; }
	atomic &operator ++ ( void ) { n++ ; return *this; }
	atomic &operator -- ( void ) { n-- ; return *this; }
	operator NumberWrapper< Number , Type , I >() const { return NumberWrapper< Number , Type , I >(n); }
	explicit operator Number () const { return n; }
};
#endif


namespace std
{
	template< typename Number , typename Type , size_t I >
	struct hash< NumberWrapper< Number , Type , I > >
	{
		size_t operator()( NumberWrapper< Number , Type , I > n ) const { return std::hash< Number >{}( n.n ); }
	};
}

template< typename Data , typename _NumberWrapper >
struct VectorWrapper : public std::vector< Data >
{
	VectorWrapper( void ){}
	VectorWrapper( size_t sz ) : std::vector< Data >( sz ){}
	VectorWrapper( size_t sz , Data d ) : std::vector< Data >( sz , d ){}

//	void resize( _NumberWrapper n )         { std::vector< Data >::resize( (size_t)(_NumberWrapper::type)n ); }
//	void resize( _NumberWrapper n , Data d ){ std::vector< Data >::resize( (size_t)(_NumberWrapper::type)n , d ); }

	typename std::vector< Data >::reference operator[]( _NumberWrapper n ){ return std::vector< Data >::operator[]( n.n ); }
	typename std::vector< Data >::const_reference operator[]( _NumberWrapper n ) const { return std::vector< Data >::operator[]( n.n ); }
};
#endif // NEW_CODE

//////////////////
// Memory Stuff //
//////////////////
size_t getPeakRSS( void );
size_t getCurrentRSS( void );

struct MemoryInfo
{
	static size_t Usage( void ){ return getCurrentRSS(); }
	static int PeakMemoryUsageMB( void ){ return (int)( getPeakRSS()>>20 ); }
};
#if defined( _WIN32 ) || defined( _WIN64 )
#include <Windows.h>
#include <Psapi.h>
inline void SetPeakMemoryMB( size_t sz )
{
	sz <<= 20;
	SIZE_T peakMemory = sz;
	HANDLE h = CreateJobObject( NULL , NULL );
	AssignProcessToJobObject( h , GetCurrentProcess() );

	JOBOBJECT_EXTENDED_LIMIT_INFORMATION jeli = { 0 };
	jeli.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_JOB_MEMORY;
	jeli.JobMemoryLimit = peakMemory;
	if( !SetInformationJobObject( h , JobObjectExtendedLimitInformation , &jeli , sizeof( jeli ) ) ) WARN( "Failed to set memory limit" );
}
#else // !_WIN32 && !_WIN64
#include <sys/time.h> 
#include <sys/resource.h> 
inline void SetPeakMemoryMB( size_t sz )
{
	sz <<= 20;
	struct rlimit rl;
	getrlimit( RLIMIT_AS , &rl );
	rl.rlim_cur = sz;
	setrlimit( RLIMIT_AS , &rl );
}
#endif // _WIN32 || _WIN64

/*
* Author:  David Robert Nadeau
* Site:    http://NadeauSoftware.com/
* License: Creative Commons Attribution 3.0 Unported License
*          http://creativecommons.org/licenses/by/3.0/deed.en_US
*/

#if defined(_WIN32) || defined( _WIN64 )
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif





/**
* Returns the peak (maximum so far) resident set size (physical
* memory use) measured in bytes, or zero if the value cannot be
* determined on this OS.
*/
inline size_t getPeakRSS( )
{
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
	/* AIX and Solaris ------------------------------------------ */
	struct psinfo psinfo;
	int fd = -1;
	if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
		return (size_t)0L;      /* Can't open? */
	if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
	{
		close( fd );
		return (size_t)0L;      /* Can't read? */
	}
	close( fd );
	return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* BSD, Linux, and OSX -------------------------------------- */
	struct rusage rusage;
	getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
	return (size_t)rusage.ru_maxrss;
#else
	return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
	/* Unknown OS ----------------------------------------------- */
	return (size_t)0L;          /* Unsupported. */
#endif
}





/**
* Returns the current resident set size (physical memory use) measured
* in bytes, or zero if the value cannot be determined on this OS.
*/
inline size_t getCurrentRSS( )
{
#if defined(_WIN32) || defined( _WIN64 )
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
	/* OSX ------------------------------------------------------ */
	struct mach_task_basic_info info;
	mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
	if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
		(task_info_t)&info, &infoCount ) != KERN_SUCCESS )
		return (size_t)0L;      /* Can't access? */
	return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
	/* Linux ---------------------------------------------------- */
	long rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;      /* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;      /* Can't read? */
	}
	fclose( fp );
	return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
	/* AIX, BSD, Solaris, and Unknown OS ------------------------ */
	return (size_t)0L;          /* Unsupported. */
#endif
}
#endif // MY_MISCELLANY_INCLUDED
