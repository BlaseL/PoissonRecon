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

#undef VERBOSE_MESSAGING

//////////////////
// OpenMP Stuff //
//////////////////
#ifdef _OPENMP
#include <omp.h>
#else // !_OPENMP
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

#ifdef VERBOSE_MESSAGING
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
#else // !VERBOSE_MESSAGING
	template< typename ... Arguments >
	std::string MakeMessageString( std::string header , std::string functionName , Arguments ... arguments )
	{
		std::stringstream stream;

		// The first line is the header, the file name , and the line number
		stream << header << " " << functionName << ": ";

		_AddToMessageStream( stream , arguments ... );

		return stream.str();
	}

	struct Exception : public std::exception
	{
		const char *what( void ) const noexcept { return _message.c_str(); }
		template< typename ... Args >
		Exception( const char *functionName , const char *format , Args ... args )
		{
			_message = MakeMessageString( "[EXCEPTION]" , functionName , format , args ... );
		}
	private:
		std::string _message;
	};
	template< typename ... Args > void Throw( const char *functionName , const char *format , Args ... args ){ throw Exception( functionName , format , args ... ); }
	template< typename ... Args >
	void Warn( const char *functionName , const char *format , Args ... args )
	{
		std::cerr << MakeMessageString( "[WARNING]" , functionName , format , args ... ) << std::endl;
	}
	template< typename ... Args >
	void ErrorOut( const char *functionName , const char *format , Args ... args )
	{
		std::cerr << MakeMessageString( "[ERROR]" , functionName , format , args ... ) << std::endl;
		exit( 0 );
	}
#endif // VERBOSE_MESSAGING
}
#ifdef VERBOSE_MESSAGING
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
#else // !VERBOSE_MESSAGING
#ifndef WARN
#define WARN( ... ) MKExceptions::Warn( __FUNCTION__ , __VA_ARGS__ )
#endif // WARN
#ifndef WARN_ONCE
#define WARN_ONCE( ... ) { static bool firstTime = true ; if( firstTime ) MKExceptions::Warn( __FUNCTION__ , __VA_ARGS__ ) ; firstTime = false; }
#endif // WARN_ONCE
#ifndef THROW
#define THROW( ... ) MKExceptions::Throw( __FUNCTION__ , __VA_ARGS__ )
#endif // THROW
#ifndef ERROR_OUT
#define ERROR_OUT( ... ) MKExceptions::ErrorOut( __FUNCTION__ , __VA_ARGS__ )
#endif // ERROR_OUT
#endif // VERBOSE_MESSAGING

#ifdef NEW_CODE

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

#if 0
// The assumption is that Kernel is the type of a function taking two arguments, the first is the index of the thread and the second is index of the iteration
struct MKThread
{
	static unsigned int Threads;

#if 1
	template< typename Kernel >
	static void parallel_thread_for( size_t start , size_t end , Kernel kernel )
	{
		size_t range = end - start;

		auto _functor = [&]( unsigned int t )
		{
			size_t _start = start + ( range * t ) / Threads;
			size_t _end = start + ( range *(t+1) ) / Threads;
			for( size_t i=_start ; i<_end ; i++ ) kernel( t , i );
		};

		std::vector< std::thread > threads( Threads );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t].join();
	}

	template< typename Kernel , typename ThreadArg >
	static void parallel_thread_for( size_t start , size_t end , std::vector< ThreadArg > &threadArgs , Kernel kernel )
	{
		size_t range = end - start;

		auto _functor = [&]( unsigned int t )
		{
			size_t _start = start + ( range * t ) / Threads;
			size_t _end = start + ( range *(t+1) ) / Threads;
			ThreadArg &threadArg = threadArgs[t];
			for( size_t i=_start ; i<_end ; i++ ) kernel( t , i , threadArg );
		};

		std::vector< std::thread > threads( Threads );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t].join();
	}

	template< typename Kernel , typename ThreadArg >
	static void parallel_thread_for( size_t start , size_t end , ThreadArg *threadArgs , Kernel kernel )
	{
		size_t range = end - start;

		auto _functor = [&]( unsigned int t )
		{
			size_t _start = start + ( range * t ) / Threads;
			size_t _end = start + ( range *(t+1) ) / Threads;
			ThreadArg &threadArg = threadArgs[t];
			for( size_t i=_start ; i<_end ; i++ ) kernel( t , i , threadArg );
		};

		std::vector< std::thread > threads( Threads );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t].join();
	}

	template< typename Kernel >
	static void parallel_thread_block_for( size_t start , size_t end , size_t blockSize , Kernel kernel )
	{
		size_t range = end - start;
		size_t blocks = ( range + blockSize - 1 ) / blockSize;

		auto _functor = [&]( unsigned int t )
		{
			for( size_t b=t ; b<blocks ; b+=Threads )
			{
				size_t blockStart = start + ( range * b ) / blocks;
				size_t blockEnd = start + ( range * (b+1) ) / blocks;
				for( size_t i=blockStart ; i<blockEnd ; i++ ) kernel( t , i );
			}
		};

		if( blocks==1 ) for( size_t i=start ; i<end ; i++ ) kernel( 0 , i );
		else if( blocks>1 )
		{
			unsigned int threadCount = (unsigned int)std::min< size_t >( Threads , blocks );
			std::vector< std::thread > threads( threadCount );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t].join();
		}
	}

	template< typename Kernel , typename ThreadArg >
	static void parallel_thread_block_for( size_t start , size_t end , size_t blockSize , std::vector< ThreadArg > &threadArgs , Kernel kernel )
	{
		size_t range = end - start;
		size_t blocks = ( range + blockSize - 1 ) / blockSize;

		auto _functor = [&]( unsigned int t )
		{
			for( size_t b=t ; b<blocks ; b+=Threads )
			{
				size_t blockStart = start + ( range * b ) / blocks;
				size_t blockEnd = start + ( range * (b+1) ) / blocks;
				ThreadArg &threadArg = threadArgs[t];
				for( size_t i=blockStart ; i<blockEnd ; i++ ) kernel( t , i , threadArg );
			}
		};

		if( blocks==1 )
		{
			ThreadArg &threadArg = threadArgs[0];
			for( size_t i=start ; i<end ; i++ ) kernel( 0 , i , threadArg );
		}
		else if( blocks>1 )
		{
			unsigned int threadCount = (unsigned int)std::min< size_t >( Threads , blocks );
			std::vector< std::thread > threads( threadCount );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t].join();
		}
	}

	template< typename Kernel , typename ThreadArg >
	static void parallel_thread_block_for( size_t start , size_t end , size_t blockSize , ThreadArg *threadArgs , Kernel kernel )
	{
		size_t range = end - start;
		size_t blocks = ( range + blockSize - 1 ) / blockSize;

		auto _functor = [&]( unsigned int t )
		{
			for( size_t b=t ; b<blocks ; b+=Threads )
			{
				size_t blockStart = start + ( range * b ) / blocks;
				size_t blockEnd = start + ( range * (b+1) ) / blocks;
				ThreadArg &threadArg = threadArgs[t];
				for( size_t i=blockStart ; i<blockEnd ; i++ ) kernel( t , i , threadArg );
			}
		};

		if( blocks==1 )
		{
			ThreadArg &threadArg = threadArgs[0];
			for( size_t i=start ; i<end ; i++ ) kernel( 0 , i , threadArg );
		}
		else if( blocks>1 )
		{
			unsigned int threadCount = (unsigned int)std::min< size_t >( Threads , blocks );
			std::vector< std::thread > threads( threadCount );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t].join();
		}
	}

	template< typename Kernel >
	static void parallel_task_block_for( size_t start , size_t end , size_t blockSize , Kernel &kernel )
	{
		std::atomic< size_t > idx = start;
		auto _functor = [&]( unsigned int t )
		{
			while( idx<end )
			{
				size_t _start = idx.fetch_add( blockSize );
				size_t _end = std::min< size_t >( _start + blockSize , end );
				for( size_t i=_start ; i<_end ; i++ ) kernel( t , i );
			}
		};

		if( end-start<=blockSize ) for( size_t i=start ; i<end ; i++ ) kernel( 0 , i ); 
		else
		{
			unsigned int threadCount = (unsigned int)std::min< size_t >( Threads , (end-start+blockSize-1)/blockSize );
			std::vector< std::thread > threads( threadCount );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t].join();
		}
	}

	template< typename Kernel , typename ThreadArg >
	static void parallel_task_block_for( size_t start , size_t end , size_t blockSize , std::vector< ThreadArg > &threadArgs , Kernel &kernel )
	{
		std::atomic< size_t > idx = start;
		auto _functor = [&]( unsigned int t )
		{
			while( idx<end )
			{
				size_t _start = idx.fetch_add( blockSize );
				size_t _end = std::min< size_t >( _start + blockSize , end );
				ThreadArg &threadArg = threadArgs[t];
				for( size_t i=_start ; i<_end ; i++ ) kernel( t , i , threadArg );
			}
		};

		if( end-start<=blockSize )
		{
			ThreadArg &threadArg = threadArgs[0];
			for( size_t i=start ; i<end ; i++ ) kernel( 0 , i , threadArg ); 
		}
		else
		{
			unsigned int threadCount = (unsigned int)std::min< size_t >( Threads , (end-start+blockSize-1)/blockSize );
			std::vector< std::thread > threads( threadCount );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t].join();
		}
	}

	template< typename Kernel , typename ThreadArg >
	static void parallel_task_block_for( size_t start , size_t end , size_t blockSize , ThreadArg *threadArgs , Kernel &kernel )
	{
		std::atomic< size_t > idx = start;
		auto _functor = [&]( unsigned int t )
		{
			while( idx<end )
			{
				size_t _start = idx.fetch_add( blockSize );
				size_t _end = std::min< size_t >( _start + blockSize , end );
				ThreadArg &threadArg = threadArgs[t];
				for( size_t i=_start ; i<_end ; i++ ) kernel( t , i , threadArg );
			}
		};

		if( end-start<=blockSize )
		{
			ThreadArg &threadArg = threadArgs[0];
			for( size_t i=start ; i<end ; i++ ) kernel( 0 , i , threadArg ); 
		}
		else
		{
			unsigned int threadCount = (unsigned int)std::min< size_t >( Threads , (end-start+blockSize-1)/blockSize );
			std::vector< std::thread > threads( threadCount );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
			for( unsigned int t=0 ; t<threadCount ; t++ ) threads[t].join();
		}
	}

	template< typename Kernel >
	static void parallel_for( size_t start , size_t end , Kernel kernel )
	{
//		return parallel_thread_for( start , end , kernel );
		return parallel_thread_block_for( start , end , 200 , kernel );
//		return parallel_task_block_for( start , end , 200 , kernel );
	}

	template< typename Kernel , typename ThreadArg >
	static void parallel_for( size_t start , size_t end , std::vector< ThreadArg > &threadArgs , Kernel kernel )
	{
//		return parallel_thread_for( start , end , threadArgs , kernel );
		return parallel_thread_block_for( start , end , 200 , threadArgs , kernel );
//		return parallel_task_block_for( start , end , 200 , threadArgs , kernel );
	}

	template< typename Kernel , typename ThreadArg >
	static void parallel_for( size_t start , size_t end , ThreadArg *threadArgs , Kernel kernel )
	{
//		return parallel_thread_for( start , end , threadArgs , kernel );
		return parallel_thread_block_for( start , end , 200 , threadArgs , kernel);
//		return parallel_task_block_for( start , end , 200 , threadArgs , kernel );
	}
#else
	template< typename Kernel , typename ... Args >
	static void parallel_thread_for( size_t start , size_t end , Kernel kernel , Args ... args )
	{
		size_t range = end - start;

		auto _functor = [&]( unsigned int t )
		{
			size_t _start = start + ( range * t ) / Threads;
			size_t _end = start + ( range *(t+1) ) / Threads;
			for( size_t i=_start ; i<_end ; i++ ) kernel( t , i , args ... );
		};

		std::vector< std::thread > threads( Threads );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t].join();
	}

	template< typename Kernel , typename ... Args >
	static void parallel_thread_block_for( size_t start , size_t end , size_t blockSize , Kernel kernel , Args ... args )
	{
		size_t range = end - start;
		size_t blocks = ( range + blockSize - 1 ) / blockSize;
		auto _functor = [&]( unsigned int t )
		{
			for( size_t b=t ; b<blocks ; b+=Threads )
			{
				size_t blockStart = start + ( range * b ) / blocks;
				size_t blockEnd = start + ( range * (b+1) ) / blocks;
				for( size_t i=blockStart ; i<blockEnd ; i++ ) kernel( t , i , args ... );
			}
		};

		std::vector< std::thread > threads( Threads );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t].join();
	}

	template< typename Kernel , typename ... Args >
	static void parallel_task_block_for( size_t start , size_t end , size_t blockSize , Kernel &kernel , Args ... args )
	{
		std::atomic< size_t > idx = start;
		auto _functor = [&]( unsigned int t )
		{
			while( idx<end )
			{
				size_t _start = idx.fetch_add( blockSize );
				size_t _end = std::min< size_t >( _start + blockSize , end );
				for( size_t i=_start ; i<_end ; i++ ) kernel( t , i , args ... );
			}
		};

		std::vector< std::thread > threads( Threads );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );
		for( unsigned int t=0 ; t<Threads ; t++ ) threads[t].join();
	}

	template< typename Kernel , typename ... Args >
	static void parallel_for( size_t start , size_t end , Kernel kernel , Args ... args )
	{
		//		return parallel_thread_for( start , end , kernel , args ... );
		return parallel_thread_block_for( start , end , 100 , kernel , args ... );
		//		return parallel_task_block_for( start , end , 200 , kernel , args ... );
	}
#endif

};
unsigned int MKThread::Threads = std::thread::hardware_concurrency();
#endif

#define NEW_THREAD_NUM

#ifdef USE_HH_THREAD_POOL
#include <memory>
inline int get_max_threads( void ) { return std::max< int >( (int)std::thread::hardware_concurrency() , 1 ); }

class ThreadPoolIndexedTask
{
public:
	using Task = std::function< void (int) >;
	ThreadPoolIndexedTask( void )
	{
		const int num_threads = get_max_threads();
		_threads.reserve( num_threads );
		for( int i=0 ; i<num_threads ; i++ ) _threads.emplace_back( &ThreadPoolIndexedTask::worker_main , this );
	}
	~ThreadPoolIndexedTask( void )
	{
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
	int num_threads( void ) const { return (int)_threads.size(); }
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


// Evaluates function(element) for each element in range by parallelizing across chunks of elements using
// a cached thread pool.  The range must support begin/end functions returning random-access iterators.
// Parallelism is disabled if the estimated cost (estimated_cycles_per_element * size(range)) is less than some
// internal threshold, or if this parallel_for_each() is already executing (nested) within another parallel_for_each().
// Exceptions within function() cause program termination as they are not caught.
// One drawback over OpenMP is that if an exception or abort occurs within function(), the stack trace will not
// include the functions that called parallel_for_each() because these lie in the stack frames of a different thread.
struct ThreadPool
{
	struct ThreadNum
	{
#ifdef FORCE_OMP
		unsigned int operator()( void ) const { return (unsigned int)omp_get_thread_num(); }
#else // !FORCE_OMP
		ThreadNum( unsigned int threadNum ) : _threadNum( threadNum ){}
		unsigned int operator()( void ) const { return _threadNum; }
	protected:
		unsigned int _threadNum;
#endif // FORCE_OMP
	};

	static size_t DefaultThreadNum;
	static size_t DefaultMaxBlocksPerThread;
	static const uint64_t k_omp_thresh = 100*1000; // number of instruction cycles above which a loop should be parallelized
	static constexpr uint64_t k_parallelism_always = k_omp_thresh;

	unsigned int threadNum( void ) const { return get_max_threads(); }

	template< typename Function >
	void parallel_for( size_t begin , size_t end , const Function &function , uint64_t estimated_cycles_per_element=k_parallelism_always )
	{
#ifdef FORCE_OMP
		ThreadNum threadNum;
#pragma omp parallel for
		for( long long i=(long long)begin ; i<(long long)end ; i++ ) function( threadNum , i );
#else // !FORCE_OMP
		const size_t num_elements = size_t( end - begin );
		uint64_t total_num_cycles = num_elements * estimated_cycles_per_element;
		const int max_num_threads = get_max_threads();
		const int num_threads = (int)std::min< size_t >( max_num_threads , num_elements );
		const bool desire_parallelism = num_threads > 1 && total_num_cycles >= k_omp_thresh;
		ThreadPoolIndexedTask* const thread_pool = desire_parallelism ? &ThreadPoolIndexedTask::default_threadpool() : nullptr;
		if( !thread_pool || thread_pool->already_active() )
		{
			// Traverse the range elements sequentially.
			ThreadNum threadNum(0);
			for( size_t index=begin ; index<end ; ++index) function( threadNum , index );
		}
		else
		{
			// Traverse the range elements in parallel.
			const size_t chunk_size = ( num_elements + num_threads - 1 ) / num_threads;
			thread_pool->execute( num_threads , [ begin , num_elements , chunk_size , &function ]( int thread_index )
			{
				size_t index_start = thread_index * chunk_size;
				size_t index_stop = begin + std::min< size_t >( ( thread_index + 1 ) * chunk_size , num_elements );
				for( size_t index = begin + index_start ; index< index_stop ; ++index ) function( thread_index , index );
			}
			);
		}
#endif // FORCE_OMP
	}
};
size_t ThreadPool::DefaultThreadNum;
size_t ThreadPool::DefaultMaxBlocksPerThread;

#else // !USE_HH_THREAD_POOL

// [WARNING] The ThreadPool is not thread safe. Specifically, every thread should have its own ThreadPool.
struct ThreadPool
{
	struct ThreadNum
	{
#ifdef FORCE_OMP
		unsigned int operator()( void ) const { return (unsigned int)omp_get_thread_num(); }
#else // !FORCE_OMP
		ThreadNum( unsigned int threadNum ) : _threadNum( threadNum ){}
		unsigned int operator()( void ) const { return _threadNum; }
	protected:
		unsigned int _threadNum;
#endif // FORCE_OMP
	};

	static unsigned int DefaultThreadNum;
#ifdef NEW_THREAD_POOL
	static size_t DefaultMaxBlocksPerThread;
	static size_t DefaultMinBlockSize;
#else // !NEW_THREAD_POOL
	static size_t DefaultBlockSize;
	static size_t DefaultMinParallelSize;
#endif // NEW_THREAD_POOL

#ifdef NEW_THREAD_POOL
#ifdef NEW_THREAD_NUM
	void parallel_for( size_t begin , size_t end , std::function< void ( const ThreadNum & , size_t ) > iterationFunction , size_t maxBlocksPerThread=DefaultMaxBlocksPerThread , size_t minBlockSize=DefaultMinBlockSize )
#else // !NEW_THREAD_NUM
	void parallel_for( size_t begin , size_t end , std::function< void ( unsigned int , size_t ) > iterationFunction , size_t maxBlocksPerThread=DefaultMaxBlocksPerThread , size_t minBlockSize=DefaultMinBlockSize )
#endif // NEW_THREAD_NUM
#else // !NEW_THREAD_POOL
#ifdef NEW_THREAD_NUM
	void parallel_for( size_t begin , size_t end , std::function< void ( const ThreadNum & , size_t ) > iterationFunction , size_t blockSize=DefaultBlockSize , size_t minParallelSize=DefaultMinParallelSize )
#else // !NEW_THREAD_NUM
	void parallel_for( size_t begin , size_t end , std::function< void ( unsigned int , size_t ) > iterationFunction , size_t blockSize=DefaultBlockSize , size_t minParallelSize=DefaultMinParallelSize )
#endif //NEW_THREAD_NUM
#endif // NEW_THREAD_POOL
	{
#ifdef FORCE_OMP
		ThreadNum threadNum;
#pragma omp parallel for num_threads( (int)_threads.size() )
		for( long long i=(long long)begin ; i<(long long)end ; i++ ) iterationFunction( threadNum , i );
#else // !FORCE_OMP
#ifdef NEW_THREAD_POOL
		size_t blocks = _threads.size() * maxBlocksPerThread;
		_blockSize = std::max< size_t >( minBlockSize , ( end - begin ) / blocks );
#else // !NEW_THREAD_POOL
		_blockSize = blockSize;
#endif // NEW_THREAD_POOL
		if( !_threads.size() ) for( size_t i=begin ; i<end ; i++ ) iterationFunction( 0 , i );
		else
		{
			auto WorkersWorking = [&]( void ){ unsigned int c=0 ; for( unsigned int t=0 ; t<_threads.size() ; t++ ) if( _workToBeDone[t] ) c++ ; return c; };
			auto WorkComplete = [&]( void ){ return WorkersWorking()==0; };
			if( begin>=end ) return;
#ifdef NEW_THREAD_POOL
			if( end-begin<_blockSize && _blockSize!=-1 ) for( size_t i=begin ; i<end ; i++ ) iterationFunction( 0 , i );
#else // !NEW_THREAD_POOL
			if( end-begin<minParallelSize && _blockSize!=-1 ) for( size_t i=begin ; i<end ; i++ ) iterationFunction( 0 , i );
#endif // NEW_THREAD_POOL
			else
			{
				_begin = begin;
				_index = begin;
				_end = end;
				_iterationFunction = iterationFunction;
				for( unsigned int t=0 ; t<_threads.size() ; t++ ) _workToBeDone[t] = 1;
				_waitingForWorkOrClose.notify_all();

				{
					std::unique_lock< std::mutex > lock( _mutex );
//					while( !WorkComplete()!=0 ){ _doneWithWork.wait( lock , WorkComplete ); }
					_doneWithWork.wait( lock , WorkComplete );
				}
			}
		}
#endif // FORCE_OMP
	}
	void setThreadNum( unsigned int threadNum )
	{
#ifdef FORCE_OMP
		_threads.resize( threadNum );
#else // !FORCE_OMP
		_close = 1;
		if( _threads.size() )
		{
			_waitingForWorkOrClose.notify_all();
			for( unsigned int t=0 ; t<_threads.size() ; t++ ) _threads[t].join();
		}
		_threads.resize( threadNum );
		_workToBeDone.resize( threadNum , 0 );
		_close = 0;
		for( unsigned int t=0 ; t<threadNum ; t++ ) _threads[t] = std::thread( _ThreadFunction , t , this );
#endif // FORCE_OMP
	}

	unsigned int threadNum( void ) const { return (unsigned int)_threads.size(); }

	ThreadPool( unsigned int threadNum=DefaultThreadNum )
	{
		_threads.resize( threadNum );
		_workToBeDone.resize( threadNum , 0 );
		_close = 0;
#ifdef FORCE_OMP
#else // !FORCE_OMP
		for( unsigned int t=0 ; t<threadNum ; t++ ) _threads[t] = std::thread( _ThreadFunction , t , this );
#endif // FORCE_OMP
	}

	~ThreadPool( void )
	{
#ifdef FORCE_OMP
#else // !FORCE_OMP
		_close = 1;
		if( _threads.size() )
		{
			_waitingForWorkOrClose.notify_all();
			for( unsigned int t=0 ; t<_threads.size() ; t++ ) _threads[t].join();
		}
#endif // FORCE_OMP
	}
protected:
	ThreadPool( const ThreadPool & ){}
	ThreadPool &operator = ( const ThreadPool & ){}
#ifdef FORCE_OMP
#else // !FORCE_OMP
	static void _ThreadFunction( unsigned int thread , ThreadPool *tPool )
	{
		unsigned int threads = (unsigned int)tPool->_threads.size();
#ifdef NEW_THREAD_NUM
		ThreadNum threadNum( thread );
#endif // NEW_THREAD_NUM

		// Wait for the first job to come in
		{
			std::unique_lock< std::mutex > lock( tPool->_mutex );
			tPool->_waitingForWorkOrClose.wait( lock , [&]{ return tPool->_workToBeDone[thread] || tPool->_close; } );
		}
		while( true )
		{
			// do the job (or terminate)
			if( tPool->_close ) break;
			if( tPool->_blockSize==-1 )
			{
				size_t range = ( tPool->_end - tPool->_begin );
				size_t begin = tPool->_begin + ( range * (thread+0) )/threads;
				size_t end   = tPool->_begin + ( range * (thread+1) )/threads;
#ifdef NEW_THREAD_NUM
				for( size_t i=begin ; i<end ; i++ ) tPool->_iterationFunction( threadNum , i );
#else // !NEW_THREAD_NUM
				for( size_t i=begin ; i<end ; i++ ) tPool->_iterationFunction( thread , i );
#endif // NEW_THREAD_NUM
			}
			else
			{
				while( tPool->_index<tPool->_end )
				{
					size_t begin = tPool->_index.fetch_add( tPool->_blockSize );
					size_t end = std::min< size_t >( begin + tPool->_blockSize , tPool->_end );
#ifdef NEW_THREAD_NUM
					for( size_t i=begin ; i<end ; i++ ) tPool->_iterationFunction( threadNum , i );
#else // !NEW_THREAD_NUM
					for( size_t i=begin ; i<end ; i++ ) tPool->_iterationFunction( thread , i );
#endif // NEW_THREAD_NUM
				}
			}

			// Notify and wait for the next job
			// [NOTE] We are locking the mutex so that the server cannot initiate new jobs until the client has started waiting
			{
				std::unique_lock< std::mutex > lock( tPool->_mutex );
				tPool->_workToBeDone[thread] = false;
				tPool->_doneWithWork.notify_all();
				tPool->_waitingForWorkOrClose.wait( lock , [&]{ return tPool->_workToBeDone[thread] || tPool->_close; } );
			}
		}
	}
#endif // FORCE_OMP

	char _close;
	std::vector< bool > _workToBeDone;
	std::mutex _mutex;
	std::condition_variable _waitingForWorkOrClose , _doneWithWork;
	std::vector< std::thread > _threads;
#ifdef NEW_THREAD_NUM
	std::function< void ( const ThreadNum & , size_t ) > _iterationFunction;
#else // !NEW_THREAD_NUM
	std::function< void ( unsigned int , size_t ) > _iterationFunction;
#endif // NEW_THREAD_NUM
	size_t _begin , _end , _blockSize;
	std::atomic< size_t > _index;
};
unsigned int ThreadPool::DefaultThreadNum = std::thread::hardware_concurrency();
#ifdef NEW_THREAD_POOL
size_t ThreadPool::DefaultMaxBlocksPerThread=100;
size_t ThreadPool::DefaultMinBlockSize=100;
#else // !NEW_THREAD_POOL
size_t ThreadPool::DefaultBlockSize = 100;
size_t ThreadPool::DefaultMinParallelSize = 100;
#endif // NEW_THREAD_POOL
#define PARALLEL_FOR( threadPool , begin , end , functor ) if( (begin)<(end) ) threadPool.parallel_for( (begin) , (end) , functor )
#endif // USE_HH_THREAD_POOL

#endif // NEW_THREADS

#ifdef NEW_CODE
#ifdef NEW_THREADS
#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#endif // _WIN32 || _WIN64
#endif // NEW_THREADS

template< typename Number >
void AddAtomic32( Number &a , const Number &b )
{
#ifdef NEW_THREADS
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
#else // !NEW_THREADS
#pragma omp atomic
	a += b;
#endif // NEW_THREADS
}

template< typename Number >
void AddAtomic64( Number &a , const Number &b )
{
#ifdef NEW_THREADS
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
#else // !NEW_THREADS
#pragma omp atomic
	a += b;
#endif // NEW_THREADS
}

template< typename Data >
void AddAtomic( Data& a , const Data& b )
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


#if 0
void AddAtomic( float& a , const float& b )
{
#ifdef NEW_THREADS
#if defined( _WIN32 ) || defined( _WIN64 )
	float current = a;
	float sum = current+b;
	long &_current = *(long *)&current;
	long &_sum = *(long *)&sum;
	while( InterlockedCompareExchange( (long*)&a , _sum , _current )!=_current ) current = a , sum = a+b;
#else // !_WIN32 && !_WIN64
	float current = a;
	float sum = current+b;
	uint32_t &_current = *(uint32_t *)&current;
	uint32_t &_sum = *(uint32_t *)&sum;
	while( __sync_val_compare_and_swap( (uint32_t *)&a , _current , _sum )!=_current ) current = a , sum = a+b;
#endif // _WIN32 || _WIN64
#else // !NEW_THREADS
#pragma omp atomic
	a += b;
#endif // NEW_THREADS
}

void AddAtomic( double& a , const double& b )
{
#ifdef NEW_THREADS
#if defined( _WIN32 ) || defined( _WIN64 )
	double current = a;
	double sum = current+b;
	__int64 &_current = *(__int64 *)&current;
	__int64 &_sum = *(__int64 *)&sum;
	while( InterlockedCompareExchange64( (__int64*)&a , _sum , _current )!=_current ) current = a , sum = a+b;
#else // !_WIN32 && !_WIN64
	double current = a;
	double sum = current+b;
	uint64_t &_current = *(uint64_t *)&current;
	uint64_t &_sum = *(uint64_t *)&sum;
	while( __sync_val_compare_and_swap( (uint64_t*)&a , _current , _sum )!=_current ) current = a , sum = a+b;
#endif // _WIN32 || _WIN64
#else // !NEW_THREADS
#pragma omp atomic
	a += b;
#endif // NEW_THREADS
}
#endif // NEW_CODE

#endif
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
