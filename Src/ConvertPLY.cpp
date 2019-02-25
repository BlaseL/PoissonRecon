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
DAMAGE.
*/
#define NEW_CODE
#define NEW_THREADS

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <limits>
#include <sstream>
#include <unordered_map>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "Geometry.h"
#include "Ply.h"
#include "PointStream.h"
#include "PointStreamData.h"

cmdLineParameter< char* > In( "in" ) , Out( "out" );
cmdLineParameter< float > Width( "width" , -1.f ) , PadRadius( "radius" , 0.f );
cmdLineReadable ASCII( "ascii" );

cmdLineReadable* params[] = { &In , &Out , &Width , &PadRadius , &ASCII , NULL };

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input polygon mesh>\n" , In.name );
	printf( "\t[--%s <ouput polygon mesh name/header>]\n" , Out.name );
	printf( "\t[--%s <chunk width>=%f]\n" , Width.name , Width.value );
	printf( "\t[--%s <padding radius>=%f]\n" , PadRadius.name , PadRadius.value );
	printf( "\t[--%s]\n" , ASCII.name );
}

void PrintBoundingBox( Point< float , 3 > min , Point< float , 3 > max )
{
	printf( "[" );
	for( unsigned int d=0 ; d<3 ; d++ ) printf( " %f" , min[d] );
	printf( " ] [" );
	for( unsigned int d=0 ; d<3 ; d++ ) printf( " %f" , max[d] );
	printf( " ]" );
}

template< typename ... VertexData >
void WritePoints( const char *fileName , int ft , const std::vector< PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > > &vertices , const std::vector< std::string > &comments )
{
	typedef PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > Vertex;
	typedef MultiPointStreamData< float , VertexData ... > StreamData;

	PLYOutputPointStreamWithData< float , 3 , StreamData > pointStream( fileName , vertices.size() , ft , StreamData::PlyWriteProperties() , StreamData::PlyReadNum );
	for( size_t i=0 ; i<vertices.size() ; i++ ) pointStream.nextPoint( vertices[i].point , vertices[i].data );
}


template< typename ... VertexData >
void WriteMesh( const char *fileName , int ft , const std::vector< PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > > &vertices , const std::vector< std::vector< long long > > &polygons , const std::vector< std::string > &comments )
{
	typedef PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > Vertex;

	if( vertices.size()>std::numeric_limits< int >::max() )
	{
		if( vertices.size()>std::numeric_limits< unsigned int >::max() ) ERROR_OUT( "more vertices than can be indexed by an unsigned int: %llu" , (unsigned long long)vertices.size() );
		WARN( "more vertices than can be indexed by an int, using unsigned int instead: %llu" , (unsigned long long)vertices.size() );
		std::vector< std::vector< unsigned int > > outPolygons;
		outPolygons.resize( polygons.size() );
		for( size_t i=0 ; i<polygons.size() ; i++ )
		{
			outPolygons[i].resize( polygons[i].size() );
			for( int j=0 ; j<polygons[i].size() ; j++ ) outPolygons[i][j] = (unsigned int)polygons[i][j];
		}
		if( !PlyWritePolygons< Vertex >( fileName , vertices , outPolygons , Vertex::PlyWriteProperties() , Vertex::PlyWriteNum , ft , comments ) ) ERROR_OUT( "Could not write mesh to: " , fileName );
	}
	else
	{
		std::vector< std::vector< int > > outPolygons;
		outPolygons.resize( polygons.size() );
		for( size_t i=0 ; i<polygons.size() ; i++ )
		{
			outPolygons[i].resize( polygons[i].size() );
			for( int j=0 ; j<polygons[i].size() ; j++ ) outPolygons[i][j] = (int)polygons[i][j];
		}
		if( !PlyWritePolygons< Vertex >( fileName , vertices , outPolygons , Vertex::PlyWriteProperties() , Vertex::PlyWriteNum , ft , comments ) ) ERROR_OUT( "Could not write mesh to: " , fileName );
	}
}

template< typename Vertex >
void GetBoundingBox( const std::vector< Vertex > &vertices , Point< float , 3 > &min , Point< float , 3 > &max )
{
	min = max = vertices[0].point;
	for( size_t i=0 ; i<vertices.size() ; i++ ) for( unsigned int d=0 ; d<3 ; d++ )
	{
		min[d] = std::min< float >( min[d] , vertices[i].point[d] );
		max[d] = std::max< float >( max[d] , vertices[i].point[d] );
	}
}

template< typename Vertex >
void GetSubPoints( const std::vector< Vertex > &vertices , Point< float , 3 > min , Point< float , 3 > max , std::vector< Vertex > &subVertices )
{
	subVertices.resize( 0 );

	for( size_t i=0 ; i<vertices.size() ; i++ )
	{
		bool inside = true;
		for( unsigned int d=0 ; d<3 ; d++ ) if( vertices[i].point[d]<min[d] || vertices[i].point[d]>=max[d] ) inside = false;
		if( inside ) subVertices.push_back( vertices[i] );
	}
}

template< typename Vertex >
void GetSubMesh( const std::vector< Vertex > &vertices , const std::vector< std::vector< long long > > &polygons , Point< float , 3 > min , Point< float , 3 > max , std::vector< Vertex > &subVertices , std::vector< std::vector< long long > > &subPolygons )
{
	subVertices.resize( 0 );
	subPolygons.resize( 0 );

	for( size_t i=0 ; i<polygons.size() ; i++ )
	{
		Point< float , 3 > center;
		for( size_t j=0 ; j<polygons[i].size() ; j++ ) center += vertices[ polygons[i][j] ].point;
		center /= (float)polygons[i].size();
		bool inside = true;
		for( unsigned int d=0 ; d<3 ; d++ ) if( center[d]<min[d] || center[d]>=max[d] ) inside = false;
		if( inside ) subPolygons.push_back( polygons[i] );
	}

	long long count = 0;
	std::unordered_map< long long , long long > vMap;
	for( size_t i=0 ; i<subPolygons.size() ; i++ ) for( size_t j=0 ; j<subPolygons[i].size() ; j++ )
	{
		auto iter = vMap.find( subPolygons[i][j] );
		if( iter==vMap.end() ) vMap[ subPolygons[i][j] ] = count++;
	}

	subVertices.resize( vMap.size() );
	for( size_t i=0 ; i<subPolygons.size() ; i++ ) for( size_t j=0 ; j<subPolygons[i].size() ; j++ )
	{
		long long oldIdx = subPolygons[i][j];
		long long newIdx = vMap[ oldIdx ];
		subVertices[ newIdx ] = vertices[ oldIdx ];
		subPolygons[i][j] = newIdx;
	}
}

template< typename ... VertexData >
int Execute( void )
{
	typedef PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > Vertex;
	std::vector< Vertex > vertices;
	std::vector< std::vector< long long > > polygons;

	int ft;
	std::vector< std::string > comments;
	PlyReadPolygons< Vertex >( In.value , vertices , polygons , Vertex::PlyReadProperties() , Vertex::PlyReadNum , ft , comments );
	printf( "Vertices / Polygons: %llu / %llu\n" , (unsigned long long)vertices.size() , (unsigned long long)polygons.size() );

	Point< float , 3 > min , max;
	GetBoundingBox( vertices , min , max );
	PrintBoundingBox( min , max ) ; printf( "\n" );

	float width = Width.value;

	if( Out.set )
	{
		if( width>0 )
		{
			size_t vCount=0 , pCount=0;
			for( unsigned int d=0 ; d<3 ; d++ ) min[d] -= width/10000.f , max[d] += width/10000.f;
			int begin[] = { (int)floor( min[0]/width ) , (int)floor( min[1]/width ) , (int)floor( min[2]/width ) };
			int end  [] = { (int)ceil ( max[0]/width ) , (int)ceil ( max[1]/width ) , (int)ceil ( max[2]/width ) };
			for( int i=begin[0] ; i<end[0] ; i++ ) for( int j=begin[1] ; j<end[1] ; j++ ) for( int k=begin[2] ; k<end[2] ; k++ )
			{
				std::vector< std::vector< long long > > subPolygons;
				std::vector< Vertex > subVertices;
				Point< float , 3 > _min( (i+0)*width - PadRadius.value , (j+0)*width - PadRadius.value , (k+0)*width - PadRadius.value );
				Point< float , 3 > _max( (i+1)*width + PadRadius.value , (j+1)*width + PadRadius.value , (k+1)*width + PadRadius.value );
				if( polygons.size() )
				{
					GetSubMesh( vertices , polygons , _min , _max , subVertices , subPolygons );
					if( subPolygons.size() )
					{
						std::stringstream stream;
						stream << Out.value << "." << i << "." << j << "." << k << ".ply";

						printf( "\t%s\n" , stream.str().c_str() );
						printf( "\t\tVertices / Polygons: %llu / %llu\n" , (unsigned long long)subVertices.size() , (unsigned long long)subPolygons.size() );
						printf( "\t\t" ) ; PrintBoundingBox( _min , _max ) ; printf( "\n" );

#if 1
						WARN( "Not writing mesh: " stream.str() );
#else
						WriteMesh( stream.str().c_str() , ASCII.set ? PLY_ASCII : ft , subVertices , subPolygons , comments );
#endif

						vCount += subVertices.size() , pCount += subPolygons.size();
					}
				}
				else
				{
					GetSubPoints( vertices , _min , _max , subVertices );
					if( subVertices.size() )
					{
						std::stringstream stream;
						stream << Out.value << "." << i << "." << j << "." << k << ".ply";

						printf( "\t%s\n" , stream.str().c_str() );
						printf( "\t\tPoints: %llu\n" , (unsigned long long)subVertices.size() );
						printf( "\t\t" ) ; PrintBoundingBox( _min , _max ) ; printf( "\n" );

#if 1
						WARN( "Not writing points: " stream.str() );
#else
						WritePoints( stream.str().c_str() , ASCII.set ? PLY_ASCII : ft , subVertices , comments );
#endif

						vCount += subVertices.size();
					}
				}
			}
			if( vCount!=vertices.size() ) WARN( "vertex counts don't match" , vertices.size() , vCount );
			if( pCount!=polygons.size() ) WARN( "polygon counts don't match" , polygons.size() , pCount );
		}
		else
		{
#if 1
			WARN( "Not writing: " , Out.value );
#else
			if( polygons.size() ) WriteMesh( Out.value , ASCII.set ? PLY_ASCII : ft , vertices , polygons , comments );
			else WritePoints( Out.value , ASCII.set ? PLY_ASCII : ft , vertices , comments );
#endif
		}
	}

	return EXIT_SUCCESS;
}
int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , &argv[1] , params );

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	typedef MultiPointStreamData< float , PointStreamValue< float > , PointStreamNormal< float , 3 > , PointStreamColor< float > > VertexData;
	typedef PlyVertexWithData< float , 3 , VertexData > Vertex;
	bool readFlags[ Vertex::PlyReadNum ];
	if( !PlyReadHeader( In.value , Vertex::PlyReadProperties() , Vertex::PlyReadNum , readFlags ) ) ERROR_OUT( "Failed to read ply header: " , In.value );

	bool hasValue  = VertexData::ValidPlyReadProperties< 0 >( readFlags + 3 );
	bool hasNormal = VertexData::ValidPlyReadProperties< 1 >( readFlags + 3 );
	bool hasColor  = VertexData::ValidPlyReadProperties< 2 >( readFlags + 3 );

	if( hasValue )
		if( hasColor )
			if( hasNormal ) return Execute< PointStreamValue< float > , PointStreamNormal< float , 3 > , PointStreamColor< float > >();
			else            return Execute< PointStreamValue< float > ,                                  PointStreamColor< float > >();
		else
			if( hasNormal ) return Execute< PointStreamValue< float > , PointStreamNormal< float , 3 >                             >();
			else            return Execute< PointStreamValue< float >                                                              >();
	else
		if( hasColor )
			if( hasNormal ) return Execute< PointStreamNormal< float , 3 > , PointStreamColor< float > >();
			else            return Execute<                                  PointStreamColor< float > >();
		else
			if( hasNormal ) return Execute< PointStreamNormal< float , 3 >                             >();
			else            return Execute<                                                            >();
}
