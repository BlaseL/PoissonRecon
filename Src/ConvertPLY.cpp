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

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <limits>
#include <sstream>
#include <unordered_map>
#include "FEMTree.h"
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "MAT.h"
#include "Geometry.h"
#include "Ply.h"
#include "PointStreamData.h"

cmdLineParameter< char* > In( "in" ) , Out( "out" );
cmdLineParameter< float > ChunkWidth( "width" , -1.f );
cmdLineReadable ASCII( "ascii" );

cmdLineReadable* params[] = { &In , &Out , &ChunkWidth , &ASCII , NULL };

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input polygon mesh>\n" , In.name );
	printf( "\t[--%s <ouput polygon mesh name/header>]\n" , Out.name );
	printf( "\t[--%s <chunk width>]\n" , ChunkWidth.name );
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
void GetSubMesh( const std::vector< Vertex > &vertices , const std::vector< std::vector< long long > > &polygons , Point< float , 3 > min , Point< float , 3 > max , std::vector< Vertex > &subVertices , std::vector< std::vector< long long > > &subPolygons )
{
printf("a\n" );
	subVertices.resize( 0 );
	subPolygons.resize( 0 );

printf( "Vertices: %d\n" , (int)vertices.size() );
	for( size_t i=0 ; i<polygons.size() ; i++ )
	{
printf( "processing polygon: %d / %d\n" , (int)i , (int)polygons.size() );
for( int j=0 ; j<polygons[i].size() ; j++ ) printf( " %d" , (int)polygons[i][j] );
printf( "\n" );
		Point< float , 3 > center;
printf( "1\n ");
		for( size_t j=0 ; j<polygons[i].size() ; j++ )
		{
printf( "processing corner %d / %d: %llu\n" , (int)j , (int)polygons[i].size() , (unsigned long long)polygons[i][j] );
			center += vertices[ polygons[i][j] ].point;
		}
printf( "2\n ");
		center /= (float)polygons[i].size();
printf( "Center: %f %f %f\n" , center[0] , center[1] , center[2] );
		bool inside = true;
		for( unsigned int d=0 ; d<3 ; d++ ) if( center[d]<min[d] || center[d]>=max[d] ) inside = false;
		if( inside ) subPolygons.push_back( polygons[i] );
	}
printf("b\n" );

	long long count = 0;
	std::unordered_map< long long , long long > vMap;
	for( size_t i=0 ; i<subPolygons.size() ; i++ ) for( size_t j=0 ; j<subPolygons[i].size() ; j++ )
	{
		auto iter = vMap.find( subPolygons[i][j] );
		if( iter==vMap.end() ) vMap[ subPolygons[i][j] ] = count++;
	}

printf("c\n" );
	subVertices.resize( vMap.size() );
	for( size_t i=0 ; i<subPolygons.size() ; i++ ) for( size_t j=0 ; j<subPolygons[i].size() ; j++ )
	{
		long long oldIdx = subPolygons[i][j];
		long long newIdx = vMap[ oldIdx ];
		subVertices[ newIdx ] = vertices[ oldIdx ];
		subPolygons[i][j] = newIdx;
	}
printf("d\n" );
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


	for( unsigned int d=0 ; d<3 ; d++ ) max[d] += ( max[d]-min[d] ) / 1000.f;

	float width = 0;
	for( unsigned int d=0 ; d<3 ; d++ ) width = std::max< float >( width , max[d]-min[d] );
	if( ChunkWidth.value>=width ) WARN( "width is large enough, outputting mesh as a single file" ) , ChunkWidth.value = -1.f;
	width = ChunkWidth.value;

	if( Out.set )
	{
		if( width>0 )
		{
			float xMin , yMin , zMin;
			int i , j , k;
			for( i=0 , xMin=min[0] ; xMin<=max[0] ; xMin += width , i++ )
				for( j=0 , yMin=min[1] ; yMin<=max[1] ; yMin += width , j++ )
					for( k=0 , zMin=min[2] ; zMin<=max[2] ; zMin += width , k++ )
					{
						std::vector< std::vector< long long > > subPolygons;
						std::vector< Vertex > subVertices;
						GetSubMesh( vertices , polygons , Point< float , 3 >( xMin , yMin , zMin ) , Point< float , 3 >( xMin+width , yMin+width , zMin+width ) , subVertices , subPolygons );
						if( subPolygons.size() )
						{
							printf( "\t[%d %d %d]Vertices / Polygons: %llu / %llu\n" , i , j , k , (unsigned long long)subVertices.size() , (unsigned long long)subPolygons.size() );

							Point< float , 3 > _min , _max;
							GetBoundingBox( subVertices , _min , _max );
							printf( "\t[%d %d %d] " , i , j , k ) ; PrintBoundingBox( _min , _max ) ; printf( "\n" );
							std::stringstream stream;
							stream << Out.value << "." << i << "." << j << "." << k << ".ply";
							WriteMesh( stream.str().c_str() , ASCII.set ? PLY_ASCII : ft , subVertices , subPolygons , comments );
						}
					}
		}
		else WriteMesh( Out.value , ASCII.set ? PLY_ASCII : ft , vertices , polygons , comments );
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
