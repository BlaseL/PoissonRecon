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

#define DIMENSION 3

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include "FEMTree.h"
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "MAT.h"
#include "Geometry.h"
#include "Ply.h"
#include "PointStreamData.h"

cmdLineParameter< char* > In( "in" ) , Out( "out" );
cmdLineParameterArray< float , 2*DIMENSION > BoundingBox( "bBox" );
cmdLineReadable ASCII( "ascii" );

cmdLineReadable* params[] = { &In , &Out , &BoundingBox , &ASCII , NULL };

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input polygon mesh>\n" , In.name );
	printf( "\t[--%s <ouput polygon mesh>]\n" , Out.name );
	printf( "\t[--%s <bounding box>]\n" , BoundingBox.name );
	printf( "\t[--%s]\n" , ASCII.name );
}

template< typename ... VertexData >
void WriteMesh( const char *fileName , int ft , const std::vector< PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , VertexData ... > > > &vertices , const std::vector< std::vector< long long > > &polygons , const std::vector< std::string > &comments )
{
	typedef PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , VertexData ... > > Vertex;
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

template< typename ... VertexData >
int Execute( void )
{
	typedef PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , VertexData ... > > Vertex;
	std::vector< Vertex > vertices;
	std::vector< std::vector< long long > > polygons;

	int ft;
	std::vector< std::string > comments;
	PlyReadPolygons< Vertex >( In.value , vertices , polygons , Vertex::PlyReadProperties() , Vertex::PlyReadNum , ft , comments );
	printf( "Vertices / Polygons: %llu / %llu\n" , (unsigned long long)vertices.size() , (unsigned long long)polygons.size() );

	Point< float , DIMENSION > min , max;
	min = max = vertices[0].point;
	for( size_t i=0 ; i<vertices.size() ; i++ ) for( unsigned int d=0 ; d<DIMENSION ; d++ )
	{
		min[d] = std::min< float >( min[d] , vertices[i].point[d] );
		max[d] = std::max< float >( max[d] , vertices[i].point[d] );
	}
	printf( "min / max: [" );
	for( unsigned int d=0 ; d<DIMENSION ; d++ ) printf( " %f" , min[d] );
	printf( " ] [" );
	for( unsigned int d=0 ; d<DIMENSION ; d++ ) printf( " %f" , max[d] );
	printf( " ]\n" );

	if( BoundingBox.set )
	{
		std::vector< std::vector< long long > > subPolygons;
		for( size_t i=0 ; i<polygons.size() ; i++ )
		{
			Point< float , DIMENSION > center;
			for( size_t j=0 ; j<polygons[i].size() ; j++ ) center += vertices[ polygons[i][j] ].point;
			center /= (float)polygons[i].size();
			bool inside = true;
			for( unsigned int d=0 ; d<DIMENSION ; d++ ) if( center[d]<BoundingBox.values[d] || center[d]>=BoundingBox.values[d+DIMENSION] ) inside = false;
			if( inside ) subPolygons.push_back( polygons[i] );
		}

		long long count = 0;
		std::unordered_map< long long , long long > vMap;
		for( size_t i=0 ; i<polygons.size() ; i++ ) for( size_t j=0 ; j<polygons[i].size() ; j++ )
		{
			auto iter = vMap.find( polygons[i][j] );
			if( iter==vMap.end() ) vMap[ polygons[i][j] ] = count++;
		}
		std::vector< Vertex > subVertices( vMap.size() );
		for( size_t i=0 ; i<subPolygons.size() ; i++ ) for( size_t j=0 ; j<subPolygons[i].size() ; j++ )
		{
			long long oldIdx = polygons[i][j];
			long long newIdx = vMap[ oldIdx ];
			subVertices[ newIdx ] = vertices[ oldIdx ];
			subPolygons[i][j] = newIdx;
		}
		if( Out.set ) WriteMesh( Out.value , ASCII.set ? PLY_ASCII : ft , subVertices , subPolygons , comments );
	}
	else
	{
		if( Out.set ) WriteMesh( Out.value , ASCII.set ? PLY_ASCII : ft , vertices , polygons , comments );
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
	typedef MultiPointStreamData< float , PointStreamValue< float > , PointStreamNormal< float , DIMENSION > , PointStreamColor< float > > VertexData;
	typedef PlyVertexWithData< float , DIMENSION , VertexData > Vertex;
	bool readFlags[ Vertex::PlyReadNum ];
	if( !PlyReadHeader( In.value , Vertex::PlyReadProperties() , Vertex::PlyReadNum , readFlags ) ) ERROR_OUT( "Failed to read ply header: " , In.value );

	bool hasValue  = VertexData::ValidPlyReadProperties< 0 >( readFlags + DIMENSION );
	bool hasNormal = VertexData::ValidPlyReadProperties< 1 >( readFlags + DIMENSION );
	bool hasColor  = VertexData::ValidPlyReadProperties< 2 >( readFlags + DIMENSION );

	if( hasValue )
		if( hasColor )
			if( hasNormal ) return Execute< PointStreamValue< float > , PointStreamNormal< float , DIMENSION > , PointStreamColor< float > >();
			else            return Execute< PointStreamValue< float > ,                                          PointStreamColor< float > >();
		else
			if( hasNormal ) return Execute< PointStreamValue< float > , PointStreamNormal< float , DIMENSION >                             >();
			else            return Execute< PointStreamValue< float >                                                                      >();
	else
		if( hasColor )
			if( hasNormal ) return Execute< PointStreamNormal< float , DIMENSION > , PointStreamColor< float > >();
			else            return Execute<                                          PointStreamColor< float > >();
		else
			if( hasNormal ) return Execute< PointStreamNormal< float , DIMENSION >                             >();
			else            return Execute<                                                                    >();
}
