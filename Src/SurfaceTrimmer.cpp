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

#undef ARRAY_DEBUG
#define DIMENSION 3

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include "FEMTree.h"
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "MAT.h"
#include "Geometry.h"
#include "Ply.h"
#include "PointStreamData.h"

MessageWriter messageWriter;


cmdLineParameter< char* >
	In( "in" ) ,
	Out( "out" );
cmdLineParameter< int >
	Smooth( "smooth" , 5 );
cmdLineParameter< float >
	Trim( "trim" ) ,
	IslandAreaRatio( "aRatio" , 0.001f );
cmdLineReadable
	PolygonMesh( "polygonMesh" ) ,
#ifdef NEW_CODE
	Long( "long" ) ,
#endif // NEW_CODE
	Verbose( "verbose" );


cmdLineReadable* params[] =
{
	&In , &Out , &Trim , &PolygonMesh , &Smooth , &IslandAreaRatio , &Verbose ,
#ifdef NEW_CODE
	&Long ,
#endif // NEW_CODE
	NULL
};

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input polygon mesh>\n" , In.name );
	printf( "\t --%s <trimming value>\n" , Trim.name );
	printf( "\t[--%s <ouput polygon mesh>]\n" , Out.name );
	printf( "\t[--%s <smoothing iterations>=%d]\n" , Smooth.name , Smooth.value );
	printf( "\t[--%s <relative area of islands>=%f]\n" , IslandAreaRatio.name , IslandAreaRatio.value );
	printf( "\t[--%s]\n" , PolygonMesh.name );
#ifdef NEW_CODE
	printf( "\t[--%s]\n" , Long.name );
#endif // NEW_CODE
	printf( "\t[--%s]\n" , Verbose.name );
}

#ifdef NEW_CODE
template< typename Index >
struct EdgeKey
{
	Index key1 , key2;
	EdgeKey( Index k1=0 , Index k2=0 ) : key1(k1) , key2(k2) {}
	bool operator == ( const EdgeKey &key ) const  { return key1==key.key1 && key2==key.key2; }
	struct Hasher{ size_t operator()( const EdgeKey &key ) const { return key.key1 ^ key.key2; } };
};
#else // !NEW_CODE
long long EdgeKey( int key1 , int key2 )
{
	if( key1<key2 ) return ( ( (long long)key1 )<<32 ) | ( (long long)key2 );
	else            return ( ( (long long)key2 )<<32 ) | ( (long long)key1 );
}
#endif // NEW_CODE

template< typename Real , typename ... VertexData >
PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , PointStreamValue< float > , VertexData ... > > InterpolateVertices( const PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , PointStreamValue< float > , VertexData ... > >& v1 , const PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , PointStreamValue< float > , VertexData ... > >& v2 , Real value )
{
	if( std::get<0>( v1.data.data ).data==std::get<0>( v2.data.data ).data ) return (v1+v2)/Real(2.);
	Real dx = ( std::get<0>( v1.data.data ).data-value ) / ( std::get<0>( v1.data.data ).data-std::get<0>( v2.data.data ).data );
	return v1*(1.f-dx) + v2*dx;
}

#ifdef NEW_CODE
template< typename Real , typename Index , typename ... VertexData >
void SmoothValues( std::vector< PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , PointStreamValue< float > , VertexData ... > > >& vertices , const std::vector< std::vector< Index > >& polygons )
#else // !NEW_CODE
template< typename Real , typename ... VertexData >
void SmoothValues( std::vector< PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , PointStreamValue< float > , VertexData ... > > >& vertices , const std::vector< std::vector< int > >& polygons )
#endif // NEW_CODE
{
	std::vector< int > count( vertices.size() );
	std::vector< Real > sums( vertices.size() , 0 );
	for( size_t i=0 ; i<polygons.size() ; i++ )
	{
		int sz = int(polygons[i].size());
		for( int j=0 ; j<sz ; j++ )
		{
			int j1 = j , j2 = (j+1)%sz;
#ifdef NEW_CODE
			Index v1 = polygons[i][j1] , v2 = polygons[i][j2];
#else // !NEW_CODE
			int v1 = polygons[i][j1] , v2 = polygons[i][j2];
#endif // NEW_CODE
			count[v1]++ , count[v2]++;
			sums[v1] += std::get< 0 >( vertices[v2].data.data ).data , sums[v2] += std::get< 0 >( vertices[v1].data.data ).data;
		}
	}
	for( size_t i=0 ; i<vertices.size() ; i++ ) std::get< 0 >( vertices[i].data.data ).data = ( sums[i] + std::get< 0 >( vertices[i].data.data ).data ) / ( count[i] + 1 );
}

#ifdef NEW_CODE
template< class Real , typename Index , typename ... VertexData >
void SplitPolygon
(
	const std::vector< Index >& polygon ,
	std::vector< PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , PointStreamValue< float > , VertexData ... > > >& vertices ,
	std::vector< std::vector< Index > >* ltPolygons , std::vector< std::vector< Index > >* gtPolygons ,
	std::vector< bool >* ltFlags , std::vector< bool >* gtFlags ,
	std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher >& vertexTable,
	Real trimValue
)
#else // !NEW_CODE
template< class Real , typename ... VertexData >
void SplitPolygon
(
	const std::vector< int >& polygon ,
	std::vector< PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , PointStreamValue< float > , VertexData ... > > >& vertices ,
	std::vector< std::vector< int > >* ltPolygons , std::vector< std::vector< int > >* gtPolygons ,
	std::vector< bool >* ltFlags , std::vector< bool >* gtFlags ,
	std::unordered_map< long long, int >& vertexTable,
	Real trimValue
)
#endif // NEW_CODE
{
	int sz = int( polygon.size() );
	std::vector< bool > gt( sz );
	int gtCount = 0;
	for( int j=0 ; j<sz ; j++ )
	{
		gt[j] = ( std::get<0>( vertices[ polygon[j] ].data.data ).data>trimValue );
		if( gt[j] ) gtCount++;
	}
	if     ( gtCount==sz ){ if( gtPolygons ) gtPolygons->push_back( polygon ) ; if( gtFlags ) gtFlags->push_back( false ); }
	else if( gtCount==0  ){ if( ltPolygons ) ltPolygons->push_back( polygon ) ; if( ltFlags ) ltFlags->push_back( false ); }
	else
	{
		int start;
		for( start=0 ; start<sz ; start++ ) if( gt[start] && !gt[(start+sz-1)%sz] ) break;

		bool gtFlag = true;
#ifdef NEW_CODE
		std::vector< Index > poly;
#else // !NEW_CODE
		std::vector< int > poly;
#endif // NEW_CODE

		// Add the initial vertex
		{
			int j1 = (start+int(sz)-1)%sz , j2 = start;
#ifdef NEW_CODE
			Index v1 = polygon[j1] , v2 = polygon[j2] , vIdx;
			typename std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher >::iterator iter = vertexTable.find( EdgeKey< Index >(v1,v2) );
#else // !NEW_CODE
			int v1 = polygon[j1] , v2 = polygon[j2];
			int vIdx;
			std::unordered_map< long long, int >::iterator iter = vertexTable.find( EdgeKey(v1,v2) );
#endif // NEW_CODE
			if( iter==vertexTable.end() )
			{
#ifdef NEW_CODE
				vertexTable[ EdgeKey< Index >(v1,v2) ] = vIdx = (Index)vertices.size();
#else // !NEW_CODE
				vertexTable[ EdgeKey( v1 , v2 ) ] = vIdx = int( vertices.size() );
#endif // NEW_CODE
				vertices.push_back( InterpolateVertices( vertices[v1] , vertices[v2] , trimValue ) );
			}
			else vIdx = iter->second;
			poly.push_back( vIdx );
		}

		for( int _j=0  ; _j<=sz ; _j++ )
		{
			int j1 = (_j+start+sz-1)%sz , j2 = (_j+start)%sz;
#ifdef NEW_CODE
			Index v1 = polygon[j1] , v2 = polygon[j2];
#else // !NEW_CODE
			int v1 = polygon[j1] , v2 = polygon[j2];
#endif // NEW_CODE
			if( gt[j2]==gtFlag ) poly.push_back( v2 );
			else
			{
#ifdef NEW_CODE
				Index vIdx;
				typename std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher >::iterator iter = vertexTable.find( EdgeKey< Index >(v1,v2) );
				if( iter==vertexTable.end() )
				{
					vertexTable[ EdgeKey< Index >(v1,v2) ] = vIdx = (Index)vertices.size();
					vertices.push_back( InterpolateVertices( vertices[v1] , vertices[v2] , trimValue ) );
				}
#else // !NEW_CODE
				int vIdx;
				std::unordered_map< long long, int >::iterator iter = vertexTable.find(EdgeKey(v1, v2));
				if( iter==vertexTable.end() )
				{
					vertexTable[ EdgeKey( v1 , v2 ) ] = vIdx = int( vertices.size() );
					vertices.push_back( InterpolateVertices( vertices[v1] , vertices[v2] , trimValue ) );
				}
#endif // NEW_CODE
				else vIdx = iter->second;
				poly.push_back( vIdx );
				if( gtFlag ){ if( gtPolygons ) gtPolygons->push_back( poly ) ; if( ltFlags ) ltFlags->push_back( true ); }
				else        { if( ltPolygons ) ltPolygons->push_back( poly ) ; if( gtFlags ) gtFlags->push_back( true ); }
				poly.clear() , poly.push_back( vIdx ) , poly.push_back( v2 );
				gtFlag = !gtFlag;
			}
		}
	}
}

#ifdef NEW_CODE
template< class Real , typename Index , class Vertex >
void Triangulate( const std::vector< Vertex >& vertices , const std::vector< std::vector< Index > >& polygons , std::vector< std::vector< Index > >& triangles )
#else // !NEW_CODE
template< class Real , class Vertex >
void Triangulate( const std::vector< Vertex >& vertices , const std::vector< std::vector< int > >& polygons , std::vector< std::vector< int > >& triangles )
#endif // NEW_CODE
{
	triangles.clear();
	for( size_t i=0 ; i<polygons.size() ; i++ )
		if( polygons.size()>3 )
		{
			std::vector< Point< Real , DIMENSION > > _vertices( polygons[i].size() );
			for( int j=0 ; j<int( polygons[i].size() ) ; j++ ) _vertices[j] = vertices[ polygons[i][j] ].point;
#ifdef NEW_CODE
			std::vector< TriangleIndex< Index > > _triangles = MinimalAreaTriangulation< Index , Real , DIMENSION >( ( ConstPointer( Point< Real , DIMENSION > ) )GetPointer( _vertices ) , _vertices.size() );
#else // !NEW_CODE
			std::vector< TriangleIndex > _triangles = MinimalAreaTriangulation< Real , DIMENSION >( ( ConstPointer( Point< Real , DIMENSION > ) )GetPointer( _vertices ) , _vertices.size() );
#endif // NEW_CODE

			// Add the triangles to the mesh
			size_t idx = triangles.size();
			triangles.resize( idx+_triangles.size() );
			for( int j=0 ; j<int(_triangles.size()) ; j++ )
			{
				triangles[idx+j].resize(3);
				for( int k=0 ; k<3 ; k++ ) triangles[idx+j][k] = polygons[i][ _triangles[j].idx[k] ];
			}
		}
		else if( polygons[i].size()==3 ) triangles.push_back( polygons[i] );
}

#ifdef NEW_CODE
template< class Real , typename Index , class Vertex >
double PolygonArea( const std::vector< Vertex >& vertices , const std::vector< Index >& polygon )
#else // !NEW_CODE
template< class Real , class Vertex >
double PolygonArea( const std::vector< Vertex >& vertices , const std::vector< int >& polygon )
#endif // NEW_CODE
{
	if( polygon.size()<3 ) return 0.;
	else if( polygon.size()==3 ) return Area( vertices[polygon[0]].point , vertices[polygon[1]].point , vertices[polygon[2]].point );
	else
	{
		Point< Real , DIMENSION > center;
		for( size_t i=0 ; i<polygon.size() ; i++ ) center += vertices[ polygon[i] ].point;
		center /= Real( polygon.size() );
		double area = 0;
		for( size_t i=0 ; i<polygon.size() ; i++ ) area += Area( center , vertices[ polygon[i] ].point , vertices[ polygon[ (i+1)%polygon.size() ] ].point );
		return area;
	}
}

#ifdef NEW_CODE
template< typename Index , class Vertex >
void RemoveHangingVertices( std::vector< Vertex >& vertices , std::vector< std::vector< Index > >& polygons )
#else // !NEW_CODE
template< class Vertex >
void RemoveHangingVertices( std::vector< Vertex >& vertices , std::vector< std::vector< int > >& polygons )
#endif // NEW_CODE
{
#ifdef NEW_CODE
	std::unordered_map< Index, Index > vMap;
#else // !NEW_CODE
	std::unordered_map< int, int > vMap;
#endif // NEW_CODE
	std::vector< bool > vertexFlags( vertices.size() , false );
	for( size_t i=0 ; i<polygons.size() ; i++ ) for( size_t j=0 ; j<polygons[i].size() ; j++ ) vertexFlags[ polygons[i][j] ] = true;
#ifdef NEW_CODE
	Index vCount = 0;
	for( Index i=0 ; i<(Index)vertices.size() ; i++ ) if( vertexFlags[i] ) vMap[i] = vCount++;
#else // !NEW_CODE
	int vCount = 0;
	for( int i=0 ; i<int(vertices.size()) ; i++ ) if( vertexFlags[i] ) vMap[i] = vCount++;
#endif // NEW_CODE
	for( size_t i=0 ; i<polygons.size() ; i++ ) for( size_t j=0 ; j<polygons[i].size() ; j++ ) polygons[i][j] = vMap[ polygons[i][j] ];

	std::vector< Vertex > _vertices( vCount );
#ifdef NEW_CODE
	for( Index i=0 ; i<(Index)vertices.size() ; i++ ) if( vertexFlags[i] ) _vertices[ vMap[i] ] = vertices[i];
#else // !NEW_CODE
	for( int i=0 ; i<int(vertices.size()) ; i++ ) if( vertexFlags[i] ) _vertices[ vMap[i] ] = vertices[i];
#endif // NEW_CODE
	vertices = _vertices;
}

#ifdef NEW_CODE
template< typename Index >
void SetConnectedComponents( const std::vector< std::vector< Index > >& polygons , std::vector< std::vector< Index > >& components )
#else // !NEW_CODE
void SetConnectedComponents( const std::vector< std::vector< int > >& polygons , std::vector< std::vector< int > >& components )
#endif // NEW_CODE
{
#ifdef NEW_CODE
	std::vector< Index > polygonRoots( polygons.size() );
	for( size_t i=0 ; i<polygons.size() ; i++ ) polygonRoots[i] = (Index)i;
	std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher > edgeTable;
#else // !NEW_CODE
	std::vector< int > polygonRoots( polygons.size() );
	for( size_t i=0 ; i<polygons.size() ; i++ ) polygonRoots[i] = int(i);
	std::unordered_map< long long, int > edgeTable;
#endif // NEW_CODE
	for( size_t i=0 ; i<polygons.size() ; i++ )
	{
		int sz = int( polygons[i].size() );
		for( int j=0 ; j<sz ; j++ )
		{
			int j1 = j , j2 = (j+1)%sz;
#ifdef NEW_CODE
			Index v1 = polygons[i][j1] , v2 = polygons[i][j2];
			EdgeKey< Index > eKey = EdgeKey< Index >(v1,v2);
			typename std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher >::iterator iter = edgeTable.find(eKey);
			if( iter==edgeTable.end() ) edgeTable[ eKey ] = (Index)i;
#else // !NEW_CODE
			int v1 = polygons[i][j1] , v2 = polygons[i][j2];
			long long eKey = EdgeKey( v1 , v2 );
			std::unordered_map< long long, int >::iterator iter = edgeTable.find(eKey);
			if( iter==edgeTable.end() ) edgeTable[ eKey ] = int(i);
#endif // NEW_CODE
			else
			{
#ifdef NEW_CODE
				Index p = iter->second;
#else // !NEW_CODE
				int p = iter->second;
#endif // NEW_CODE
				while( polygonRoots[p]!=p )
				{
#ifdef NEW_CODE
					Index temp = polygonRoots[p];
					polygonRoots[p] = (Index)i;
#else // !NEW_CODE
					int temp = polygonRoots[p];
					polygonRoots[p] = int(i);
#endif // NEW_CODE
					p = temp;
				}
#ifdef NEW_CODE
				polygonRoots[p] = (Index)i;
#else // !NEW_CODE
				polygonRoots[p] = int(i);
#endif // NEW_CODE
			}
		}
	}
	for( size_t i=0 ; i<polygonRoots.size() ; i++ )
	{
#ifdef NEW_CODE
		Index p = (Index)i;
#else // !NEW_CODE
		int p = int(i);
#endif // NEW_CODE
		while( polygonRoots[p]!=p ) p = polygonRoots[p];
#ifdef NEW_CODE
		Index root = p;
		p = (Index)i;
#else // !NEW_CODE
		int root = p;
		p = int(i);
#endif // NEW_CODE
		while( polygonRoots[p]!=p )
		{
#ifdef NEW_CODE
			Index temp = polygonRoots[p];
#else // !NEW_CODE
			int temp = polygonRoots[p];
#endif // NEW_CODE
			polygonRoots[p] = root;
			p = temp;
		}
	}
	int cCount = 0;
#ifdef NEW_CODE
	std::unordered_map< Index , Index > vMap;
	for( Index i=0 ; i<(Index)polygonRoots.size() ; i++ ) if( polygonRoots[i]==i ) vMap[i] = cCount++;
#else // !NEW_CODE
	std::unordered_map< int , int > vMap;
	for( int i=0 ; i<int(polygonRoots.size()) ; i++ ) if( polygonRoots[i]==i ) vMap[i] = cCount++;
#endif //  NEW_CODE
	components.resize( cCount );
#ifdef NEW_CODE
	for( Index i=0 ; i<(Index)polygonRoots.size() ; i++ ) components[ vMap[ polygonRoots[i] ] ].push_back(i);
#else // !NEW_CODE
	for( int i=0 ; i<int(polygonRoots.size()) ; i++ ) components[ vMap[ polygonRoots[i] ] ].push_back(i);
#endif // NEW_CODE
}

#ifdef NEW_CODE
template< typename Index , typename ... VertexData >
#else // !NEW_CODE
template< typename ... VertexData >
#endif // NEW_CODE
int Execute( void )
{
	typedef PlyVertexWithData< float , DIMENSION , MultiPointStreamData< float , PointStreamValue< float > , VertexData ... > > Vertex;
	float min , max;
	std::vector< Vertex > vertices;
#ifdef NEW_CODE
	std::vector< std::vector< Index > > polygons;
#else // !NEW_CODE
	std::vector< std::vector< int > > polygons;
#endif // NEW_CODE

	int ft;
	std::vector< std::string > comments;
	PlyReadPolygons< Vertex >( In.value , vertices , polygons , Vertex::PlyReadProperties() , Vertex::PlyReadNum , ft , comments );

#ifdef NEW_CODE
	for( int i=0 ; i<Smooth.value ; i++ ) SmoothValues< float , Index >( vertices , polygons );
#else // !NEW_CODE
	for( int i=0 ; i<Smooth.value ; i++ ) SmoothValues< float >( vertices , polygons );
#endif // NEW_CODE
	min = max = std::get< 0 >( vertices[0].data.data ).data;
	for( size_t i=0 ; i<vertices.size() ; i++ ) min = std::min< float >( min , std::get< 0 >( vertices[i].data.data ).data ) , max = std::max< float >( max , std::get< 0 >( vertices[i].data.data ).data );

#ifdef NEW_CODE
	std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher > vertexTable;
	std::vector< std::vector< Index > > ltPolygons , gtPolygons;
#else // !NEW_CODE
	std::unordered_map< long long, int > vertexTable;
	std::vector< std::vector< int > > ltPolygons , gtPolygons;
#endif // NEW_CODE
	std::vector< bool > ltFlags , gtFlags;

	messageWriter( comments , "*********************************************\n" );
	messageWriter( comments , "*********************************************\n" );
	messageWriter( comments , "** Running Surface Trimmer (Version %s) **\n" , VERSION );
	messageWriter( comments , "*********************************************\n" );
	messageWriter( comments , "*********************************************\n" );
	char str[1024];
	for( int i=0 ; params[i] ; i++ )
		if( params[i]->set )
		{
			params[i]->writeValue( str );
			if( strlen( str ) ) messageWriter( comments , "\t--%s %s\n" , params[i]->name , str );
			else                messageWriter( comments , "\t--%s\n" , params[i]->name );
		}
	if( Verbose.set ) printf( "Value Range: [%f,%f]\n" , min , max );

	double t=Time();
	for( size_t i=0 ; i<polygons.size() ; i++ ) SplitPolygon( polygons[i] , vertices , &ltPolygons , &gtPolygons , &ltFlags , &gtFlags , vertexTable , Trim.value );
	if( IslandAreaRatio.value>0 )
	{
#ifdef NEW_CODE
		std::vector< std::vector< Index > > _ltPolygons , _gtPolygons;
		std::vector< std::vector< Index > > ltComponents , gtComponents;
#else // !NEW_CODE
		std::vector< std::vector< int > > _ltPolygons , _gtPolygons;
		std::vector< std::vector< int > > ltComponents , gtComponents;
#endif // NEW_CODE
		SetConnectedComponents( ltPolygons , ltComponents );
		SetConnectedComponents( gtPolygons , gtComponents );
		std::vector< double > ltAreas( ltComponents.size() , 0. ) , gtAreas( gtComponents.size() , 0. );
		std::vector< bool > ltComponentFlags( ltComponents.size() , false ) , gtComponentFlags( gtComponents.size() , false );
		double area = 0.;
		for( size_t i=0 ; i<ltComponents.size() ; i++ )
		{
			for( size_t j=0 ; j<ltComponents[i].size() ; j++ )
			{
#ifdef NEW_CODE
				ltAreas[i] += PolygonArea< float , Index , Vertex >( vertices , ltPolygons[ ltComponents[i][j] ] );
#else // !NEW_CODE
				ltAreas[i] += PolygonArea< float , Vertex >( vertices , ltPolygons[ ltComponents[i][j] ] );
#endif // NEW_CODE
				ltComponentFlags[i] = ( ltComponentFlags[i] || ltFlags[ ltComponents[i][j] ] );
			}
			area += ltAreas[i];
		}
		for( size_t i=0 ; i<gtComponents.size() ; i++ )
		{
			for( size_t j=0 ; j<gtComponents[i].size() ; j++ )
			{
#ifdef NEW_CODE
				gtAreas[i] += PolygonArea< float , Index , Vertex >( vertices , gtPolygons[ gtComponents[i][j] ] );
#else // !NEW_CODE
				gtAreas[i] += PolygonArea< float , Vertex >( vertices , gtPolygons[ gtComponents[i][j] ] );
#endif // NEW_CODE
				gtComponentFlags[i] = ( gtComponentFlags[i] || gtFlags[ gtComponents[i][j] ] );
			}
			area += gtAreas[i];
		}
		for( size_t i=0 ; i<ltComponents.size() ; i++ )
		{
			if( ltAreas[i]<area*IslandAreaRatio.value && ltComponentFlags[i] ) for( size_t j=0 ; j<ltComponents[i].size() ; j++ ) _gtPolygons.push_back( ltPolygons[ ltComponents[i][j] ] );
			else                                                               for( size_t j=0 ; j<ltComponents[i].size() ; j++ ) _ltPolygons.push_back( ltPolygons[ ltComponents[i][j] ] );
		}
		for( size_t i=0 ; i<gtComponents.size() ; i++ )
		{
			if( gtAreas[i]<area*IslandAreaRatio.value && gtComponentFlags[i] ) for( size_t j=0 ; j<gtComponents[i].size() ; j++ ) _ltPolygons.push_back( gtPolygons[ gtComponents[i][j] ] );
			else                                                               for( size_t j=0 ; j<gtComponents[i].size() ; j++ ) _gtPolygons.push_back( gtPolygons[ gtComponents[i][j] ] );
		}
		ltPolygons = _ltPolygons , gtPolygons = _gtPolygons;
	}
	if( !PolygonMesh.set )
	{
		{
#ifdef NEW_CODE
			std::vector< std::vector< Index > > polys = ltPolygons;
			Triangulate< float , Index , Vertex >( vertices , ltPolygons , polys ) , ltPolygons = polys;
#else // !NEW_CODE
			std::vector< std::vector< int > > polys = ltPolygons;
			Triangulate< float , Vertex >( vertices , ltPolygons , polys ) , ltPolygons = polys;
#endif // NEW_CODE
		}
		{
#ifdef NEW_CODE
			std::vector< std::vector< Index > > polys = gtPolygons;
			Triangulate< float  , Index , Vertex >( vertices , gtPolygons , polys ) , gtPolygons = polys;
#else // !NEW_CODE
			std::vector< std::vector< int > > polys = gtPolygons;
			Triangulate< float  , Vertex >( vertices , gtPolygons , polys ) , gtPolygons = polys;
#endif // NEW_CODE
		}
	}

	RemoveHangingVertices( vertices , gtPolygons );
	char comment[1024];
	sprintf( comment , "#Trimmed In: %9.1f (s)" , Time()-t );
	comments.push_back( comment );
	if( Out.set )
		if( !PlyWritePolygons< Vertex >( Out.value , vertices , gtPolygons , Vertex::PlyWriteProperties() , Vertex::PlyWriteNum , ft , comments ) )
			ERROR_OUT( "Could not write mesh to: " , Out.value );
	
	return EXIT_SUCCESS;
}
int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , &argv[1] , params );
	messageWriter.echoSTDOUT = Verbose.set;

	if( !In.set || !Trim.set )
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

	if( !hasValue ) ERROR_OUT( "Ply file does not contain values" );

#ifdef NEW_CODE
	if( Long.set )
		if( hasColor )
			if( hasNormal ) return Execute< long long , PointStreamNormal< float , DIMENSION > , PointStreamColor< float > >();
			else            return Execute< long long ,                                          PointStreamColor< float > >();
		else
			if( hasNormal ) return Execute< long long , PointStreamNormal< float , DIMENSION >                             >();
			else            return Execute< long long                                                                      >();
	else
		if( hasColor )
			if( hasNormal ) return Execute< int , PointStreamNormal< float , DIMENSION > , PointStreamColor< float > >();
			else            return Execute< int ,                                          PointStreamColor< float > >();
		else
			if( hasNormal ) return Execute< int , PointStreamNormal< float , DIMENSION >                             >();
			else            return Execute< int                                                                      >();
#else // !NEW_CODE
	if( hasColor )
		if( hasNormal ) return Execute< PointStreamNormal< float , DIMENSION > , PointStreamColor< float > >();
		else            return Execute<                                          PointStreamColor< float > >();
	else
		if( hasNormal ) return Execute< PointStreamNormal< float , DIMENSION >                             >();
		else            return Execute<                                                                    >();
#endif // NEW_CODE
}
