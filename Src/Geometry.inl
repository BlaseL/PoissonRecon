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

#include <stdio.h>
#include "MyMiscellany.h"

template< class Real > Real Random( void ){ return Real( rand() )/RAND_MAX; }

template< class Real , int Dim >
Point< Real , Dim > RandomBallPoint( void )
{
	Point< Real , Dim > p;
	while(1)
	{
		for( int d=0 ; d<Dim ; d++ ) p[d] = Real( 1.0-2.0*Random< Real >() );
		double l=SquareLength(p);
		if( SquareLength( p )<=1 ) return p;
	}
}
template< class Real , int Dim >
Point< Real , Dim > RandomSpherePoint( void )
{
	Point< Real , Dim > p = RandomBallPoint< Real , Dim >();
	return p / (Real)Length( p );
}

#ifdef NEW_CODE
#else // !NEW_CODE
///////////////////
// Triangulation //
///////////////////
template< class Real , unsigned int Dim >
long long Triangulation< Real , Dim >::EdgeIndex( unsigned int p1 , unsigned int p2 )
{
	if( p1>p2 ) return ((long long)(p1)<<32) | ((long long)(p2));
	else        return ((long long)(p2)<<32) | ((long long)(p1));
}

template< class Real , unsigned int Dim >
unsigned int Triangulation< Real , Dim >::factor( unsigned int tIndex , unsigned int &p1 , unsigned int &p2 , unsigned int &p3 ) const
{
	if( triangles[tIndex].eIndex[0]<0 || triangles[tIndex].eIndex[1]<0 || triangles[tIndex].eIndex[2]<0 ) return 0;
	if( edges[triangles[tIndex].eIndex[0]].tIndex[0]==tIndex ) p1=edges[triangles[tIndex].eIndex[0]].pIndex[0];
	else                                                       p1=edges[triangles[tIndex].eIndex[0]].pIndex[1];
	if( edges[triangles[tIndex].eIndex[1]].tIndex[0]==tIndex ) p2=edges[triangles[tIndex].eIndex[1]].pIndex[0];
	else                                                       p2=edges[triangles[tIndex].eIndex[1]].pIndex[1];
	if( edges[triangles[tIndex].eIndex[2]].tIndex[0]==tIndex ) p3=edges[triangles[tIndex].eIndex[2]].pIndex[0];
	else                                                       p3=edges[triangles[tIndex].eIndex[2]].pIndex[1];
	return 1;
}
template< class Real , unsigned int Dim > Real Triangulation< Real , Dim >::area( unsigned int p1 , unsigned int p2 , unsigned int p3 ) const { return Area( points[p1] , points[p2] , points[p3] ); }
template< class Real , unsigned int Dim >
Real Triangulation< Real , Dim >::area( unsigned int tIndex ) const
{
	unsigned int p1 , p2 , p3;
	factor( tIndex , p1 , p2 , p3 );
	return area(p1,p2,p3);
}
template< class Real , unsigned int Dim >
Real Triangulation< Real , Dim >::area( void ) const
{
	Real a=0;
	for( int i=0 ; i<(int)triangles.size() ; i++ ) a += area(i);
	return a;
}
template< class Real , unsigned int Dim >
unsigned int Triangulation< Real , Dim >::addTriangle( unsigned int p1 , unsigned int p2 , unsigned int p3 )
{
	std::unordered_map<long long, int>::iterator iter;
	unsigned int tIdx , eIdx , p[] = { p1 , p2 , p3 };
	triangles.push_back( TriangulationTriangle() );
	tIdx = (int)triangles.size()-1;

	for( int i=0 ; i<3 ; i++ )
	{
		long long e = EdgeIndex( p[i] , p[(i+1)%3] );
		iter = edgeMap.find(e);
		if( iter==edgeMap.end() )
		{
			TriangulationEdge edge;
			edge.pIndex[0] = p[i];
			edge.pIndex[1] = p[(i+1)%3];
			edges.push_back(edge);
			eIdx = (int)edges.size()-1;
			edgeMap[e] = eIdx;
			edges[eIdx].tIndex[0]=tIdx;
		}
		else
		{
			eIdx = edgeMap[e];
			if( edges[eIdx].pIndex[0]==p[i] )
			{
				if( edges[eIdx].tIndex[0]<0 ) edges[eIdx].tIndex[0] = tIdx;
				else{ printf( "Edge Triangle in use 1\n" ) ; return 0; }
			}
			else
			{
				if( edges[eIdx].tIndex[1]<0 ) edges[eIdx].tIndex[1] = tIdx;
				else{ printf( "Edge Triangle in use 2\n") ; return 0; }
			}

		}
		triangles[tIdx].eIndex[i] = eIdx;
	}
	return tIdx;
}
template< class Real , unsigned int Dim >
unsigned int Triangulation< Real , Dim >::flipMinimize( unsigned int eIndex )
{
	Real oldArea,newArea;
	int oldP[3] , oldQ[3] , newP[3] , newQ[3];
	TriangulationEdge newEdge;

	if( edges[eIndex].tIndex[0]<0 || edges[eIndex].tIndex[1]<0 ) return 0;

	if( !factor( edges[eIndex].tIndex[0] , oldP[0] , oldP[1] , oldP[2] ) ) return 0;
	if( !factor( edges[eIndex].tIndex[1] , oldQ[0] , oldQ[1] , oldQ[2] ) ) return 0;

	oldArea = area( oldP[0] , oldP[1] , oldP[2] ) + area( oldQ[0] , oldQ[1] , oldQ[2] );
	unsigned int idxP , idxQ;
	for( idxP=0 ; idxP<3 ; idxP++ )
	{
		unsigned int i;
		for( i=0 ; i<3 ; i++ ) if( oldP[idxP]==oldQ[i] ) break;
		if(i==3) break;
	}
	for( idxQ=0 ; idxQ<3 ; idxQ++ )
	{
		unsigned int i;
		for( i=0 ; i<3 ; i++ ) if( oldP[i]==oldQ[idxQ] ) break;
		if( i==3 ) break;
	}
	if(idxP==3 || idxQ==3) return 0;
	newP[0]=oldP[idxP];
	newP[1]=oldP[(idxP+1)%3];
	newP[2]=oldQ[idxQ];
	newQ[0]=oldQ[idxQ];
	newQ[1]=oldP[(idxP+2)%3];
	newQ[2]=oldP[idxP];

	newArea = area( newP[0] , newP[1] , newP[2] ) + area( newQ[0] , newQ[1] , newQ[2] );
	if( oldArea<=newArea ) return 0;

	// Remove the entry in the hash-table for the old edge
	edgeMap.erase( EdgeIndex( edges[eIndex].pIndex[0] , edges[eIndex].pIndex[1] ) );
	// Set the new edge so that the zero-side is newQ
	edges[eIndex].pIndex[0] = newP[0];
	edges[eIndex].pIndex[1] = newQ[0];
	// Insert the entry into the hash-table for the new edge
	edgeMap[EdgeIndex(newP[0],newQ[0])] = eIndex;
	// Update the triangle information
	for( unsigned int i=0 ; i<3 ; i++ )
	{
		int idx;
		idx = edgeMap[EdgeIndex(newQ[i],newQ[(i+1)%3])];
		triangles[edges[eIndex].tIndex[0]].eIndex[i] = idx;
		if(idx!=eIndex)
		{
			if( edges[idx].tIndex[0]==edges[eIndex].tIndex[1] ) edges[idx].tIndex[0] = edges[eIndex].tIndex[0];
			if( edges[idx].tIndex[1]==edges[eIndex].tIndex[1] ) edges[idx].tIndex[1] = edges[eIndex].tIndex[0];
		}

		idx = edgeMap[EdgeIndex(newP[i],newP[(i+1)%3])];
		triangles[edges[eIndex].tIndex[1]].eIndex[i]=idx;
		if( idx!=eIndex )
		{
			if( edges[idx].tIndex[0]==edges[eIndex].tIndex[0] ) edges[idx].tIndex[0]=edges[eIndex].tIndex[1];
			if( edges[idx].tIndex[1]==edges[eIndex].tIndex[0] ) edges[idx].tIndex[1]=edges[eIndex].tIndex[1];
		}
	}
	return 1;
}
#endif // NEW_CODE

/////////////////////////
// CoredVectorMeshData //
/////////////////////////
#ifdef NEW_CODE
template< class Vertex , typename Index >
#ifdef NEW_THREADS
CoredVectorMeshData< Vertex , Index >::CoredVectorMeshData( void ) { oocPointIndex = polygonIndex = threadIndex = 0 ; polygons.resize( std::thread::hardware_concurrency() ); }
#else // !NEW_THREADS
CoredVectorMeshData< Vertex , Index >::CoredVectorMeshData( void ) { oocPointIndex = polygonIndex = threadIndex = 0 ; polygons.resize( omp_get_max_threads() ); }
#endif // NEW_THREADS
template< class Vertex , typename Index >
void CoredVectorMeshData< Vertex , Index >::resetIterator ( void ) { oocPointIndex = polygonIndex = threadIndex = 0; }
template< class Vertex , typename Index >
Index CoredVectorMeshData< Vertex , Index >::addOutOfCorePoint( const Vertex& p )
{
	oocPoints.push_back(p);
	return ( Index )oocPoints.size()-1;
}
template< class Vertex , typename Index >
#ifdef NEW_THREADS
Index CoredVectorMeshData< Vertex , Index >::addOutOfCorePoint_s( unsigned int thread , const Vertex& p )
#else // !NEW_THREADS
Index CoredVectorMeshData< Vertex , Index >::addOutOfCorePoint_s( const Vertex& p )
#endif // NEW_THREADS
{
	size_t sz;
#ifdef NEW_THREADS
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock(m);
#else // !NEW_THREADS
#pragma omp critical (CoredVectorMeshData_addOutOfCorePoint_s )
	{
#endif // NEW_THREADS
		sz = oocPoints.size();
		oocPoints.push_back(p);
	}
	return (Index)sz;
}
#ifdef NEW_THREADS
template< class Vertex , typename Index >
void CoredVectorMeshData< Vertex , Index >::addPolygon_s( unsigned int thread , const std::vector< Index >& polygon )
{
	polygons[ thread ].push_back( polygon );
}
#else // !NEW_THREADS
template< class Vertex , typename Index >
void CoredVectorMeshData< Vertex , Index >::addPolygon_s( const std::vector< Index >& polygon )
{
	polygons[ omp_get_thread_num() ].push_back( polygon );
}
#endif // NEW_THREADS
template< class Vertex , typename Index >
#ifdef NEW_THREADS
void CoredVectorMeshData< Vertex , Index >::addPolygon_s( unsigned int thread , const std::vector< CoredVertexIndex< Index > >& vertices )
#else // !NEW_THREADS
void CoredVectorMeshData< Vertex , Index >::addPolygon_s( const std::vector< CoredVertexIndex< Index > >& vertices )
#endif // NEW_THREADS
{
	std::vector< Index > polygon( vertices.size() );
	for( int i=0 ; i<(int)vertices.size() ; i++ ) 
		if( vertices[i].inCore ) polygon[i] =  vertices[i].idx;
		else                     polygon[i] = -vertices[i].idx-1;
#ifdef NEW_THREADS
	return addPolygon_s( thread , polygon );
#else // !NEW_THREADS
	return addPolygon_s( polygon );
#endif // NEW_THREADS
}
template< class Vertex , typename Index >
Index CoredVectorMeshData< Vertex , Index >::nextOutOfCorePoint( Vertex& p )
{
	if( oocPointIndex<(Index)oocPoints.size() )
	{
		p=oocPoints[oocPointIndex++];
		return 1;
	}
	else return 0;
}
template< class Vertex , typename Index >
Index CoredVectorMeshData< Vertex , Index >::nextPolygon( std::vector< CoredVertexIndex< Index > >& vertices )
{
	while( true )
	{
		if( threadIndex<(int)polygons.size() )
		{
			if( polygonIndex<(Index)( polygons[threadIndex].size() ) )
			{
				std::vector< Index >& polygon = polygons[threadIndex][ polygonIndex++ ];
				vertices.resize( polygon.size() );
				for( int i=0 ; i<int(polygon.size()) ; i++ )
					if( polygon[i]<0 ) vertices[i].idx = -polygon[i]-1 , vertices[i].inCore = false;
					else               vertices[i].idx =  polygon[i]   , vertices[i].inCore = true;
				return 1;
			}
			else threadIndex++ , polygonIndex = 0;
		}
		else return 0;
	}
}
template< class Vertex , typename Index >
size_t CoredVectorMeshData< Vertex , Index >::outOfCorePointCount( void ){ return oocPoints.size(); }
template< class Vertex , typename Index >
size_t CoredVectorMeshData< Vertex , Index >::polygonCount( void )
{
	size_t count = 0;
	for( size_t i=0 ; i<polygons.size() ; i++ ) count += polygons[i].size();
	return count;
}
#else // !NEW_CODE
template< class Vertex >
CoredVectorMeshData< Vertex >::CoredVectorMeshData( void ) { oocPointIndex = polygonIndex = threadIndex = 0 ; polygons.resize( omp_get_max_threads() ); }
template< class Vertex >
void CoredVectorMeshData< Vertex >::resetIterator ( void ) { oocPointIndex = polygonIndex = threadIndex = 0; }
template< class Vertex >
int CoredVectorMeshData< Vertex >::addOutOfCorePoint( const Vertex& p )
{
	oocPoints.push_back(p);
	return int(oocPoints.size())-1;
}
template< class Vertex >
#ifdef NEW_THREADS
int CoredVectorMeshData< Vertex >::addOutOfCorePoint_s( unsigned int thread , const Vertex& p )
#else // !NEW_THREADS
int CoredVectorMeshData< Vertex >::addOutOfCorePoint_s( const Vertex& p )
#endif // NEW_THREADS
{
	size_t sz;
#ifdef NEW_THREADS
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock(m);
#else // !NEW_THREADS
#pragma omp critical (CoredVectorMeshData_addOutOfCorePoint_s )
	{
#endif // NEW_THREADS
		sz = oocPoints.size();
		oocPoints.push_back(p);
	}
	return (int)sz;
}
#ifdef NEW_THREADS
template< class Vertex >
void CoredVectorMeshData< Vertex >::addPolygon_s( unsigned int thread , const std::vector< int >& polygon )
{
	polygons[ thread ].push_back( polygon );
}
#else // !NEW_THREADS
template< class Vertex >
void CoredVectorMeshData< Vertex >::addPolygon_s( const std::vector< int >& polygon )
{
	polygons[ omp_get_thread_num() ].push_back( polygon );
}
#endif // NEW_THREADS
#ifdef NEW_THREADS
template< class Vertex >
void CoredVectorMeshData< Vertex >::addPolygon_s( unsigned int thread , const std::vector< CoredVertexIndex >& vertices )
#else // !NEW_THREADS
template< class Vertex >
void CoredVectorMeshData< Vertex >::addPolygon_s( const std::vector< CoredVertexIndex >& vertices )
#endif // NEW_THREADS
{
	std::vector< int > polygon( vertices.size() );
	for( int i=0 ; i<(int)vertices.size() ; i++ ) 
		if( vertices[i].inCore ) polygon[i] =  vertices[i].idx;
		else                     polygon[i] = -vertices[i].idx-1;
#ifdef NEW_THREADS
	return addPolygon_s( thread , polygon );
#else // !NEW_THREADS
	return addPolygon_s( polygon );
#endif // NEW_THREADS
}
template< class Vertex >
int CoredVectorMeshData< Vertex >::nextOutOfCorePoint( Vertex& p )
{
	if( oocPointIndex<int(oocPoints.size()) )
	{
		p=oocPoints[oocPointIndex++];
		return 1;
	}
	else return 0;
}
template< class Vertex >
int CoredVectorMeshData< Vertex >::nextPolygon( std::vector< CoredVertexIndex >& vertices )
{
	while( true )
	{
		if( threadIndex<(int)polygons.size() )
		{
			if( polygonIndex<int( polygons[threadIndex].size() ) )
			{
				std::vector< int >& polygon = polygons[threadIndex][ polygonIndex++ ];
				vertices.resize( polygon.size() );
				for( int i=0 ; i<int(polygon.size()) ; i++ )
					if( polygon[i]<0 ) vertices[i].idx = -polygon[i]-1 , vertices[i].inCore = false;
					else               vertices[i].idx =  polygon[i]   , vertices[i].inCore = true;
					return 1;
			}
			else threadIndex++ , polygonIndex = 0;
		}
		else return 0;
	}
}
template< class Vertex >
int CoredVectorMeshData< Vertex >::outOfCorePointCount(void){return int(oocPoints.size());}
template< class Vertex >
int CoredVectorMeshData< Vertex >::polygonCount( void )
{
	int count = 0;
	for( int i=0 ; i<polygons.size() ; i++ ) count += (int)polygons[i].size();
	return count;
}
#endif // NEW_CODE

///////////////////////
// CoredFileMeshData //
///////////////////////
#ifdef NEW_CODE
template< class Vertex , typename Index >
CoredFileMeshData< Vertex , Index >::CoredFileMeshData( const char* fileHeader )
{
	threadIndex = 0;
	oocPoints = 0;
#ifdef NEW_THREADS
	polygons.resize( std::thread::hardware_concurrency() );
#else // !NEW_THREADS
	polygons.resize( omp_get_max_threads() );
#endif // NEW_THREADS
	for( unsigned int i=0 ; i<polygons.size() ; i++ ) polygons[i] = 0;

	oocPointFile = new BufferedReadWriteFile( NULL , fileHeader );
#ifdef NEW_THREADS
	polygonFiles.resize( std::thread::hardware_concurrency() );
#else // !NEW_THREADS
	polygonFiles.resize( omp_get_max_threads() );
#endif // NEW_THREADS
	for( unsigned int i=0 ; i<polygonFiles.size() ; i++ ) polygonFiles[i] = new BufferedReadWriteFile( NULL , fileHeader );
}
template< class Vertex , typename Index >
CoredFileMeshData< Vertex , Index >::~CoredFileMeshData( void )
{
	delete oocPointFile;
	for( unsigned int i=0 ; i<polygonFiles.size() ; i++ ) delete polygonFiles[i];
}
template< class Vertex , typename Index >
void CoredFileMeshData< Vertex , Index >::resetIterator ( void )
{
	oocPointFile->reset();
	threadIndex = 0;
	for( unsigned int i=0 ; i<polygonFiles.size() ; i++ ) polygonFiles[i]->reset();
}
template< class Vertex , typename Index >
Index CoredFileMeshData< Vertex , Index >::addOutOfCorePoint( const Vertex& p )
{
	oocPointFile->write( &p , sizeof( Vertex ) );
	oocPoints++;
	return oocPoints-1;
}
template< class Vertex , typename Index >
#ifdef NEW_THREADS
Index CoredFileMeshData< Vertex , Index >::addOutOfCorePoint_s( unsigned int thread , const Vertex& p )
#else // !NEW_THREADS
Index CoredFileMeshData< Vertex , Index >::addOutOfCorePoint_s( const Vertex& p )
#endif // NEW_THREADS
{
	Index sz;
#ifdef NEW_THREADS
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock(m);
#else // !NEW_THREADS
#pragma omp critical (CoredFileMeshData_addOutOfCorePoint_s)
	{
#endif // NEW_THREADS
		sz = oocPoints;
		oocPointFile->write( &p , sizeof( Vertex ) );
		oocPoints++;
	}
	return sz;
}
template< class Vertex , typename Index >
#ifdef NEW_THREADS
void CoredFileMeshData< Vertex , Index >::addPolygon_s( unsigned int thread , const std::vector< Index >& vertices )
#else // !NEW_THREADS
void CoredFileMeshData< Vertex , Index >::addPolygon_s( const std::vector< Index >& vertices )
#endif // NEW_THREADS
{
	unsigned int vSize = (unsigned int)vertices.size();
#ifdef NEW_THREADS
#else // !NEW_THREADS
	unsigned int thread = omp_get_thread_num();
#endif // NEW_THREADS
	polygonFiles[thread]->write( &vSize , sizeof(unsigned int) );
	polygonFiles[thread]->write( &vertices[0] , sizeof(Index) * vSize );
	polygons[thread]++;
}
template< class Vertex , typename Index >
#ifdef NEW_THREADS
void CoredFileMeshData< Vertex , Index >::addPolygon_s( unsigned int thread , const std::vector< CoredVertexIndex< Index > >& vertices )
#else // !NEW_THREADS
void CoredFileMeshData< Vertex , Index >::addPolygon_s( const std::vector< CoredVertexIndex< Index > >& vertices )
#endif // NEW_THREADS
{
	std::vector< Index > polygon( vertices.size() );
	for( unsigned int i=0 ; i<(unsigned int)vertices.size() ; i++ ) 
		if( vertices[i].inCore ) polygon[i] =  vertices[i].idx;
		else                     polygon[i] = -vertices[i].idx-1;
#ifdef NEW_THREADS
	return addPolygon_s( thread , polygon );
#else // !NEW_THREADS
	return addPolygon_s( polygon );
#endif // NEW_THREADS
}
template< class Vertex , typename Index >
Index CoredFileMeshData< Vertex , Index >::nextOutOfCorePoint( Vertex& p )
{
	if( oocPointFile->read( &p , sizeof( Vertex ) ) ) return 1;
	else return 0;
}
template< class Vertex , typename Index >
Index CoredFileMeshData< Vertex , Index >::nextPolygon( std::vector< CoredVertexIndex< Index > >& vertices )
{
	while( true )
	{
		if( threadIndex<(unsigned int)polygonFiles.size() )
		{
			unsigned int pSize;
			if( polygonFiles[threadIndex]->read( &pSize , sizeof(unsigned int) ) )
			{
				std::vector< Index > polygon( pSize );
				if( polygonFiles[threadIndex]->read( &polygon[0] , sizeof(Index)*pSize ) )
				{
					vertices.resize( pSize );
					for( unsigned int i=0 ; i<(unsigned int)polygon.size() ; i++ )
						if( polygon[i]<0 ) vertices[i].idx = -polygon[i]-1 , vertices[i].inCore = false;
						else               vertices[i].idx =  polygon[i]   , vertices[i].inCore = true;
					return 1;
				}
				ERROR_OUT( "Failed to read polygon from file" );
			}
			else threadIndex++;
		}
		else return 0;
	}
}
template< class Vertex , typename Index >
size_t CoredFileMeshData< Vertex , Index >::outOfCorePointCount( void ){ return oocPoints; }
template< class Vertex , typename Index >
size_t CoredFileMeshData< Vertex , Index >::polygonCount( void )
{
	size_t count = 0;
	for( size_t i=0 ; i<polygons.size() ; i++ ) count += polygons[i];
	return count;
}
#else // !NEW_CODE
template< class Vertex >
CoredFileMeshData< Vertex >::CoredFileMeshData( const char* fileHeader )
{
	threadIndex = 0;
	oocPoints = 0;
	polygons.resize( omp_get_max_threads() );
	for( unsigned int i=0 ; i<polygons.size() ; i++ ) polygons[i] = 0;

	oocPointFile = new BufferedReadWriteFile( NULL , fileHeader );
	polygonFiles.resize( omp_get_max_threads() );
	for( unsigned int i=0 ; i<polygonFiles.size() ; i++ ) polygonFiles[i] = new BufferedReadWriteFile( NULL , fileHeader );
}
template< class Vertex >
CoredFileMeshData< Vertex >::~CoredFileMeshData( void )
{
	delete oocPointFile;
	for( unsigned int i=0 ; i<polygonFiles.size() ; i++ ) delete polygonFiles[i];
}
template< class Vertex >
void CoredFileMeshData< Vertex >::resetIterator ( void )
{
	oocPointFile->reset();
	threadIndex = 0;
	for( unsigned int i=0 ; i<polygonFiles.size() ; i++ ) polygonFiles[i]->reset();
}
template< class Vertex >
int CoredFileMeshData< Vertex >::addOutOfCorePoint( const Vertex& p )
{
	oocPointFile->write( &p , sizeof( Vertex ) );
	oocPoints++;
	return oocPoints-1;
}
template< class Vertex >
int CoredFileMeshData< Vertex >::addOutOfCorePoint_s( const Vertex& p )
{
	int sz;
#pragma omp critical (CoredFileMeshData_addOutOfCorePoint_s)
	{
		sz = oocPoints;
		oocPointFile->write( &p , sizeof( Vertex ) );
		oocPoints++;
	}
	return sz;
}
template< class Vertex >
void CoredFileMeshData< Vertex >::addPolygon_s( const std::vector< int >& vertices )
{
	int vSize = (int)vertices.size();
	int thread = omp_get_thread_num();
	polygonFiles[thread]->write( &vSize , sizeof(int) );
	polygonFiles[thread]->write( &vertices[0] , sizeof(int) * vSize );
	polygons[thread]++;
}
template< class Vertex >
void CoredFileMeshData< Vertex >::addPolygon_s( const std::vector< CoredVertexIndex >& vertices )
{
	std::vector< int > polygon( vertices.size() );
	for( int i=0 ; i<(int)vertices.size() ; i++ ) 
		if( vertices[i].inCore ) polygon[i] =  vertices[i].idx;
		else                     polygon[i] = -vertices[i].idx-1;
	return addPolygon_s( polygon );
}
template< class Vertex >
int CoredFileMeshData< Vertex >::nextOutOfCorePoint( Vertex& p )
{
	if( oocPointFile->read( &p , sizeof( Vertex ) ) ) return 1;
	else return 0;
}
template< class Vertex >
int CoredFileMeshData< Vertex >::nextPolygon( std::vector< CoredVertexIndex >& vertices )
{
	while( true )
	{
		if( threadIndex<(int)polygonFiles.size() )
		{
			int pSize;
			if( polygonFiles[threadIndex]->read( &pSize , sizeof(int) ) )
			{
				std::vector< int > polygon( pSize );
				if( polygonFiles[threadIndex]->read( &polygon[0] , sizeof(int)*pSize ) )
				{
					vertices.resize( pSize );
					for( int i=0 ; i<int(polygon.size()) ; i++ )
						if( polygon[i]<0 ) vertices[i].idx = -polygon[i]-1 , vertices[i].inCore = false;
						else               vertices[i].idx =  polygon[i]   , vertices[i].inCore = true;
						return 1;
				}
				ERROR_OUT( "Failed to read polygon from file" );
			}
			else threadIndex++;
		}
		else return 0;
	}
}
template< class Vertex >
int CoredFileMeshData< Vertex >::outOfCorePointCount( void ){ return oocPoints; }
template< class Vertex >
int CoredFileMeshData< Vertex >::polygonCount( void )
{
	int count = 0;
	for( int i=0 ; i<polygons.size() ; i++ ) count += polygons[i];
	return count;
}
#endif // NEW_CODE
/////////////
// Simplex //
/////////////
template< class Real , unsigned int Dim , unsigned int K >
void Simplex< Real , Dim , K >::split( Point< Real , Dim > pNormal , Real pOffset , std::vector< Simplex >& back , std::vector< Simplex >& front ) const
{
	Real values[K+1];
	bool frontSet = false , backSet = false;

	// Evaluate the hyper-plane's function at the vertices and mark if strictly front/back vertices have been found
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		values[k] = Point< Real , Dim >::Dot( p[k] , pNormal ) - pOffset;
		backSet |= ( values[k]<0 ) , frontSet |= ( values[k]>0 );
	}

	// If all the vertices are behind or on, or all the vertices are in front or on, we are done.
	if( !frontSet ){ back.push_back( *this ) ; return; }
	if( !backSet ){ front.push_back( *this ) ; return; }

	// Pick some intersection of the hyper-plane with a simplex edge
	unsigned int v1 , v2;
	Point< Real , Dim > midPoint;
	{
		for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=i+1 ; j<=K ; j++ ) if( values[i]*values[j]<0 )
		{
			v1 = i , v2 = j;
			Real t1 = values[i] / ( values[i] - values[j] ) , t2 = (Real)( 1. - t1 );
			midPoint = p[j]*t1 + p[i]*t2;
		}
	}
	// Iterate over each face of the simplex, split it with the hyper-plane and connect the sub-simplices to the mid-point
	for( unsigned int i=0 ; i<=K ; i++ )
	{
		if( i!=v1 && i!=v2 ) continue;
		Simplex< Real , Dim , K-1 > f;		// The face
		Simplex< Real , Dim , K > s;		// The sub-simplex
		for( unsigned int j=0 , idx=0 ; j<=K ; j++ )	if( j!=i ) f[idx++] = p[j];
		std::vector< Simplex< Real , Dim , K-1 > > _back , _front;
		f.split( pNormal , pOffset , _back , _front );
		s[i] = midPoint;

		for( unsigned int j=0 ; j<_back.size() ; j++ ) 
		{
			for( unsigned int k=0 ; k<K ; k++ ) s[ k<i ? k : k+1 ] = _back[j][k];
			back.push_back( s );
		}

		for( unsigned int j=0 ; j<_front.size() ; j++ ) 
		{
			for( unsigned int k=0 ; k<K ; k++ ) s[ k<i ? k : k+1 ] = _front[j][k];
			front.push_back( s );
		}
	}
}