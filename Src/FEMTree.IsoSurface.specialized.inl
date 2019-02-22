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

#include <sstream>
#include <iomanip>
#include <unordered_map>
#include "MyMiscellany.h"
#include "MarchingCubes.h"
#include "MAT.h"

#ifdef NEW_CODE
#define NEW_HASH
#endif // NEW_CODE


#ifdef NEW_HASH
#include <sstream>
#endif // NEW_HASH

// Specialized iso-surface extraction
template< class Real , class Vertex >
struct IsoSurfaceExtractor< 3 , Real , Vertex >
{
	static const unsigned int Dim = 3;
	typedef typename FEMTree< Dim , Real >::LocalDepth LocalDepth;
	typedef typename FEMTree< Dim , Real >::LocalOffset LocalOffset;
	typedef typename FEMTree< Dim , Real >::ConstOneRingNeighborKey ConstOneRingNeighborKey;
	typedef typename FEMTree< Dim , Real >::ConstOneRingNeighbors ConstOneRingNeighbors;
#ifdef NEW_CODE
	typedef RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > TreeNode;
#else // !NEW_CODE
	typedef RegularTreeNode< Dim , FEMTreeNodeData > TreeNode;
#endif // NEW_CODE
	template< unsigned int WeightDegree > using DensityEstimator = typename FEMTree< Dim , Real >::template DensityEstimator< WeightDegree >;
	template< typename FEMSigPack , unsigned int PointD > using _Evaluator = typename FEMTree< Dim , Real >::template _Evaluator< FEMSigPack , PointD >;
protected:
#ifdef NEW_THREADS
	static std::mutex _pointInsertionMutex;
#endif // NEW_THREADS
#ifdef NEW_HASH
	//////////
	// _Key //
	//////////
	struct _Key
	{
		int idx[Dim];

		_Key( void ){ for( unsigned int d=0 ; d<Dim ; d++ ) idx[d] = 0; }

		int &operator[]( int i ){ return idx[i]; }
		const int &operator[]( int i ) const { return idx[i]; }

		bool operator == ( const _Key &key ) const
		{
			for( unsigned int d=0 ; d<Dim ; d++ ) if( idx[d]!=key[d] ) return false;
			return true;
		}
		bool operator != ( const _Key &key ) const { return !operator==( key ); }

		std::string to_string( void ) const
		{
			std::stringstream stream;
			stream << "(";
			for( unsigned int d=0 ; d<Dim ; d++ )
			{
				if( d ) stream << ",";
				stream << idx[d];
			}
			stream << ")";
			return stream.str();
		}

		struct Hasher
		{
			size_t operator()( const _Key &i ) const
			{
				size_t hash = 0;
				for( unsigned int d=0 ; d<Dim ; d++ ) hash ^= i.idx[d];
				return hash;
			}
		};
	};
#endif // NEW_HASH

	//////////////
	// _IsoEdge //
	//////////////
	struct _IsoEdge
	{
#ifdef NEW_HASH
		_Key vertices[2];
		_IsoEdge( void ) {}
		_IsoEdge( _Key v1 , _Key v2 ){ vertices[0] = v1 , vertices[1] = v2; }
		_Key &operator[]( int idx ){ return vertices[idx]; }
		const _Key &operator[]( int idx ) const { return vertices[idx]; }
#else // !NEW_HASH
		long long edges[2];
		_IsoEdge( void ){ edges[0] = edges[1] = 0; }
		_IsoEdge( long long v1 , long long v2 ){ edges[0] = v1 , edges[1] = v2; }
		long long& operator[]( int idx ){ return edges[idx]; }
		const long long& operator[]( int idx ) const { return edges[idx]; }
#endif // NEW_HASH
	};

	////////////////
	// _FaceEdges //
	////////////////
	struct _FaceEdges{ _IsoEdge edges[2] ; int count; };

	///////////////
	// SliceData //
	///////////////
	class SliceData
	{
#ifdef NEW_CODE
		typedef RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > TreeOctNode;
#else // !NEW_CODE
		typedef RegularTreeNode< Dim , FEMTreeNodeData > TreeOctNode;
#endif // NEW_CODE
	public:
		template< unsigned int Indices >
		struct  _Indices
		{
#ifdef NEW_CODE
			node_index_type idx[Indices];
			_Indices( void ){ for( unsigned int i=0 ; i<Indices ; i++ ) idx[i] = -1; }
#else // !NEW_CODE
			int idx[Indices];
			_Indices( void ){ memset( idx , -1 , sizeof( int ) * Indices ); }
#endif // NEW_CODE
#ifdef NEW_CODE
			node_index_type& operator[] ( int i ) { return idx[i]; }
			const node_index_type& operator[] ( int i ) const { return idx[i]; }
#else // !NEW_CODE
			int& operator[] ( int i ) { return idx[i]; }
			const int& operator[] ( int i ) const { return idx[i]; }
#endif // NEW_CODE
		};
		typedef _Indices< HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() > SquareCornerIndices;
		typedef _Indices< HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() > SquareEdgeIndices;
		typedef _Indices< HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() > SquareFaceIndices;

		struct SliceTableData
		{
			Pointer( SquareCornerIndices ) cTable;
			Pointer( SquareEdgeIndices   ) eTable;
			Pointer( SquareFaceIndices   ) fTable;
#ifdef NEW_CODE
			node_index_type nodeOffset;
			node_index_type cCount , eCount , fCount;
			node_index_type nodeCount;
			SliceTableData( void ){ fCount = eCount = cCount = 0 , _oldNodeCount = 0 , cTable = NullPointer( SquareCornerIndices ) , eTable = NullPointer( SquareEdgeIndices ) , fTable = NullPointer( SquareFaceIndices ) , _cMap = _eMap = _fMap = NullPointer( node_index_type ) , _processed = NullPointer( char ); }
#else // !NEW_CODE
			int cCount , eCount , fCount , nodeOffset , nodeCount;
			SliceTableData( void ){ fCount = eCount = cCount = _oldNodeCount = 0 , cTable = NullPointer( SquareCornerIndices ) , eTable = NullPointer( SquareEdgeIndices ) , fTable = NullPointer( SquareFaceIndices ) , _cMap = _eMap = _fMap = NullPointer( int ) , _processed = NullPointer( char ); }
#endif // NEW_CODE
			void clear( void ){ DeletePointer( cTable ) ; DeletePointer( eTable ) ; DeletePointer( fTable ) ; DeletePointer( _cMap ) ; DeletePointer( _eMap ) ; DeletePointer( _fMap ) ; DeletePointer( _processed ) ; fCount = eCount = cCount = 0; }
			~SliceTableData( void ){ clear(); }

#ifdef NEW_CODE
			SquareCornerIndices& cornerIndices( const TreeOctNode* node )             { return cTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareCornerIndices& cornerIndices( node_index_type idx )                 { return cTable[ idx - nodeOffset ]; }
			const SquareCornerIndices& cornerIndices( const TreeOctNode* node ) const { return cTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareCornerIndices& cornerIndices( node_index_type idx )     const { return cTable[ idx - nodeOffset ]; }
			SquareEdgeIndices& edgeIndices( const TreeOctNode* node )                 { return eTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareEdgeIndices& edgeIndices( node_index_type idx )                     { return eTable[ idx - nodeOffset ]; }
			const SquareEdgeIndices& edgeIndices( const TreeOctNode* node )     const { return eTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareEdgeIndices& edgeIndices( node_index_type idx )         const { return eTable[ idx - nodeOffset ]; }
			SquareFaceIndices& faceIndices( const TreeOctNode* node )                 { return fTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareFaceIndices& faceIndices( node_index_type idx )                     { return fTable[ idx - nodeOffset ]; }
			const SquareFaceIndices& faceIndices( const TreeOctNode* node )     const { return fTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareFaceIndices& faceIndices( node_index_type idx )         const { return fTable[ idx - nodeOffset ]; }
#else // !NEW_CODE
			SquareCornerIndices& cornerIndices( const TreeOctNode* node ) { return cTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareCornerIndices& cornerIndices( int idx ) { return cTable[ idx - nodeOffset ]; }
			const SquareCornerIndices& cornerIndices( const TreeOctNode* node ) const { return cTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareCornerIndices& cornerIndices( int idx ) const { return cTable[ idx - nodeOffset ]; }
			SquareEdgeIndices& edgeIndices( const TreeOctNode* node ) { return eTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareEdgeIndices& edgeIndices( int idx ) { return eTable[ idx - nodeOffset ]; }
			const SquareEdgeIndices& edgeIndices( const TreeOctNode* node ) const { return eTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareEdgeIndices& edgeIndices( int idx ) const { return eTable[ idx - nodeOffset ]; }
			SquareFaceIndices& faceIndices( const TreeOctNode* node ) { return fTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareFaceIndices& faceIndices( int idx ) { return fTable[ idx - nodeOffset ]; }
			const SquareFaceIndices& faceIndices( const TreeOctNode* node ) const { return fTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareFaceIndices& faceIndices( int idx ) const { return fTable[ idx - nodeOffset ]; }
#endif // NEW_CODE

		protected:
#ifdef NEW_CODE
			Pointer( node_index_type ) _cMap;
			Pointer( node_index_type ) _eMap;
			Pointer( node_index_type ) _fMap;
			Pointer( char ) _processed;
			node_index_type _oldNodeCount;
#else // !NEW_CODE
			Pointer( int ) _cMap;
			Pointer( int ) _eMap;
			Pointer( int ) _fMap;
			Pointer( char ) _processed;
			int _oldNodeCount;
#endif // NEW_CODE
			friend SliceData;
		};
		struct XSliceTableData
		{
			Pointer( SquareCornerIndices ) eTable;
			Pointer( SquareEdgeIndices ) fTable;
#ifdef NEW_CODE
			node_index_type nodeOffset;
			node_index_type fCount , eCount;
			node_index_type nodeCount;
			XSliceTableData( void ){ fCount = eCount = 0 , _oldNodeCount = 0 , eTable = NullPointer( SquareCornerIndices ) , fTable = NullPointer( SquareEdgeIndices ) , _eMap = _fMap = NullPointer( node_index_type ); }
#else // !NEW_CODE
			int fCount , eCount , nodeOffset , nodeCount;
			XSliceTableData( void ){ fCount = eCount = _oldNodeCount = 0 , eTable = NullPointer( SquareCornerIndices ) , fTable = NullPointer( SquareEdgeIndices ) , _eMap = _fMap = NullPointer( int ); }
#endif // NEW_CODE
			~XSliceTableData( void ){ clear(); }
			void clear( void ) { DeletePointer( fTable ) ; DeletePointer( eTable ) ; DeletePointer( _eMap ) ; DeletePointer( _fMap ) ; fCount = eCount = 0; }

#ifdef NEW_CODE
			SquareCornerIndices& edgeIndices( const TreeOctNode* node )             { return eTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareCornerIndices& edgeIndices( node_index_type idx )                 { return eTable[ idx - nodeOffset ]; }
			const SquareCornerIndices& edgeIndices( const TreeOctNode* node ) const { return eTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareCornerIndices& edgeIndices( node_index_type idx )     const { return eTable[ idx - nodeOffset ]; }
			SquareEdgeIndices& faceIndices( const TreeOctNode* node )               { return fTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareEdgeIndices& faceIndices( node_index_type idx )                   { return fTable[ idx - nodeOffset ]; }
			const SquareEdgeIndices& faceIndices( const TreeOctNode* node )   const { return fTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareEdgeIndices& faceIndices( node_index_type idx )       const { return fTable[ idx - nodeOffset ]; }
#else // !NEW_CODE
			SquareCornerIndices& edgeIndices( const TreeOctNode* node ) { return eTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareCornerIndices& edgeIndices( int idx ) { return eTable[ idx - nodeOffset ]; }
			const SquareCornerIndices& edgeIndices( const TreeOctNode* node ) const { return eTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareCornerIndices& edgeIndices( int idx ) const { return eTable[ idx - nodeOffset ]; }
			SquareEdgeIndices& faceIndices( const TreeOctNode* node ) { return fTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			SquareEdgeIndices& faceIndices( int idx ) { return fTable[ idx - nodeOffset ]; }
			const SquareEdgeIndices& faceIndices( const TreeOctNode* node ) const { return fTable[ node->nodeData.nodeIndex - nodeOffset ]; }
			const SquareEdgeIndices& faceIndices( int idx ) const { return fTable[ idx - nodeOffset ]; }
#endif // NEW_CODE
		protected:
#ifdef NEW_CODE
			Pointer( node_index_type ) _eMap;
			Pointer( node_index_type ) _fMap;
			node_index_type _oldNodeCount;
#else // !NEW_CODE
			Pointer( int ) _eMap;
			Pointer( int ) _fMap;
			int _oldNodeCount;
#endif // NEW_CODE
			friend SliceData;
		};
		template< unsigned int D , unsigned int ... Ks > struct HyperCubeTables{};
		template< unsigned int D , unsigned int K >
		struct HyperCubeTables< D , K >
		{
			static unsigned int CellOffset[ HyperCube::Cube< D >::template ElementNum< K >() ][ HyperCube::Cube< D >::template IncidentCubeNum< K >() ];
			static unsigned int IncidentElementCoIndex[ HyperCube::Cube< D >::template ElementNum< K >() ][ HyperCube::Cube< D >::template IncidentCubeNum< K >() ];
			static unsigned int CellOffsetAntipodal[ HyperCube::Cube< D >::template ElementNum< K >() ];
			static typename HyperCube::Cube< D >::template IncidentCubeIndex< K > IncidentCube[ HyperCube::Cube< D >::template ElementNum< K >() ];
			static typename HyperCube::Direction Directions[ HyperCube::Cube< D >::template ElementNum< K >() ][ D ];
			static void SetTables( void )
			{
				for( typename HyperCube::Cube< D >::template Element< K > e ; e<HyperCube::Cube< D >::template ElementNum< K >() ; e++ )
				{
					for( typename HyperCube::Cube< D >::template IncidentCubeIndex< K > i ; i<HyperCube::Cube< D >::template IncidentCubeNum< K >() ; i++ )
					{
						CellOffset[e.index][i.index] = HyperCube::Cube< D >::CellOffset( e , i );
						IncidentElementCoIndex[e.index][i.index] = HyperCube::Cube< D >::IncidentElement( e , i ).coIndex();
					}
					CellOffsetAntipodal[e.index] = HyperCube::Cube< D >::CellOffset( e , HyperCube::Cube< D >::IncidentCube( e ).antipodal() );
					IncidentCube[ e.index ] = HyperCube::Cube< D >::IncidentCube( e );
					e.directions( Directions[e.index] );
				}
			}
		};
		template< unsigned int D , unsigned int K1 , unsigned int K2 >
		struct HyperCubeTables< D , K1 , K2 >
		{
			static typename HyperCube::Cube< D >::template Element< K2 > OverlapElements[ HyperCube::Cube< D >::template ElementNum< K1 >() ][ HyperCube::Cube< D >::template OverlapElementNum< K1 , K2 >() ];
			static bool Overlap[ HyperCube::Cube< D >::template ElementNum< K1 >() ][ HyperCube::Cube< D >::template ElementNum< K2 >() ];
			static void SetTables( void )
			{
				for( typename HyperCube::Cube< D >::template Element< K1 > e ; e<HyperCube::Cube< D >::template ElementNum< K1 >() ; e++ )
				{
					for( typename HyperCube::Cube< D >::template Element< K2 > _e ; _e<HyperCube::Cube< D >::template ElementNum< K2 >() ; _e++ )
						Overlap[e.index][_e.index] = HyperCube::Cube< D >::Overlap( e , _e );
					HyperCube::Cube< D >::OverlapElements( e , OverlapElements[e.index] );
				}
				if( !K2 ) HyperCubeTables< D , K1 >::SetTables();
			}
		};

		template< unsigned int D=Dim , unsigned int K1=Dim , unsigned int K2=Dim > static typename std::enable_if<                 K2!=0 >::type SetHyperCubeTables( void )
		{
			HyperCubeTables< D , K1 , K2 >::SetTables() ; SetHyperCubeTables< D , K1 , K2-1 >();
		}
		template< unsigned int D=Dim , unsigned int K1=Dim , unsigned int K2=Dim > static typename std::enable_if<        K1!=0 && K2==0 >::type SetHyperCubeTables( void )
		{
			HyperCubeTables< D , K1 , K2 >::SetTables(); SetHyperCubeTables< D , K1-1 , D >();
		}
		template< unsigned int D=Dim , unsigned int K1=Dim , unsigned int K2=Dim > static typename std::enable_if< D!=1 && K1==0 && K2==0 >::type SetHyperCubeTables( void )
		{
			HyperCubeTables< D , K1 , K2 >::SetTables() ; SetHyperCubeTables< D-1 , D-1 , D-1 >();
		}
		template< unsigned int D=Dim , unsigned int K1=Dim , unsigned int K2=Dim > static typename std::enable_if< D==1 && K1==0 && K2==0 >::type SetHyperCubeTables( void )
		{
			HyperCubeTables< D , K1 , K2 >::SetTables();
		}

#ifdef NEW_THREADS
		static void SetSliceTableData( ThreadPool &tp , const SortedTreeNodes< Dim >& sNodes , SliceTableData* sData0 , XSliceTableData* xData , SliceTableData* sData1 , int depth , int offset )
#else // !NEW_THREADS
		static void SetSliceTableData( const SortedTreeNodes< Dim >& sNodes , SliceTableData* sData0 , XSliceTableData* xData , SliceTableData* sData1 , int depth , int offset )
#endif // NEW_THREADS
		{
			// [NOTE] This is structure is purely for determining adjacency and is independent of the FEM degree
			typedef typename FEMTree< Dim , Real >::ConstOneRingNeighborKey ConstOneRingNeighborKey;
			if( offset<0 || offset>((size_t)1<<depth) ) return;
			if( sData0 )
			{
#ifdef NEW_CODE
				std::pair< node_index_type , node_index_type > span( sNodes.begin( depth , offset-1 ) , sNodes.end( depth , offset ) );
#else // !NEW_CODE
				std::pair< int , int > span( sNodes.begin( depth , offset-1 ) , sNodes.end( depth , offset ) );
#endif // NEW_CODE
				sData0->nodeOffset = span.first , sData0->nodeCount = span.second - span.first;
			}
			if( sData1 )
			{
#ifdef NEW_CODE
				std::pair< node_index_type , node_index_type > span( sNodes.begin( depth , offset ) , sNodes.end( depth , offset+1 ) );
#else // !NEW_CODE
				std::pair< int , int > span( sNodes.begin( depth , offset ) , sNodes.end( depth , offset+1 ) );
#endif // NEW_CODE
				sData1->nodeOffset = span.first , sData1->nodeCount = span.second - span.first;
			}
			if( xData )
			{
#ifdef NEW_CODE
				std::pair< node_index_type , node_index_type > span( sNodes.begin( depth , offset ) , sNodes.end( depth , offset ) );
#else // !NEW_CODE
				std::pair< int , int > span( sNodes.begin( depth , offset ) , sNodes.end( depth , offset ) );
#endif // NEW_CODE
				xData->nodeOffset = span.first , xData->nodeCount = span.second - span.first;
			}
			SliceTableData* sData[] = { sData0 , sData1 };
			for( int i=0 ; i<2 ; i++ ) if( sData[i] )
			{
				if( sData[i]->nodeCount>sData[i]->_oldNodeCount )
				{
					DeletePointer( sData[i]->_cMap ) ; DeletePointer( sData[i]->_eMap ) ; DeletePointer( sData[i]->_fMap );
					DeletePointer( sData[i]->cTable ) ; DeletePointer( sData[i]->eTable ) ; DeletePointer( sData[i]->fTable );
					DeletePointer( sData[i]->_processed );
#ifdef NEW_CODE
					sData[i]->_cMap = NewPointer< node_index_type >( sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() );
					sData[i]->_eMap = NewPointer< node_index_type >( sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() );
					sData[i]->_fMap = NewPointer< node_index_type >( sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() );
#else // !NEW_CODE
					sData[i]->_cMap = NewPointer< int >( sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() );
					sData[i]->_eMap = NewPointer< int >( sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() );
					sData[i]->_fMap = NewPointer< int >( sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() );
#endif // NEW_CODE
					sData[i]->_processed = NewPointer< char >( sData[i]->nodeCount );
					sData[i]->cTable = NewPointer< typename SliceData::SquareCornerIndices >( sData[i]->nodeCount );
					sData[i]->eTable = NewPointer< typename SliceData::SquareEdgeIndices >( sData[i]->nodeCount );
					sData[i]->fTable = NewPointer< typename SliceData::SquareFaceIndices >( sData[i]->nodeCount );
					sData[i]->_oldNodeCount = sData[i]->nodeCount;
				}
#ifdef NEW_CODE
				memset( sData[i]->_cMap , 0 , sizeof(node_index_type) * sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() );
				memset( sData[i]->_eMap , 0 , sizeof(node_index_type) * sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() );
				memset( sData[i]->_fMap , 0 , sizeof(node_index_type) * sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() );
#else // !NEW_CODE
				memset( sData[i]->_cMap , 0 , sizeof(int) * sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() );
				memset( sData[i]->_eMap , 0 , sizeof(int) * sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() );
				memset( sData[i]->_fMap , 0 , sizeof(int) * sData[i]->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() );
#endif // NEW_CODE
				memset( sData[i]->_processed , 0 , sizeof(char) * sData[i]->nodeCount );
			}
			if( xData )
			{
				if( xData->nodeCount>xData->_oldNodeCount )
				{
					DeletePointer( xData->_eMap ) ; DeletePointer( xData->_fMap );
					DeletePointer( xData->eTable ) ; DeletePointer( xData->fTable );
#ifdef NEW_CODE
					xData->_eMap = NewPointer< node_index_type >( xData->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() );
					xData->_fMap = NewPointer< node_index_type >( xData->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() );
#else // !NEW_CODE
					xData->_eMap = NewPointer< int >( xData->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() );
					xData->_fMap = NewPointer< int >( xData->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() );
#endif // NEW_CODE
					xData->eTable = NewPointer< typename SliceData::SquareCornerIndices >( xData->nodeCount );
					xData->fTable = NewPointer< typename SliceData::SquareEdgeIndices >( xData->nodeCount );
					xData->_oldNodeCount = xData->nodeCount;
				}
#ifdef NEW_CODE
				memset( xData->_eMap , 0 , sizeof(node_index_type) * xData->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() );
				memset( xData->_fMap , 0 , sizeof(node_index_type) * xData->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() );
#else // !NEW_CODE
				memset( xData->_eMap , 0 , sizeof(int) * xData->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() );
				memset( xData->_fMap , 0 , sizeof(int) * xData->nodeCount * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() );
#endif // NEW_CODE
			}
#ifdef NEW_THREADS
			std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
			std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
			for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth );

			typedef typename FEMTree< Dim , Real >::ConstOneRingNeighbors ConstNeighbors;

			// Process the corners
			// z: which side of the cell	\in {0,1}
			// zOff: which neighbor			\in {-1,0,1}
			auto ProcessCorners = []( SliceTableData& sData , const ConstNeighbors& neighbors , HyperCube::Direction zDir , int zOff )
			{
				const TreeOctNode* node = neighbors.neighbors[1][1][1+zOff];
#ifdef NEW_CODE
				node_index_type i = node->nodeData.nodeIndex;
#else // !NEW_CODE
				int i = node->nodeData.nodeIndex;
#endif // NEW_CODE
				// Iterate over the corners in the face
				for( typename HyperCube::Cube< Dim-1 >::template Element< 0 > _c ; _c<HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; _c++ )
				{
					bool owner = true;

					typename HyperCube::Cube< Dim >::template Element< 0 > c( zDir , _c.index );																	// Corner-in-cube index
					typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 0 > my_ic = HyperCubeTables< Dim , 0 >::IncidentCube[c.index];						// The index of the node relative to the corner
					for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 0 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 0 >() ; ic++ )	// Iterate over the nodes adjacent to the corner
					{
						// Get the index of cube relative to the corner neighbors
						unsigned int xx = HyperCubeTables< Dim , 0 >::CellOffset[c.index][ic.index] + zOff;
						// If the neighbor exists and comes before, they own the corner
						if( neighbors.neighbors.data[xx] && ic<my_ic ){ owner = false ; break; }
					}
					if( owner )
					{
#ifdef NEW_CODE
						node_index_type myCount = ( i-sData.nodeOffset ) * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() + _c.index;
#else // !NEW_CODE
						int myCount = (i - sData.nodeOffset) * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() + _c.index;
#endif // NEW_CODE
						sData._cMap[ myCount ] = 1;
						// Set the corner pointer for all cubes incident on the corner
						for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 0 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 0 >() ; ic++ )	// Iterate over the nodes adjacent to the corner
						{
							unsigned int xx = HyperCubeTables< Dim , 0 >::CellOffset[c.index][ic.index] + zOff;
							// If the neighbor exits, sets its corner
							if( neighbors.neighbors.data[xx] ) sData.cornerIndices( neighbors.neighbors.data[xx] )[ HyperCubeTables< Dim , 0 >::IncidentElementCoIndex[c.index][ic.index] ] = myCount;
						}
					}
				}
			};
			// Process the in-plane edges
			auto ProcessIEdges = []( SliceTableData& sData , const ConstNeighbors& neighbors , HyperCube::Direction zDir , int zOff )
			{
				const TreeOctNode* node = neighbors.neighbors[1][1][1+zOff];
#ifdef NEW_CODE
				node_index_type i = node->nodeData.nodeIndex;
#else // !NEW_CODE
				int i = node->nodeData.nodeIndex;
#endif // NEW_CODE
				// Iterate over the edges in the face
				for( typename HyperCube::Cube< Dim-1 >::template Element< 1 > _e ; _e<HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; _e++ )
				{
					bool owner = true;

					// The edge in the cube
					typename HyperCube::Cube< Dim >::template Element< 1 > e( zDir , _e.index );
					// The index of the cube relative to the edge
					typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > my_ic = HyperCubeTables< Dim , 1 >::IncidentCube[e.index];
					// Iterate over the cubes incident on the edge
					for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 1 >() ; ic++ )
					{
						// Get the indices of the cube relative to the center
						unsigned int xx = HyperCubeTables< Dim , 1 >::CellOffset[e.index][ic.index] + zOff;
						// If the neighbor exists and comes before, they own the corner
						if( neighbors.neighbors.data[xx] && ic<my_ic ){ owner = false ; break; }
					}
					if( owner )
					{
#ifdef NEW_CODE
						node_index_type myCount = ( i - sData.nodeOffset ) * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() + _e.index;
#else // !NEW_CODE
						int myCount = ( i - sData.nodeOffset ) * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() + _e.index;
#endif // NEW_CODE
						sData._eMap[ myCount ] = 1;
						// Set all edge indices
						for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 1 >() ; ic++ )
						{
							unsigned int xx = HyperCubeTables< Dim , 1 >::CellOffset[e.index][ic.index] + zOff;
							// If the neighbor exists, set the index
							if( neighbors.neighbors.data[xx] ) sData.edgeIndices( neighbors.neighbors.data[xx] )[ HyperCubeTables< Dim , 1 >::IncidentElementCoIndex[e.index][ic.index] ] = myCount;
						}
					}
				}
			};
			// Process the cross-plane edges
			auto ProcessXEdges = []( XSliceTableData& xData , const ConstNeighbors& neighbors )
			{
				const TreeOctNode* node = neighbors.neighbors[1][1][1];
#ifdef NEW_CODE
				node_index_type i = node->nodeData.nodeIndex;
#else // !NEW_CODE
				int i = node->nodeData.nodeIndex;
#endif // NEW_CODE
				for( typename HyperCube::Cube< Dim-1 >::template Element< 0 > _c ; _c<HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; _c++ )
				{
					bool owner = true;

					typename HyperCube::Cube< Dim >::template Element< 1 > e( HyperCube::CROSS , _c.index );
					typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > my_ic = HyperCubeTables< Dim , 1 >::IncidentCube[e.index];

					for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 1 >() ; ic++ )
					{
						unsigned int xx = HyperCubeTables< Dim , 1 >::CellOffset[e.index][ic.index];
						if( neighbors.neighbors.data[xx] && ic<my_ic ){ owner = false ; break; }
					}
					if( owner )
					{
#ifdef NEW_CODE
						node_index_type myCount = ( i - xData.nodeOffset ) * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() + _c.index;
#else // !NEW_CODE
						int myCount = ( i - xData.nodeOffset ) * HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() + _c.index;
#endif // NEW_CODE
						xData._eMap[ myCount ] = 1;

						// Set all edge indices
						for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 1 >() ; ic++ )
						{
							unsigned int xx = HyperCubeTables< Dim , 1 >::CellOffset[e.index][ic.index];
							if( neighbors.neighbors.data[xx] ) xData.edgeIndices( neighbors.neighbors.data[xx] )[ HyperCubeTables< Dim , 1 >::IncidentElementCoIndex[e.index][ic.index] ] = myCount;
						}
					}
				}
			};
			// Process the in-plane faces
			auto ProcessIFaces = []( SliceTableData& sData , const ConstNeighbors& neighbors , HyperCube::Direction zDir , int zOff )
			{
				const TreeOctNode* node = neighbors.neighbors[1][1][1+zOff];
#ifdef NEW_CODE
				node_index_type i = node->nodeData.nodeIndex;
#else // !NEW_CODE
				int i = node->nodeData.nodeIndex;
#endif // NEW_CODE
				for( typename HyperCube::Cube< Dim-1 >::template Element< 2 > _f ; _f<HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() ; _f++ )
				{
					bool owner = true;

					typename HyperCube::Cube< Dim >::template Element< 2 > f( zDir , _f.index );				
					typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 2 > my_ic = HyperCubeTables< Dim , 2 >::IncidentCube[f.index];

					for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 2 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 2 >() ; ic++ )
					{
						unsigned int xx = HyperCubeTables< Dim , 2 >::CellOffset[f.index][ic.index] + zOff;
						if( neighbors.neighbors.data[xx] && ic<my_ic ){ owner = false ; break; }
					}
					if( owner )
					{
#ifdef NEW_CODE
						node_index_type myCount = ( i - sData.nodeOffset ) * HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() + _f.index;
#else // !NEW_CODE
						int myCount = ( i - sData.nodeOffset ) * HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() + _f.index;
#endif // NEW_CODE
						sData._fMap[ myCount ] = 1;

						// Set all the face indices
						for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 2 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 2 >() ; ic++ )
						{
							unsigned int xx = HyperCubeTables< Dim , 2 >::CellOffset[f.index][ic.index] + zOff;
							if( neighbors.neighbors.data[xx] ) sData.faceIndices( neighbors.neighbors.data[xx] )[ HyperCubeTables< Dim , 2 >::IncidentElementCoIndex[f.index][ic.index] ] = myCount;
						}
					}
				}
			};

			// Process the cross-plane faces
			auto ProcessXFaces = []( XSliceTableData& xData , const ConstNeighbors& neighbors )
			{
				const TreeOctNode* node = neighbors.neighbors[1][1][1];
#ifdef NEW_CODE
				node_index_type i = node->nodeData.nodeIndex;
#else // !NEW_CODE
				int i = node->nodeData.nodeIndex;
#endif // NEW_CODE
				for( typename HyperCube::Cube< Dim-1 >::template Element< 1 > _e ; _e<HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; _e++ )
				{
					bool owner = true;

					typename HyperCube::Cube< Dim >::template Element< 2 > f( HyperCube::CROSS , _e.index );				
					typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 2 > my_ic = HyperCubeTables< Dim , 2 >::IncidentCube[f.index];

					for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 2 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 2 >() ; ic++ )
					{
						unsigned int xx = HyperCubeTables< Dim , 2 >::CellOffset[f.index][ic.index];
						if( neighbors.neighbors.data[xx] && ic<my_ic ){ owner = false ; break; }
					}
					if( owner )
					{
#ifdef NEW_CODE
						node_index_type myCount = ( i - xData.nodeOffset ) * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() + _e.index;
#else // !NEW_CODE
						int myCount = ( i - xData.nodeOffset ) * HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() + _e.index;
#endif // NEW_CODE
						xData._fMap[ myCount ] = 1;

						// Set all the face indices
						for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 2 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 2 >() ; ic++ )
						{
							unsigned int xx = HyperCubeTables< Dim , 2 >::CellOffset[f.index][ic.index];
							if( neighbors.neighbors.data[xx] ) xData.faceIndices( neighbors.neighbors.data[xx] )[ HyperCubeTables< Dim , 2 >::IncidentElementCoIndex[f.index][ic.index] ] = myCount;
						}
					}
				}
			};

			// Try and get at the nodes outside of the slab through the neighbor key
#ifdef NEW_THREADS
			tp.parallel_for( sNodes.begin(depth,offset) , sNodes.end(depth,offset) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for schedule( guided )
#ifdef NEW_CODE
			for( node_index_type i=sNodes.begin(depth,offset) ; i<sNodes.end(depth,offset) ; i++ )
#else // !NEW_CODE
			for( int i=sNodes.begin(depth,offset) ; i<sNodes.end(depth,offset) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
			{
#ifdef NEW_THREADS
				ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#else // !NEW_THREADS
				ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
				const TreeOctNode* node = sNodes.treeNodes[i];
				ConstNeighbors& neighbors = neighborKey.getNeighbors( node );
				for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) if( !IsActiveNode< Dim >( neighbors.neighbors[i][j][k] ) ) neighbors.neighbors[i][j][k] = NULL;

				if( sData0 )
				{
					ProcessCorners( *sData0 , neighbors , HyperCube::BACK , 0 ) , ProcessIEdges( *sData0 , neighbors , HyperCube::BACK , 0 ) , ProcessIFaces( *sData0 , neighbors , HyperCube::BACK , 0 );
					const TreeOctNode* _node = neighbors.neighbors[1][1][0];
					if( _node )
					{
						ProcessCorners( *sData0 , neighbors , HyperCube::FRONT , -1 ) , ProcessIEdges( *sData0 , neighbors , HyperCube::FRONT , -1 ) , ProcessIFaces( *sData0 , neighbors , HyperCube::FRONT , -1 );
						sData0->_processed[ _node->nodeData.nodeIndex - sNodes.begin(depth,offset-1) ] = 1;
					}
				}
				if( sData1 )
				{
					ProcessCorners( *sData1 , neighbors , HyperCube::FRONT , 0 ) , ProcessIEdges( *sData1 , neighbors , HyperCube::FRONT , 0 ) , ProcessIFaces( *sData1 , neighbors , HyperCube::FRONT , 0 );
					const TreeOctNode* _node = neighbors.neighbors[1][1][2];
					if( _node )
					{
						ProcessCorners( *sData1 , neighbors , HyperCube::BACK , 1 ) , ProcessIEdges( *sData1 , neighbors , HyperCube::BACK , 1 ) , ProcessIFaces( *sData1, neighbors , HyperCube::BACK , 1 );
						sData1->_processed[ _node->nodeData.nodeIndex - sNodes.begin(depth,offset+1) ] = true;
					}
				}
				if( xData ) ProcessXEdges( *xData , neighbors ) , ProcessXFaces( *xData , neighbors );
			}
#ifdef NEW_THREADS
			);
#endif // NEW_THREADS
			if( sData0 )
			{
#ifdef NEW_CODE
				node_index_type off = sNodes.begin(depth,offset-1);
				node_index_type size = sNodes.end(depth,offset-1) - sNodes.begin(depth,offset-1);
#else // !NEW_CODE
				int off = sNodes.begin(depth,offset-1) , size = sNodes.end(depth,offset-1) - sNodes.begin(depth,offset-1);
#endif // NEW_CODE
#ifdef NEW_THREADS
				tp.parallel_for( 0 , size , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for schedule( guided )
#ifdef NEW_CODE
				for( node_index_type i=0 ; i<size ; i++ )
#else // !NEW_CODE
				for( int i=0 ; i<size ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
				{
					if( !sData0->_processed[i] )
					{
#ifdef NEW_THREADS
						ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#else // !NEW_THREADS
						ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
						const TreeOctNode* node = sNodes.treeNodes[i+off];
						ConstNeighbors& neighbors = neighborKey.getNeighbors( node );
						for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) if( !IsActiveNode< Dim >( neighbors.neighbors[i][j][k] ) ) neighbors.neighbors[i][j][k] = NULL;
						ProcessCorners( *sData0 , neighbors , HyperCube::FRONT , 0 ) , ProcessIEdges( *sData0 , neighbors , HyperCube::FRONT , 0 ) , ProcessIFaces( *sData0 , neighbors , HyperCube::FRONT , 0 );
					}
				}
#ifdef NEW_THREADS
				);
#endif // NEW_THREADS
			}
			if( sData1 )
			{
#ifdef NEW_CODE
				node_index_type off = sNodes.begin(depth,offset+1);
				node_index_type size = sNodes.end(depth,offset+1) - sNodes.begin(depth,offset+1);
#else // !NEW_CODE
				int off = sNodes.begin(depth,offset+1) , size = sNodes.end(depth,offset+1) - sNodes.begin(depth,offset+1);
#endif // NEW_CODE
#ifdef NEW_THREADS
				tp.parallel_for( 0 , size , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for schedule( guided )
#ifdef NEW_CODE
				for( node_index_type i=0 ; i<size ; i++ )
#else // !NEW_CODE
				for( int i=0 ; i<size ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
				{
					if( !sData1->_processed[i] )
					{
#ifdef NEW_THREADS
						ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#else // !NEW_THREADS
						ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
						const TreeOctNode* node = sNodes.treeNodes[i+off];
						ConstNeighbors& neighbors = neighborKey.getNeighbors( node );
						for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) if( !IsActiveNode< Dim >( neighbors.neighbors[i][j][k] ) ) neighbors.neighbors[i][j][k] = NULL;
						ProcessCorners( *sData1 , neighbors , HyperCube::BACK , 0 ) , ProcessIEdges( *sData1 , neighbors , HyperCube::BACK , 0 ) , ProcessIFaces( *sData1 , neighbors , HyperCube::BACK , 0 );
					}
				}
#ifdef NEW_THREADS
				);
#endif // NEW_THREADS
			}

			auto SetICounts = [&]( SliceTableData& sData )
			{
#ifdef NEW_CODE
				node_index_type cCount = 0 , eCount = 0 , fCount = 0;

				for( node_index_type i=0 ; i<sData.nodeCount * (node_index_type)HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; i++ ) if( sData._cMap[i] ) sData._cMap[i] = cCount++;
				for( node_index_type i=0 ; i<sData.nodeCount * (node_index_type)HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; i++ ) if( sData._eMap[i] ) sData._eMap[i] = eCount++;
				for( node_index_type i=0 ; i<sData.nodeCount * (node_index_type)HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() ; i++ ) if( sData._fMap[i] ) sData._fMap[i] = fCount++;
#ifdef NEW_THREADS
				tp.parallel_for( 0 , sData.nodeCount , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
				for( node_index_type i=0 ; i<sData.nodeCount ; i++ )
#endif // NEW_THREADS
				{
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; j++ ) sData.cTable[i][j] = sData._cMap[ sData.cTable[i][j] ];
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; j++ ) sData.eTable[i][j] = sData._eMap[ sData.eTable[i][j] ];
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() ; j++ ) sData.fTable[i][j] = sData._fMap[ sData.fTable[i][j] ];
				}
#ifdef NEW_THREADS
				);
#endif // NEW_THREADS
#else // !NEW_CODE
				int cCount = 0 , eCount = 0 , fCount = 0;

				for( int i=0 ; i<sData.nodeCount * (int)HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; i++ ) if( sData._cMap[i] ) sData._cMap[i] = cCount++;
				for( int i=0 ; i<sData.nodeCount * (int)HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; i++ ) if( sData._eMap[i] ) sData._eMap[i] = eCount++;
				for( int i=0 ; i<sData.nodeCount * (int)HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() ; i++ ) if( sData._fMap[i] ) sData._fMap[i] = fCount++;
#pragma omp parallel for
				for( int i=0 ; i<sData.nodeCount ; i++ )
				{
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; j++ ) sData.cTable[i][j] = sData._cMap[ sData.cTable[i][j] ];
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; j++ ) sData.eTable[i][j] = sData._eMap[ sData.eTable[i][j] ];
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 2 >() ; j++ ) sData.fTable[i][j] = sData._fMap[ sData.fTable[i][j] ];
				}
#endif // NEW_CODE
				sData.cCount = cCount , sData.eCount = eCount , sData.fCount = fCount;
			};
			auto SetXCounts = [&]( XSliceTableData& xData )
			{
#ifdef NEW_CODE
				node_index_type eCount = 0 , fCount = 0;

				for( node_index_type i=0 ; i<xData.nodeCount * (node_index_type)HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; i++ ) if( xData._eMap[i] ) xData._eMap[i] = eCount++;
				for( node_index_type i=0 ; i<xData.nodeCount * (node_index_type)HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; i++ ) if( xData._fMap[i] ) xData._fMap[i] = fCount++;
#ifdef NEW_THREADS
				tp.parallel_for( 0 , xData.nodeCount , [&]( unsigned int , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
				for( node_index_type i=0 ; i<xData.nodeCount ; i++ )
#endif // NEW_THREADS
				{
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; j++ ) xData.eTable[i][j] = xData._eMap[ xData.eTable[i][j] ];
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; j++ ) xData.fTable[i][j] = xData._fMap[ xData.fTable[i][j] ];
				}
#ifdef NEW_THREADS
				);
#endif // NEW_THREADS
#else // !NEW_CODE
				int eCount = 0 , fCount = 0;

				for( int i=0 ; i<xData.nodeCount * (int)HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; i++ ) if( xData._eMap[i] ) xData._eMap[i] = eCount++;
				for( int i=0 ; i<xData.nodeCount * (int)HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; i++ ) if( xData._fMap[i] ) xData._fMap[i] = fCount++;
#pragma omp parallel for
				for( int i=0 ; i<xData.nodeCount ; i++ )
				{
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; j++ ) xData.eTable[i][j] = xData._eMap[ xData.eTable[i][j] ];
					for( unsigned int j=0 ; j<HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; j++ ) xData.fTable[i][j] = xData._fMap[ xData.fTable[i][j] ];
				}
#endif // NEW_CODE
				xData.eCount = eCount , xData.fCount = fCount;
			};

			if( sData0 ) SetICounts( *sData0 );
			if( sData1 ) SetICounts( *sData1 );
			if( xData  ) SetXCounts( *xData  );
		}
	};

	//////////////////
	// _SliceValues //
	//////////////////
	struct _SliceValues
	{
		typename SliceData::SliceTableData sliceData;
		Pointer( Real ) cornerValues ; Pointer( Point< Real , Dim > ) cornerGradients ; Pointer( char ) cornerSet;
#ifdef NEW_HASH
		Pointer( _Key ) edgeKeys ; Pointer( char ) edgeSet;
#else // !NEW_HASH
		Pointer( long long ) edgeKeys ; Pointer( char ) edgeSet;
#endif // NEW_HASH
		Pointer( _FaceEdges ) faceEdges ; Pointer( char ) faceSet;
		Pointer( char ) mcIndices;
#ifdef NEW_HASH
		std::unordered_map< _Key , std::vector< _IsoEdge > , typename _Key::Hasher > faceEdgeMap;
#ifdef NEW_CODE
		std::unordered_map< _Key , std::pair< node_index_type, Vertex > , typename _Key::Hasher > edgeVertexMap;
#else // !NEW_CODE
		std::unordered_map< _Key , std::pair< int , Vertex > , typename _Key::Hasher > edgeVertexMap;
#endif // NEW_CODE
		std::unordered_map< _Key , _Key , typename _Key::Hasher > vertexPairMap;
		std::vector< std::vector< std::pair< _Key , std::vector< _IsoEdge > > > > faceEdgeKeyValues;
#ifdef NEW_CODE
		std::vector< std::vector< std::pair< _Key , std::pair< node_index_type , Vertex > > > > edgeVertexKeyValues;
#else // !NEW_CODE
		std::vector< std::vector< std::pair< _Key , std::pair< int , Vertex > > > > edgeVertexKeyValues;
#endif // NEW_CODE
		std::vector< std::vector< std::pair< _Key , _Key > > > vertexPairKeyValues;
#else // !NEW_HASH
		std::unordered_map< long long , std::vector< _IsoEdge > > faceEdgeMap;
#ifdef NEW_CODE
		std::unordered_map< long long , std::pair< node_index_type, Vertex > > edgeVertexMap;
#else // !NEW_CODE
		std::unordered_map< long long , std::pair< int, Vertex > > edgeVertexMap;
#endif // NEW_CODE
		std::unordered_map< long long , long long > vertexPairMap;
		std::vector< std::vector< std::pair< long long , std::vector< _IsoEdge > > > > faceEdgeKeyValues;
#ifdef NEW_CODE
		std::vector< std::vector< std::pair< long long , std::pair< node_index_type , Vertex > > > > edgeVertexKeyValues;
#else // !NEW_CODE
		std::vector< std::vector< std::pair< long long , std::pair< int , Vertex > > > > edgeVertexKeyValues;
#endif // NEW_CODE
		std::vector< std::vector< std::pair< long long , long long > > > vertexPairKeyValues;
#endif // NEW_HASH

#ifdef NEW_THREADS
		_SliceValues( ThreadPool &tp )
#else // !NEW_THREADS
		_SliceValues( void )
#endif // NEW_THREADS
		{
#ifdef NEW_CODE
			_oldCCount = _oldECount = _oldFCount = 0;
			_oldNCount = 0;
#else // !NEW_CODE
			_oldCCount = _oldECount = _oldFCount = _oldNCount = 0;
#endif // NEW_CODE
			cornerValues = NullPointer( Real ) ; cornerGradients = NullPointer( Point< Real , Dim > ) ; cornerSet = NullPointer( char );
#ifdef NEW_HASH
			edgeKeys = NullPointer( _Key ) ; edgeSet = NullPointer( char );
#else // !NEW_HASH
			edgeKeys = NullPointer( long long ) ; edgeSet = NullPointer( char );
#endif // NEW_HASH
			faceEdges = NullPointer( _FaceEdges ) ; faceSet = NullPointer( char );
			mcIndices = NullPointer( char );
#ifdef NEW_THREADS
			edgeVertexKeyValues.resize( tp.threadNum() );
			vertexPairKeyValues.resize( tp.threadNum() );
			faceEdgeKeyValues.resize( tp.threadNum() );
#else // !NEW_THREADS
			edgeVertexKeyValues.resize( omp_get_max_threads() );
			vertexPairKeyValues.resize( omp_get_max_threads() );
			faceEdgeKeyValues.resize( omp_get_max_threads() );
#endif // NEW_THREADS
		}
		~_SliceValues( void )
		{
#ifdef NEW_CODE
			_oldCCount = _oldECount = _oldFCount = 0;
			_oldNCount = 0;
#else // !NEW_CODE
			_oldCCount = _oldECount = _oldFCount = _oldNCount = 0;
#endif // NEW_CODE
			FreePointer( cornerValues ) ; FreePointer( cornerGradients ) ; FreePointer( cornerSet );
			FreePointer( edgeKeys ) ; FreePointer( edgeSet );
			FreePointer( faceEdges ) ; FreePointer( faceSet );
			FreePointer( mcIndices );
		}
		void setEdgeVertexMap( void )
		{
#ifdef NEW_CODE
			for( node_index_type i=0 ; i<(node_index_type)edgeVertexKeyValues.size() ; i++ )
#else // !NEW_CODE
			for( int i=0 ; i<edgeVertexKeyValues.size() ; i++ )
#endif // NEW_CODE
			{
				for( int j=0 ; j<edgeVertexKeyValues[i].size() ; j++ ) edgeVertexMap[ edgeVertexKeyValues[i][j].first ] = edgeVertexKeyValues[i][j].second;
				edgeVertexKeyValues[i].clear();
			}
		}
		void setVertexPairMap( void )
		{
#ifdef NEW_CODE
			for( node_index_type i=0 ; i<(node_index_type)vertexPairKeyValues.size() ; i++ )
#else // !NEW_CODE
			for( int i=0 ; i<vertexPairKeyValues.size() ; i++ )
#endif // NEW_CODE
			{
				for( int j=0 ; j<vertexPairKeyValues[i].size() ; j++ )
				{
					vertexPairMap[ vertexPairKeyValues[i][j].first ] = vertexPairKeyValues[i][j].second;
					vertexPairMap[ vertexPairKeyValues[i][j].second ] = vertexPairKeyValues[i][j].first;
				}
				vertexPairKeyValues[i].clear();
			}
		}
		void setFaceEdgeMap( void )
		{
#ifdef NEW_CODE
			for( node_index_type i=0 ; i<(node_index_type)faceEdgeKeyValues.size() ; i++ )
#else // !NEW_CODE
			for( int i=0 ; i<faceEdgeKeyValues.size() ; i++ )
#endif // NEW_CODE
			{
				for( int j=0 ; j<faceEdgeKeyValues[i].size() ; j++ )
				{
					auto iter = faceEdgeMap.find( faceEdgeKeyValues[i][j].first );
					if( iter==faceEdgeMap.end() ) faceEdgeMap[ faceEdgeKeyValues[i][j].first ] = faceEdgeKeyValues[i][j].second;
					else for( int k=0 ; k<faceEdgeKeyValues[i][j].second.size() ; k++ ) iter->second.push_back( faceEdgeKeyValues[i][j].second[k] );
				}
				faceEdgeKeyValues[i].clear();
			}
		}
		void reset( bool nonLinearFit )
		{
			faceEdgeMap.clear() , edgeVertexMap.clear() , vertexPairMap.clear();
#ifdef NEW_CODE
			for( node_index_type i=0 ; i<(node_index_type)edgeVertexKeyValues.size() ; i++ ) edgeVertexKeyValues[i].clear();
			for( node_index_type i=0 ; i<(node_index_type)vertexPairKeyValues.size() ; i++ ) vertexPairKeyValues[i].clear();
			for( node_index_type i=0 ; i<(node_index_type)faceEdgeKeyValues.size() ; i++ ) faceEdgeKeyValues[i].clear();
#else // !NEW_CODE
			for( int i=0 ; i<edgeVertexKeyValues.size() ; i++ ) edgeVertexKeyValues[i].clear();
			for( int i=0 ; i<vertexPairKeyValues.size() ; i++ ) vertexPairKeyValues[i].clear();
			for( int i=0 ; i<faceEdgeKeyValues.size() ; i++ ) faceEdgeKeyValues[i].clear();
#endif // NEW_CODE

			if( _oldNCount<sliceData.nodeCount )
			{
				_oldNCount = sliceData.nodeCount;
				FreePointer( mcIndices );
				if( sliceData.nodeCount>0 ) mcIndices = AllocPointer< char >( _oldNCount );
			}
			if( _oldCCount<sliceData.cCount )
			{
				_oldCCount = sliceData.cCount;
				FreePointer( cornerValues ) ; FreePointer( cornerGradients ) ; FreePointer( cornerSet );
				if( sliceData.cCount>0 )
				{
					cornerValues = AllocPointer< Real >( _oldCCount );
					if( nonLinearFit ) cornerGradients = AllocPointer< Point< Real , Dim > >( _oldCCount );
					cornerSet = AllocPointer< char >( _oldCCount );
				}
			}
			if( _oldECount<sliceData.eCount )
			{
				_oldECount = sliceData.eCount;
				FreePointer( edgeKeys ) ; FreePointer( edgeSet );
#ifdef NEW_HASH
				edgeKeys = AllocPointer< _Key >( _oldECount );
#else // !NEW_HASH
				edgeKeys = AllocPointer< long long >( _oldECount );
#endif // NEW_HASH
				edgeSet = AllocPointer< char >( _oldECount );
			}
			if( _oldFCount<sliceData.fCount )
			{
				_oldFCount = sliceData.fCount;
				FreePointer( faceEdges ) ; FreePointer( faceSet );
				faceEdges = AllocPointer< _FaceEdges >( _oldFCount );
				faceSet = AllocPointer< char >( _oldFCount );
			}

			if( sliceData.cCount>0 ) memset( cornerSet , 0 , sizeof( char ) * sliceData.cCount );
			if( sliceData.eCount>0 ) memset(   edgeSet , 0 , sizeof( char ) * sliceData.eCount );
			if( sliceData.fCount>0 ) memset(   faceSet , 0 , sizeof( char ) * sliceData.fCount );
		}
	protected:
#ifdef NEW_CODE
		node_index_type _oldCCount , _oldECount , _oldFCount;
		node_index_type _oldNCount;
#else // !NEW_CODE
		int _oldCCount , _oldECount , _oldFCount , _oldNCount;
#endif // NEW_CODE
	};

	///////////////////
	// _XSliceValues //
	///////////////////
	struct _XSliceValues
	{
		typename SliceData::XSliceTableData xSliceData;
#ifdef NEW_HASH
		Pointer( _Key ) edgeKeys ; Pointer( char ) edgeSet;
#else // !NEW_HASH
		Pointer( long long ) edgeKeys ; Pointer( char ) edgeSet;
#endif // NEW_HASH
		Pointer( _FaceEdges ) faceEdges ; Pointer( char ) faceSet;
#ifdef NEW_HASH
		std::unordered_map< _Key , std::vector< _IsoEdge > , typename _Key::Hasher > faceEdgeMap;
#ifdef NEW_CODE
		std::unordered_map< _Key , std::pair< node_index_type , Vertex > , typename _Key::Hasher > edgeVertexMap;
#else // !NEW_CODE
		std::unordered_map< _Key , std::pair< int , Vertex > , typename _Key::Hasher > edgeVertexMap;
#endif // NEW_CODE
		std::unordered_map< _Key , _Key , typename _Key::Hasher > vertexPairMap;
#ifdef NEW_CODE
		std::vector< std::vector< std::pair< _Key , std::pair< node_index_type , Vertex > > > > edgeVertexKeyValues;
#else // !NEW_CODE
		std::vector< std::vector< std::pair< _Key , std::pair< int , Vertex > > > > edgeVertexKeyValues;
#endif // NEW_CODE
		std::vector< std::vector< std::pair< _Key , _Key > > > vertexPairKeyValues;
		std::vector< std::vector< std::pair< _Key , std::vector< _IsoEdge > > > > faceEdgeKeyValues;
#else // !NEW_HASH
		std::unordered_map< long long , std::vector< _IsoEdge > > faceEdgeMap;
#ifdef NEW_CODE
		std::unordered_map< long long , std::pair< node_index_type, Vertex > > edgeVertexMap;
#else // !NEW_CODE
		std::unordered_map< long long , std::pair< int, Vertex > > edgeVertexMap;
#endif // NEW_CODE
		std::unordered_map< long long , long long > vertexPairMap;
#ifdef NEW_CODE
		std::vector< std::vector< std::pair< long long , std::pair< node_index_type , Vertex > > > > edgeVertexKeyValues;
#else // !NEW_CODE
		std::vector< std::vector< std::pair< long long , std::pair< int , Vertex > > > > edgeVertexKeyValues;
#endif // NEW_CODE
		std::vector< std::vector< std::pair< long long , long long > > > vertexPairKeyValues;
		std::vector< std::vector< std::pair< long long , std::vector< _IsoEdge > > > > faceEdgeKeyValues;
#endif // NEW_HASH

#ifdef NEW_THREADS
		_XSliceValues( ThreadPool &tp )
#else // !NEW_THREADS
		_XSliceValues( void )
#endif // NEW_THREADS
		{
			_oldECount = _oldFCount = 0;
#ifdef NEW_HASH
			edgeKeys = NullPointer( _Key ) ; edgeSet = NullPointer( char );
#else // !NEW_HASH
			edgeKeys = NullPointer( long long ) ; edgeSet = NullPointer( char );
#endif // NEW_HASH
			faceEdges = NullPointer( _FaceEdges ) ; faceSet = NullPointer( char );
#ifdef NEW_THREADS
			edgeVertexKeyValues.resize( tp.threadNum() );
			vertexPairKeyValues.resize( tp.threadNum() );
			faceEdgeKeyValues.resize( tp.threadNum() );
#else // !NEW_THREADS
			edgeVertexKeyValues.resize( omp_get_max_threads() );
			vertexPairKeyValues.resize( omp_get_max_threads() );
			faceEdgeKeyValues.resize( omp_get_max_threads() );
#endif // NEW_THREADS
		}
		~_XSliceValues( void )
		{
			_oldECount = _oldFCount = 0;
			FreePointer( edgeKeys ) ; FreePointer( edgeSet );
			FreePointer( faceEdges ) ; FreePointer( faceSet );
		}
		void setEdgeVertexMap( void )
		{
#ifdef NEW_CODE
			for( node_index_type i=0 ; i<(node_index_type)edgeVertexKeyValues.size() ; i++ )
#else // !NEW_CODE
			for( int i=0 ; i<edgeVertexKeyValues.size() ; i++ )
#endif // NEW_CODE
			{
				for( int j=0 ; j<edgeVertexKeyValues[i].size() ; j++ ) edgeVertexMap[ edgeVertexKeyValues[i][j].first ] = edgeVertexKeyValues[i][j].second;
				edgeVertexKeyValues[i].clear();
			}
		}
		void setVertexPairMap( void )
		{
#ifdef NEW_CODE
			for( node_index_type i=0 ; i<(node_index_type)vertexPairKeyValues.size() ; i++ )
#else // !NEW_CODE
			for( int i=0 ; i<vertexPairKeyValues.size() ; i++ )
#endif // NEW_CODE
			{
				for( int j=0 ; j<vertexPairKeyValues[i].size() ; j++ )
				{
					vertexPairMap[ vertexPairKeyValues[i][j].first ] = vertexPairKeyValues[i][j].second;
					vertexPairMap[ vertexPairKeyValues[i][j].second ] = vertexPairKeyValues[i][j].first;
				}
				vertexPairKeyValues[i].clear();
			}
		}
		void setFaceEdgeMap( void )
		{
#ifdef NEW_CODE
			for( node_index_type i=0 ; i<(node_index_type)faceEdgeKeyValues.size() ; i++ )
#else // !NEW_CODE
			for( int i=0 ; i<faceEdgeKeyValues.size() ; i++ )
#endif // NEW_CODE
			{
				for( int j=0 ; j<faceEdgeKeyValues[i].size() ; j++ )
				{
					auto iter = faceEdgeMap.find( faceEdgeKeyValues[i][j].first );
					if( iter==faceEdgeMap.end() ) faceEdgeMap[ faceEdgeKeyValues[i][j].first ] = faceEdgeKeyValues[i][j].second;
					else for( int k=0 ; k<faceEdgeKeyValues[i][j].second.size() ; k++ ) iter->second.push_back( faceEdgeKeyValues[i][j].second[k] );
				}
				faceEdgeKeyValues[i].clear();
			}
		}
		void reset( void )
		{
			faceEdgeMap.clear() , edgeVertexMap.clear() , vertexPairMap.clear();
#ifdef NEW_CODE
			for( node_index_type i=0 ; i<(node_index_type)edgeVertexKeyValues.size() ; i++ ) edgeVertexKeyValues[i].clear();
			for( node_index_type i=0 ; i<(node_index_type)vertexPairKeyValues.size() ; i++ ) vertexPairKeyValues[i].clear();
			for( node_index_type i=0 ; i<(node_index_type)faceEdgeKeyValues.size() ; i++ ) faceEdgeKeyValues[i].clear();
#else // !NEW_CODE
			for( int i=0 ; i<edgeVertexKeyValues.size() ; i++ ) edgeVertexKeyValues[i].clear();
			for( int i=0 ; i<vertexPairKeyValues.size() ; i++ ) vertexPairKeyValues[i].clear();
			for( int i=0 ; i<faceEdgeKeyValues.size() ; i++ ) faceEdgeKeyValues[i].clear();
#endif // NEW_CODE

			if( _oldECount<xSliceData.eCount )
			{
				_oldECount = xSliceData.eCount;
				FreePointer( edgeKeys ) ; FreePointer( edgeSet );
#ifdef NEW_HASH
				edgeKeys = AllocPointer< _Key >( _oldECount );
#else // !NEW_HASH
				edgeKeys = AllocPointer< long long >( _oldECount );
#endif // NEW_HASH
				edgeSet = AllocPointer< char >( _oldECount );
			}
			if( _oldFCount<xSliceData.fCount )
			{
				_oldFCount = xSliceData.fCount;
				FreePointer( faceEdges ) ; FreePointer( faceSet );
				faceEdges = AllocPointer< _FaceEdges >( _oldFCount );
				faceSet = AllocPointer< char >( _oldFCount );
			}
			if( xSliceData.eCount>0 ) memset( edgeSet , 0 , sizeof( char ) * xSliceData.eCount );
			if( xSliceData.fCount>0 ) memset( faceSet , 0 , sizeof( char ) * xSliceData.fCount );
		}

	protected:
#ifdef NEW_CODE
		node_index_type _oldECount , _oldFCount;
#else // !NEW_CODE
		int _oldECount , _oldFCount;
#endif // NEW_CODE
	};

	/////////////////
	// _SlabValues //
	/////////////////
	struct _SlabValues
	{
	protected:
#ifdef NEW_THREADS
		std::vector< _XSliceValues > _xSliceValues;
		std::vector< _SliceValues > _sliceValues;
#else // !NEW_THREADS
		_XSliceValues _xSliceValues[2];
		_SliceValues _sliceValues[2];
#endif // NEW_THREADS
	public:
#ifdef NEW_THREADS
		_SlabValues( ThreadPool &tp )
		{
			_xSliceValues.reserve( 2 );
			_sliceValues.reserve( 2 );
			for( unsigned int i=0 ; i<2 ; i++ ) _xSliceValues.emplace_back( tp ) , _sliceValues.emplace_back( tp );
		}
#endif // NEW_THREADS
		_SliceValues& sliceValues( int idx ){ return _sliceValues[idx&1]; }
		const _SliceValues& sliceValues( int idx ) const { return _sliceValues[idx&1]; }
		_XSliceValues& xSliceValues( int idx ){ return _xSliceValues[idx&1]; }
		const _XSliceValues& xSliceValues( int idx ) const { return _xSliceValues[idx&1]; }
	};

#ifdef NEW_THREADS
	template< unsigned int ... FEMSigs >
	static void _SetSliceIsoCorners( ThreadPool &tp , const FEMTree< Dim , Real >& tree , ConstPointer( Real ) coefficients , ConstPointer( Real ) coarseCoefficients , Real isoValue , LocalDepth depth , int slice ,         std::vector< _SlabValues >& slabValues , const _Evaluator< UIntPack< FEMSigs ... > , 1 >& evaluator )
	{
		if( slice>0          ) _SetSliceIsoCorners< FEMSigs ... >( tp , tree , coefficients , coarseCoefficients , isoValue , depth , slice , HyperCube::FRONT , slabValues , evaluator );
		if( slice<(1<<depth) ) _SetSliceIsoCorners< FEMSigs ... >( tp , tree , coefficients , coarseCoefficients , isoValue , depth , slice , HyperCube::BACK  , slabValues , evaluator );
	}
#else // !NEW_THREADS
	template< unsigned int ... FEMSigs >
	static void _SetSliceIsoCorners( const FEMTree< Dim , Real >& tree , ConstPointer( Real ) coefficients , ConstPointer( Real ) coarseCoefficients , Real isoValue , LocalDepth depth , int slice ,         std::vector< _SlabValues >& slabValues , const _Evaluator< UIntPack< FEMSigs ... > , 1 >& evaluator )
	{
		if( slice>0          ) _SetSliceIsoCorners< FEMSigs ... >( tree , coefficients , coarseCoefficients , isoValue , depth , slice , HyperCube::FRONT , slabValues , evaluator );
		if( slice<(1<<depth) ) _SetSliceIsoCorners< FEMSigs ... >( tree , coefficients , coarseCoefficients , isoValue , depth , slice , HyperCube::BACK  , slabValues , evaluator );
	}
#endif // NEW_THREADS
	template< unsigned int ... FEMSigs >
#ifdef NEW_THREADS
	static void _SetSliceIsoCorners( ThreadPool &tp , const FEMTree< Dim , Real >& tree , ConstPointer( Real ) coefficients , ConstPointer( Real ) coarseCoefficients , Real isoValue , LocalDepth depth , int slice , HyperCube::Direction zDir , std::vector< _SlabValues >& slabValues , const _Evaluator< UIntPack< FEMSigs ... > , 1 >& evaluator )
#else // !NEW_THREADS
	static void _SetSliceIsoCorners( const FEMTree< Dim , Real >& tree , ConstPointer( Real ) coefficients , ConstPointer( Real ) coarseCoefficients , Real isoValue , LocalDepth depth , int slice , HyperCube::Direction zDir , std::vector< _SlabValues >& slabValues , const _Evaluator< UIntPack< FEMSigs ... > , 1 >& evaluator )
#endif // NEW_THREADS
	{
		static const unsigned int FEMDegrees[] = { FEMSignature< FEMSigs >::Degree ... };
		_SliceValues& sValues = slabValues[depth].sliceValues( slice );
		bool useBoundaryEvaluation = false;
		for( int d=0 ; d<Dim ; d++ ) if( FEMDegrees[d]==0 || ( FEMDegrees[d]==1 && sValues.cornerGradients ) ) useBoundaryEvaluation = true;
		std::vector< ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > > > neighborKeys( omp_get_max_threads() );
		std::vector< ConstCornerSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > > > bNeighborKeys( omp_get_max_threads() );
		if( useBoundaryEvaluation ) for( size_t i=0 ; i<neighborKeys.size() ; i++ ) bNeighborKeys[i].set( tree._localToGlobal( depth ) );
		else                        for( size_t i=0 ; i<neighborKeys.size() ; i++ )  neighborKeys[i].set( tree._localToGlobal( depth ) );
#ifdef NEW_THREADS
		tp.parallel_for( tree._sNodesBegin(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) , tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=tree._sNodesBegin(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i<tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i++ )
#else // !NEW_CODE
		for( int i=tree._sNodesBegin(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i<tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i] ) )
			{
				Real squareValues[ HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ];
#ifdef NEW_THREADS
				ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey = neighborKeys[ thread ];
				ConstCornerSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bNeighborKey = bNeighborKeys[ thread ];
#else // !NEW_THREADS
				ConstPointSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& neighborKey = neighborKeys[ omp_get_thread_num() ];
				ConstCornerSupportKey< UIntPack< FEMSignature< FEMSigs >::Degree ... > >& bNeighborKey = bNeighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
				TreeNode* leaf = tree._sNodes.treeNodes[i];
				if( !IsActiveNode< Dim >( leaf->children ) )
				{
					const typename SliceData::SquareCornerIndices& cIndices = sValues.sliceData.cornerIndices( leaf );

					bool isInterior = tree._isInteriorlySupported( UIntPack< FEMSignature< FEMSigs >::Degree ... >() , leaf->parent );
					if( useBoundaryEvaluation ) bNeighborKey.getNeighbors( leaf );
					else                         neighborKey.getNeighbors( leaf );

					for( typename HyperCube::Cube< Dim-1 >::template Element< 0 > _c ; _c<HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; _c++ )
					{
						typename HyperCube::Cube< Dim >::template Element< 0 > c( zDir , _c.index );
#ifdef NEW_CODE
						node_index_type vIndex = cIndices[_c.index];
#else // !NEW_CODE
						int vIndex = cIndices[_c.index];
#endif // NEW_CODE
						if( !sValues.cornerSet[vIndex] )
						{
							if( sValues.cornerGradients )
							{
								CumulativeDerivativeValues< Real , Dim , 1 > p;
								if( useBoundaryEvaluation ) p = tree.template _getCornerValues< Real , 1 >( bNeighborKey , leaf , c.index , coefficients , coarseCoefficients , evaluator , tree._maxDepth , isInterior );
								else                        p = tree.template _getCornerValues< Real , 1 >(  neighborKey , leaf , c.index , coefficients , coarseCoefficients , evaluator , tree._maxDepth , isInterior );
								sValues.cornerValues[vIndex] = p[0] , sValues.cornerGradients[vIndex] = Point< Real , Dim >( p[1] , p[2] , p[3] );
							}
							else
							{
								if( useBoundaryEvaluation ) sValues.cornerValues[vIndex] = tree.template _getCornerValues< Real , 0 >( bNeighborKey , leaf , c.index , coefficients , coarseCoefficients , evaluator , tree._maxDepth , isInterior )[0];
								else                        sValues.cornerValues[vIndex] = tree.template _getCornerValues< Real , 0 >(  neighborKey , leaf , c.index , coefficients , coarseCoefficients , evaluator , tree._maxDepth , isInterior )[0];
							}
							sValues.cornerSet[vIndex] = 1;
						}
						squareValues[_c.index] = sValues.cornerValues[ vIndex ];
						TreeNode* node = leaf;
						LocalDepth _depth = depth;
						int _slice = slice;
						while( tree._isValidSpaceNode( node->parent ) && (node-node->parent->children)==c.index )
						{
							node = node->parent , _depth-- , _slice >>= 1;
							_SliceValues& _sValues = slabValues[_depth].sliceValues( _slice );
							const typename SliceData::SquareCornerIndices& _cIndices = _sValues.sliceData.cornerIndices( node );
#ifdef NEW_CODE
							node_index_type _vIndex = _cIndices[_c.index];
#else // !NEW_CODE
							int _vIndex = _cIndices[_c.index];
#endif // NEW_CODE
							_sValues.cornerValues[_vIndex] = sValues.cornerValues[vIndex];
							if( _sValues.cornerGradients ) _sValues.cornerGradients[_vIndex] = sValues.cornerGradients[vIndex];
							_sValues.cornerSet[_vIndex] = 1;
						}
					}
					sValues.mcIndices[ i - sValues.sliceData.nodeOffset ] = HyperCube::Cube< Dim-1 >::MCIndex( squareValues , isoValue );
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
	/////////////////
	// _VertexData //
	/////////////////
	class _VertexData
	{
	public:
#ifdef NEW_HASH
		static _Key EdgeIndex( const TreeNode* node , typename HyperCube::Cube< Dim >::template Element< 1 > e , int maxDepth )
		{
			_Key key;
			const HyperCube::Direction* x = SliceData::template HyperCubeTables< Dim , 1 >::Directions[ e.index ];
			int d , off[Dim];
			node->depthAndOffset( d , off );
			for( int dd=0 ; dd<Dim ; dd++ )
			{
				if( x[dd]==HyperCube::CROSS )
				{
					key[(dd+0)%3] = (int)BinaryNode::CornerIndex( maxDepth+1 , d+1 , off[(dd+0)%3]<<1 , 1 );
					key[(dd+1)%3] = (int)BinaryNode::CornerIndex( maxDepth+1 , d   , off[(dd+1)%3] , x[(dd+1)%3]==HyperCube::BACK ? 0 : 1 );
					key[(dd+2)%3] = (int)BinaryNode::CornerIndex( maxDepth+1 , d   , off[(dd+2)%3] , x[(dd+2)%3]==HyperCube::BACK ? 0 : 1 );
				}
			}
			return key;
		}

		static _Key FaceIndex( const TreeNode* node , typename HyperCube::Cube< Dim >::template Element< Dim-1 > f , int maxDepth )
		{
			_Key key;
			const HyperCube::Direction* x = SliceData::template HyperCubeTables< Dim , 2 >::Directions[ f.index ];
			int d , o[Dim];
			node->depthAndOffset( d , o );
			for( int dd=0 ; dd<Dim ; dd++ )
				if( x[dd]==HyperCube::CROSS ) key[dd] = (int)BinaryNode::CornerIndex( maxDepth+1 , d+1 , o[dd]<<1 , 1 );
				else                          key[dd] = (int)BinaryNode::CornerIndex( maxDepth+1 , d   , o[dd]    , x[dd]==HyperCube::BACK ? 0 : 1 );
			return key;
		}
#else // !NEW_HASH
#pragma message( "[WARNING] Replace me with a smarter hashing function" )
		static const int VERTEX_COORDINATE_SHIFT = ( sizeof( long long ) * 8 ) / Dim;
		static long long Index( const int index[Dim] ){ long long idx=0 ; for( int dd=0 ; dd<Dim ; dd++ ) idx |= ( ( long long )index[dd] )<<(dd*VERTEX_COORDINATE_SHIFT) ; return idx; }

		static long long EdgeIndex( const TreeNode* node , typename HyperCube::Cube< Dim >::template Element< 1 > e , int maxDepth , int idx[Dim] )
		{
			const HyperCube::Direction* x = SliceData::template HyperCubeTables< Dim , 1 >::Directions[ e.index ];
			int d , off[Dim];
			node->depthAndOffset( d , off );
			for( int dd=0 ; dd<Dim ; dd++ )
			{
				if( x[dd]==HyperCube::CROSS )
				{
					idx[(dd+0)%3] = BinaryNode::CornerIndex( maxDepth+1 , d+1 , off[(dd+0)%3]<<1 , 1 );
					idx[(dd+1)%3] = BinaryNode::CornerIndex( maxDepth+1 , d   , off[(dd+1)%3] , x[(dd+1)%3]==HyperCube::BACK ? 0 : 1 );
					idx[(dd+2)%3] = BinaryNode::CornerIndex( maxDepth+1 , d   , off[(dd+2)%3] , x[(dd+2)%3]==HyperCube::BACK ? 0 : 1 );
				}
			}
			return Index( idx );
		}
		static long long EdgeIndex( const TreeNode* node , typename HyperCube::Cube< Dim >::template Element< 1 > e , int maxDepth ){ int idx[Dim] ; return EdgeIndex( node , e , maxDepth , idx ); }

		static long long FaceIndex( const TreeNode* node , typename HyperCube::Cube< Dim >::template Element< Dim-1 > f , int maxDepth , int idx[Dim] )
		{
			const HyperCube::Direction* x = SliceData::template HyperCubeTables< Dim , 2 >::Directions[ f.index ];
			int d , o[Dim];
			node->depthAndOffset( d , o );
			for( int dd=0 ; dd<Dim ; dd++ )
				if( x[dd]==HyperCube::CROSS ) idx[dd] = BinaryNode::CornerIndex( maxDepth+1 , d+1 , o[dd]<<1 , 1 );
				else                          idx[dd] = BinaryNode::CornerIndex( maxDepth+1 , d   , o[dd]    , x[dd]==HyperCube::BACK ? 0 : 1 );
			return Index( idx );
		}
		static long long FaceIndex( const TreeNode* node , typename HyperCube::Cube< Dim >::template Element< Dim-1 > f , int maxDepth ){ int idx[Dim] ; return FaceIndex( node , f , maxDepth , idx ); }
#endif // NEW_HASH
	};

#ifdef NEW_THREADS
	template< unsigned int WeightDegree , typename Data , unsigned int DataSig >
	static void _SetSliceIsoVertices( ThreadPool &tp , const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , LocalDepth depth , int slice , node_index_type& vOffset , CoredMeshData< Vertex , node_index_type >& mesh , std::vector< _SlabValues >& slabValues , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
	{
		if( slice>0          ) _SetSliceIsoVertices< WeightDegree , Data , DataSig >( tp , tree , pointEvaluator , densityWeights , data , isoValue , depth , slice , HyperCube::FRONT , vOffset , mesh , slabValues , SetVertex );
		if( slice<(1<<depth) ) _SetSliceIsoVertices< WeightDegree , Data , DataSig >( tp , tree , pointEvaluator , densityWeights , data , isoValue , depth , slice , HyperCube::BACK  , vOffset , mesh , slabValues , SetVertex );
	}
#else // !NEW_THREADS
	template< unsigned int WeightDegree , typename Data , unsigned int DataSig >
#ifdef NEW_CODE
	static void _SetSliceIsoVertices( const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , LocalDepth depth , int slice , node_index_type& vOffset , CoredMeshData< Vertex , node_index_type >& mesh , std::vector< _SlabValues >& slabValues , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
#else // !NEW_CODE
	static void _SetSliceIsoVertices( const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , LocalDepth depth , int slice , int& vOffset , CoredMeshData< Vertex >& mesh , std::vector< _SlabValues >& slabValues , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
#endif // NEW_CODE
	{
		if( slice>0          ) _SetSliceIsoVertices< WeightDegree , Data , DataSig >( tree , pointEvaluator , densityWeights , data , isoValue , depth , slice , HyperCube::FRONT , vOffset , mesh , slabValues , SetVertex );
		if( slice<(1<<depth) ) _SetSliceIsoVertices< WeightDegree , Data , DataSig >( tree , pointEvaluator , densityWeights , data , isoValue , depth , slice , HyperCube::BACK  , vOffset , mesh , slabValues , SetVertex );
	}
#endif // NEW_THREADS
	template< unsigned int WeightDegree , typename Data , unsigned int DataSig >
#ifdef NEW_CODE
#ifdef NEW_THREADS
	static void _SetSliceIsoVertices( ThreadPool &tp , const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , LocalDepth depth , int slice , HyperCube::Direction zDir , node_index_type& vOffset , CoredMeshData< Vertex , node_index_type >& mesh , std::vector< _SlabValues >& slabValues , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
#else // !NEW_THREADS
	static void _SetSliceIsoVertices( const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , LocalDepth depth , int slice , HyperCube::Direction zDir , node_index_type& vOffset , CoredMeshData< Vertex , node_index_type >& mesh , std::vector< _SlabValues >& slabValues , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
#endif // NEW_THREADS
#else // !NEW_CODE
	static void _SetSliceIsoVertices( const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , LocalDepth depth , int slice , HyperCube::Direction zDir , int& vOffset , CoredMeshData< Vertex >& mesh , std::vector< _SlabValues >& slabValues , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
#endif // NEW_CODE
	{
		static const unsigned int DataDegree = FEMSignature< DataSig >::Degree;
		_SliceValues& sValues = slabValues[depth].sliceValues( slice );
		// [WARNING] In the case Degree=2, these two keys are the same, so we don't have to maintain them separately.
#ifdef NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
		std::vector< ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > > > weightKeys( tp.threadNum() );
		std::vector< ConstPointSupportKey< IsotropicUIntPack< Dim , DataDegree > > > dataKeys( tp.threadNum() );
#else // !NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
		std::vector< ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > > > weightKeys( omp_get_max_threads() );
		std::vector< ConstPointSupportKey< IsotropicUIntPack< Dim , DataDegree > > > dataKeys( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( tree._localToGlobal( depth ) ) , weightKeys[i].set( tree._localToGlobal( depth ) ) , dataKeys[i].set( tree._localToGlobal( depth ) );
#ifdef NEW_THREADS
		tp.parallel_for( tree._sNodesBegin(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) , tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=tree._sNodesBegin(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i<tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i++ )
#else // !NEW_CODE
		for( int i=tree._sNodesBegin(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i<tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i] ) )
			{
#ifdef NEW_THREADS
				ConstOneRingNeighborKey& neighborKey =  neighborKeys[ thread ];
				ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey = weightKeys[ thread ];
				ConstPointSupportKey< IsotropicUIntPack< Dim , DataDegree > >& dataKey = dataKeys[ thread ];
#else // !NEW_THREADS
				ConstOneRingNeighborKey& neighborKey =  neighborKeys[ omp_get_thread_num() ];
				ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey = weightKeys[ omp_get_thread_num() ];
				ConstPointSupportKey< IsotropicUIntPack< Dim , DataDegree > >& dataKey = dataKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
				TreeNode* leaf = tree._sNodes.treeNodes[i];
				if( !IsActiveNode< Dim >( leaf->children ) )
				{
#ifdef NEW_CODE
					node_index_type idx = ( i - sValues.sliceData.nodeOffset );
#else // !NEW_CODE
					int idx = i - sValues.sliceData.nodeOffset;
#endif // NEW_CODE
					const typename SliceData::SquareEdgeIndices& eIndices = sValues.sliceData.edgeIndices( leaf );
					if( HyperCube::Cube< Dim-1 >::HasMCRoots( sValues.mcIndices[idx] ) )
					{
						neighborKey.getNeighbors( leaf );
						if( densityWeights ) weightKey.getNeighbors( leaf );
						if( data ) dataKey.getNeighbors( leaf );

						for( typename HyperCube::Cube< Dim-1 >::template Element< 1 > _e ; _e<HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; _e++ )
							if( HyperCube::Cube< 1 >::HasMCRoots( HyperCube::Cube< Dim-1 >::ElementMCIndex( _e , sValues.mcIndices[idx] ) ) )
							{
								typename HyperCube::Cube< Dim >::template Element< 1 > e( zDir , _e.index );
#ifdef NEW_CODE
								node_index_type vIndex = eIndices[_e.index];
#else // !NEW_CODE
								int vIndex = eIndices[_e.index];
#endif // NEW_CODE
								if( !sValues.edgeSet[vIndex] )
								{
									Vertex vertex;
#ifdef NEW_HASH
									_Key key = _VertexData::EdgeIndex( leaf , e , tree._localToGlobal( tree._maxDepth ) );
#else // !NEW_HASH
									long long key = _VertexData::EdgeIndex( leaf , e , tree._localToGlobal( tree._maxDepth ) );
#endif // NEW_HASH
									_GetIsoVertex< WeightDegree , Data , DataSig >( tree , pointEvaluator , densityWeights , data , isoValue , weightKey , dataKey , leaf , _e , zDir , sValues , vertex , SetVertex );
									bool stillOwner = false;
#ifdef NEW_CODE
									std::pair< node_index_type , Vertex > hashed_vertex;
#else // !NEW_CODE
									std::pair< int , Vertex > hashed_vertex;
#endif // NEW_CODE
#ifdef NEW_THREADS
									{
										std::lock_guard< std::mutex > lock( _pointInsertionMutex );
										if( !sValues.edgeSet[vIndex] )
										{
											mesh.addOutOfCorePoint( vertex );
											sValues.edgeSet[ vIndex ] = 1;
#ifdef NEW_CODE
											hashed_vertex = std::pair< node_index_type , Vertex >( vOffset , vertex );
#else // !NEW_CODE
											hashed_vertex = std::pair< int , Vertex >( vOffset , vertex );
#endif // NEW_CODE
											sValues.edgeKeys[ vIndex ] = key;
											vOffset++;
											stillOwner = true;
										}
									}
#else // !NEW_THREADS
#pragma omp critical (add_point_access)
									if( !sValues.edgeSet[vIndex] )
									{
										mesh.addOutOfCorePoint( vertex );
										sValues.edgeSet[ vIndex ] = 1;
#ifdef NEW_CODE
										hashed_vertex = std::pair< node_index_type , Vertex >( vOffset , vertex );
#else // !NEW_CODE
										hashed_vertex = std::pair< int , Vertex >( vOffset , vertex );
#endif // NEW_CODE
										sValues.edgeKeys[ vIndex ] = key;
										vOffset++;
										stillOwner = true;
									}
#endif // NEW_THREADS
#ifdef NEW_CODE
#ifdef NEW_THREADS
#ifdef NEW_HASH
									if( stillOwner ) sValues.edgeVertexKeyValues[ thread ].push_back( std::pair< _Key , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#else // !NEW_HASH
									if( stillOwner ) sValues.edgeVertexKeyValues[ thread ].push_back( std::pair< long long , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#endif // NEW_HASH
#else // !NEW_THREADS
#ifdef NEW_HASH
									if( stillOwner ) sValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< _Key , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#else // !NEW_HASH
									if( stillOwner ) sValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#endif // NEW_HASH
#endif // NEW_THREADS
#else // !NEW_CODE
									if( stillOwner ) sValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::pair< int , Vertex > >( key , hashed_vertex ) );
#endif // NEW_CODE
									if( stillOwner )
									{
										// We only need to pass the iso-vertex down if the edge it lies on is adjacent to a coarser leaf
										auto IsNeeded = [&]( unsigned int depth )
										{
											bool isNeeded = false;
											typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > my_ic = SliceData::template HyperCubeTables< Dim , 1 >::IncidentCube[e.index];
											for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 1 >() ; ic++ ) if( ic!=my_ic )
											{
												unsigned int xx = SliceData::template HyperCubeTables< Dim , 1 >::CellOffset[e.index][ic.index];
												isNeeded |= !tree._isValidSpaceNode( neighborKey.neighbors[ tree._localToGlobal( depth ) ].neighbors.data[xx] );
											}
											return isNeeded;
										};
										if( IsNeeded( depth ) )
										{
											const typename HyperCube::Cube< Dim >::template Element< Dim-1 > *f = SliceData::template HyperCubeTables< Dim , 1 , Dim-1 >::OverlapElements[e.index];
											for( int k=0 ; k<2 ; k++ )
											{
												TreeNode* node = leaf;
												LocalDepth _depth = depth;
												int _slice = slice;
												while( tree._isValidSpaceNode( node->parent ) && SliceData::template HyperCubeTables< Dim , 2 , 0 >::Overlap[f[k].index][(unsigned int)(node-node->parent->children) ] )
												{
													node = node->parent , _depth-- , _slice >>= 1;
													_SliceValues& _sValues = slabValues[_depth].sliceValues( _slice );
#ifdef NEW_CODE
#ifdef NEW_THREADS
#ifdef NEW_HASH
													_sValues.edgeVertexKeyValues[ thread ].push_back( std::pair< _Key , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#else // !NEW_HASH
													_sValues.edgeVertexKeyValues[ thread ].push_back( std::pair< long long , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#endif // NEW_HASH
#else // !NEW_THREADS
#ifdef NEW_HASH
													_sValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< _Key , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#else // !NEW_HASH
													_sValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#endif // NEW_HASH
#endif // NEW_THREADS
#else // !NEW_CODE
													_sValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::pair< int , Vertex > >( key , hashed_vertex ) );
#endif // NEW_CODE
													if( !IsNeeded( _depth ) ) break;
												}
											}
										}
									}
								}
							}
					}
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}

	////////////////////
	// Iso-Extraction //
	////////////////////
	template< unsigned int WeightDegree , typename Data , unsigned int DataSig >
#ifdef NEW_CODE
#ifdef NEW_THREADS
	static void _SetXSliceIsoVertices( ThreadPool &tp , const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , LocalDepth depth , int slab , node_index_type &vOffset , CoredMeshData< Vertex , node_index_type >& mesh , std::vector< _SlabValues >& slabValues , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
#else // !NEW_THREADS
	static void _SetXSliceIsoVertices( const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , LocalDepth depth , int slab , node_index_type &vOffset , CoredMeshData< Vertex , node_index_type >& mesh , std::vector< _SlabValues >& slabValues , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
#endif // NEW_THREADS
#else // !NEW_CODE
	static void _SetXSliceIsoVertices( const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , LocalDepth depth , int slab , int& vOffset , CoredMeshData< Vertex >& mesh , std::vector< _SlabValues >& slabValues , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
#endif // NEW_CODE
	{
		static const unsigned int DataDegree = FEMSignature< DataSig >::Degree;
		_SliceValues& bValues = slabValues[depth].sliceValues ( slab   );
		_SliceValues& fValues = slabValues[depth].sliceValues ( slab+1 );
		_XSliceValues& xValues = slabValues[depth].xSliceValues( slab   );

		// [WARNING] In the case Degree=2, these two keys are the same, so we don't have to maintain them separately.
#ifdef NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
		std::vector< ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > > > weightKeys( tp.threadNum() );
		std::vector< ConstPointSupportKey< IsotropicUIntPack< Dim , DataDegree > > > dataKeys( tp.threadNum() );
#else // !NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
		std::vector< ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > > > weightKeys( omp_get_max_threads() );
		std::vector< ConstPointSupportKey< IsotropicUIntPack< Dim , DataDegree > > > dataKeys( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( tree._localToGlobal( depth ) ) , weightKeys[i].set( tree._localToGlobal( depth ) ) , dataKeys[i].set( tree._localToGlobal( depth ) );
#ifdef NEW_THREADS
		tp.parallel_for( tree._sNodesBegin(depth,slab) , tree._sNodesEnd(depth,slab) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=tree._sNodesBegin(depth,slab) ; i<tree._sNodesEnd(depth,slab) ; i++ )
#else // !NEW_CODE
		for( int i=tree._sNodesBegin(depth,slab) ; i<tree._sNodesEnd(depth,slab) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i] ) )
			{
#ifdef NEW_THREADS
				ConstOneRingNeighborKey& neighborKey =  neighborKeys[ thread ];
				ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey = weightKeys[ thread ];
				ConstPointSupportKey< IsotropicUIntPack< Dim , DataDegree > >& dataKey = dataKeys[ thread ];
#else // !NEW_THREADS
				ConstOneRingNeighborKey& neighborKey =  neighborKeys[ omp_get_thread_num() ];
				ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey = weightKeys[ omp_get_thread_num() ];
				ConstPointSupportKey< IsotropicUIntPack< Dim , DataDegree > >& dataKey = dataKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
				TreeNode* leaf = tree._sNodes.treeNodes[i];
				if( !IsActiveNode< Dim >( leaf->children ) )
				{
					unsigned char mcIndex = ( bValues.mcIndices[ i - bValues.sliceData.nodeOffset ] ) | ( fValues.mcIndices[ i - fValues.sliceData.nodeOffset ] )<<4;
					const typename SliceData::SquareCornerIndices& eIndices = xValues.xSliceData.edgeIndices( leaf );
					if( HyperCube::Cube< Dim >::HasMCRoots( mcIndex ) )
					{
						neighborKey.getNeighbors( leaf );
						if( densityWeights ) weightKey.getNeighbors( leaf );
						if( data ) dataKey.getNeighbors( leaf );
						for( typename HyperCube::Cube< Dim-1 >::template Element< 0 > _c ; _c<HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; _c++ )
						{
							typename HyperCube::Cube< Dim >::template Element< 1 > e( HyperCube::CROSS , _c.index );
							unsigned int _mcIndex = HyperCube::Cube< Dim >::ElementMCIndex( e , mcIndex );
							if( HyperCube::Cube< 1 >::HasMCRoots( _mcIndex ) )
							{
#ifdef NEW_CODE
								node_index_type vIndex = eIndices[_c.index];
#else // !NEW_CODE
								int vIndex = eIndices[_c.index];
#endif // NEW_CODE
								if( !xValues.edgeSet[vIndex] )
								{
									Vertex vertex;
#ifdef NEW_HASH
									_Key key = _VertexData::EdgeIndex( leaf , e.index , tree._localToGlobal( tree._maxDepth ) );
#else // !NEW_HASH
									long long key = _VertexData::EdgeIndex( leaf , e.index , tree._localToGlobal( tree._maxDepth ) );
#endif // NEW_HASH
									_GetIsoVertex< WeightDegree , Data , DataSig >( tree , pointEvaluator , densityWeights , data , isoValue , weightKey , dataKey , leaf , _c , bValues , fValues , vertex , SetVertex );
									bool stillOwner = false;
#ifdef NEW_CODE
									std::pair< node_index_type , Vertex > hashed_vertex;
#else // !NEW_CODE
									std::pair< int , Vertex > hashed_vertex;
#endif // NEW_CODE
#ifdef NEW_THREADS
									{
										std::lock_guard< std::mutex > lock( _pointInsertionMutex );
										if( !xValues.edgeSet[vIndex] )
										{
											mesh.addOutOfCorePoint( vertex );
											xValues.edgeSet[ vIndex ] = 1;
#ifdef NEW_CODE
											hashed_vertex = std::pair< node_index_type , Vertex >( vOffset , vertex );
#else // !NEW_CODE
											hashed_vertex = std::pair< int , Vertex >( vOffset , vertex );
#endif // NEW_CODE
											xValues.edgeKeys[ vIndex ] = key;
											vOffset++;
											stillOwner = true;
										}
									}
#else // !NEW_THREADS
#pragma omp critical (add_point_access)
									if( !xValues.edgeSet[vIndex] )
									{
										mesh.addOutOfCorePoint( vertex );
										xValues.edgeSet[ vIndex ] = 1;
#ifdef NEW_CODE
										hashed_vertex = std::pair< node_index_type , Vertex >( vOffset , vertex );
#else // !NEW_CODE
										hashed_vertex = std::pair< int , Vertex >( vOffset , vertex );
#endif // NEW_CODE
										xValues.edgeKeys[ vIndex ] = key;
										vOffset++;
										stillOwner = true;
									}
#endif // NEW_THREADS
#ifdef NEW_CODE
#ifdef NEW_THREADS
#ifdef NEW_HASH
									if( stillOwner ) xValues.edgeVertexKeyValues[ thread ].push_back( std::pair< _Key , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#else // !NEW_HASH
									if( stillOwner ) xValues.edgeVertexKeyValues[ thread ].push_back( std::pair< long long , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#endif // NEW_HASH
#else // !NEW_THREADS
#ifdef NEW_HASH
									if( stillOwner ) xValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< _Key , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#else // !NEW_HASH
									if( stillOwner ) xValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#endif // NEW_HASH
#endif // NEW_THREADS
#else // !NEW_CODE
									if( stillOwner ) xValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::pair< int , Vertex > >( key , hashed_vertex ) );
#endif // NEW_CODE
									if( stillOwner )
									{
										// We only need to pass the iso-vertex down if the edge it lies on is adjacent to a coarser leaf
										auto IsNeeded = [&]( unsigned int depth )
										{
											bool isNeeded = false;
											typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > my_ic = SliceData::template HyperCubeTables< Dim , 1 >::IncidentCube[e.index];
											for( typename HyperCube::Cube< Dim >::template IncidentCubeIndex< 1 > ic ; ic<HyperCube::Cube< Dim >::template IncidentCubeNum< 1 >() ; ic++ ) if( ic!=my_ic )
											{
												unsigned int xx = SliceData::template HyperCubeTables< Dim , 1 >::CellOffset[e.index][ic.index];
												isNeeded |= !tree._isValidSpaceNode( neighborKey.neighbors[ tree._localToGlobal( depth ) ].neighbors.data[xx] );
											}
											return isNeeded;
										};
										if( IsNeeded( depth ) )
										{
											const typename HyperCube::Cube< Dim >::template Element< Dim-1 > *f = SliceData::template HyperCubeTables< Dim , 1 , Dim-1 >::OverlapElements[e.index];
											for( int k=0 ; k<2 ; k++ )
											{
												TreeNode* node = leaf;
												LocalDepth _depth = depth;
												int _slab = slab;
												while( tree._isValidSpaceNode( node->parent ) && SliceData::template HyperCubeTables< Dim , 2 , 0 >::Overlap[f[k].index][(unsigned int)(node-node->parent->children) ] )
												{
													node = node->parent , _depth-- , _slab >>= 1;
													_XSliceValues& _xValues = slabValues[_depth].xSliceValues( _slab );
#ifdef NEW_CODE
#ifdef NEW_THREADS
#ifdef NEW_HASH
													_xValues.edgeVertexKeyValues[ thread ].push_back( std::pair< _Key , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#else // !NEW_HASH
													_xValues.edgeVertexKeyValues[ thread ].push_back( std::pair< long long , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#endif // NEW_HASH
#else // !NEW_THREADS
#ifdef NEW_HASH
													_xValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< _Key , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#else // !NEW_HASH
													_xValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::pair< node_index_type , Vertex > >( key , hashed_vertex ) );
#endif // NEW_HASH
#endif // NEW_THREADS
#else // !NEW_CODE
													_xValues.edgeVertexKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::pair< int , Vertex > >( key , hashed_vertex ) );
#endif // NEW_CODE
													if( !IsNeeded( _depth ) ) break;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
#ifdef NEW_THREADS
	static void _CopyFinerSliceIsoEdgeKeys( ThreadPool &tp , const FEMTree< Dim , Real >& tree , LocalDepth depth , int slice , std::vector< _SlabValues >& slabValues )
	{
		if( slice>0          ) _CopyFinerSliceIsoEdgeKeys( tp , tree , depth , slice , HyperCube::FRONT , slabValues );
		if( slice<(1<<depth) ) _CopyFinerSliceIsoEdgeKeys( tp , tree , depth , slice , HyperCube::BACK  , slabValues );
	}
#else // !NEW_THREADS
	static void _CopyFinerSliceIsoEdgeKeys( const FEMTree< Dim , Real >& tree , LocalDepth depth , int slice , std::vector< _SlabValues >& slabValues )
	{
		if( slice>0          ) _CopyFinerSliceIsoEdgeKeys( tree , depth , slice , HyperCube::FRONT , slabValues );
		if( slice<(1<<depth) ) _CopyFinerSliceIsoEdgeKeys( tree , depth , slice , HyperCube::BACK  , slabValues );
	}
#endif // NEW_THREADS
#ifdef NEW_THREADS
	static void _CopyFinerSliceIsoEdgeKeys( ThreadPool &tp , const FEMTree< Dim , Real >& tree , LocalDepth depth , int slice , HyperCube::Direction zDir , std::vector< _SlabValues >& slabValues )
#else // !NEW_THREADS
	static void _CopyFinerSliceIsoEdgeKeys( const FEMTree< Dim , Real >& tree , LocalDepth depth , int slice , HyperCube::Direction zDir , std::vector< _SlabValues >& slabValues )
#endif // NEW_THREADS
	{
		_SliceValues& pSliceValues = slabValues[depth  ].sliceValues(slice   );
		_SliceValues& cSliceValues = slabValues[depth+1].sliceValues(slice<<1);
		typename SliceData::SliceTableData& pSliceData = pSliceValues.sliceData;
		typename SliceData::SliceTableData& cSliceData = cSliceValues.sliceData;
#ifdef NEW_THREADS
		tp.parallel_for( tree._sNodesBegin(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) , tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=tree._sNodesBegin(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i<tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i++ )
#else // !NEW_CODE
		for( int i=tree._sNodesBegin(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i<tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i] ) ) if( IsActiveNode< Dim >( tree._sNodes.treeNodes[i]->children ) )
			{
#ifdef NEW_THREADS
#else // !NEW_THREADS
				int thread = omp_get_thread_num();
#endif // NEW_THREADS
				typename SliceData::SquareEdgeIndices& pIndices = pSliceData.edgeIndices( i );
				// Copy the edges that overlap the coarser edges
				for( typename HyperCube::Cube< Dim-1 >::template Element< 1 > _e ; _e<HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; _e++ )
				{
#ifdef NEW_CODE
					node_index_type pIndex = pIndices[_e.index];
#else // !NEW_CODE
					int pIndex = pIndices[_e.index];
#endif // NEW_CODE
					if( !pSliceValues.edgeSet[ pIndex ] )
					{
						typename HyperCube::Cube< Dim >::template Element< 1 > e( zDir , _e.index );
						const typename HyperCube::Cube< Dim >::template Element< 0 > *c = SliceData::template HyperCubeTables< Dim , 1 , 0 >::OverlapElements[e.index];
						// [SANITY CHECK]
						//						if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i]->children + c[0].index )!=tree._isValidSpaceNode( tree._sNodes.treeNodes[i]->children + c[1].index ) ) ERROR_OUT( "Finer edges should both be valid or invalid" );
						if( !tree._isValidSpaceNode( tree._sNodes.treeNodes[i]->children + c[0].index ) || !tree._isValidSpaceNode( tree._sNodes.treeNodes[i]->children + c[1].index ) ) continue;

#ifdef NEW_CODE
						node_index_type cIndex1 = cSliceData.edgeIndices( tree._sNodes.treeNodes[i]->children + c[0].index )[_e.index];
						node_index_type cIndex2 = cSliceData.edgeIndices( tree._sNodes.treeNodes[i]->children + c[1].index )[_e.index];
#else // !NEW_CODE
						int cIndex1 = cSliceData.edgeIndices( tree._sNodes.treeNodes[i]->children + c[0].index )[_e.index];
						int cIndex2 = cSliceData.edgeIndices( tree._sNodes.treeNodes[i]->children + c[1].index )[_e.index];
#endif // NEW_CODE
						if( cSliceValues.edgeSet[cIndex1] != cSliceValues.edgeSet[cIndex2] )
						{
#ifdef NEW_HASH
							_Key key;
#else // !NEW_HASH
							long long key;
#endif // NEW_HASH
							if( cSliceValues.edgeSet[cIndex1] ) key = cSliceValues.edgeKeys[cIndex1];
							else                                key = cSliceValues.edgeKeys[cIndex2];
							pSliceValues.edgeKeys[pIndex] = key;
							pSliceValues.edgeSet[pIndex] = 1;
						}
						else if( cSliceValues.edgeSet[cIndex1] && cSliceValues.edgeSet[cIndex2] )
						{
#ifdef NEW_HASH
							_Key key1 = cSliceValues.edgeKeys[cIndex1] , key2 = cSliceValues.edgeKeys[cIndex2];
							pSliceValues.vertexPairKeyValues[ thread ].push_back( std::pair< _Key , _Key >( key1 , key2 ) );
#else // !NEW_HASH
							long long key1 = cSliceValues.edgeKeys[cIndex1] , key2 = cSliceValues.edgeKeys[cIndex2];
							pSliceValues.vertexPairKeyValues[ thread ].push_back( std::pair< long long , long long >( key1 , key2 ) );
#endif // NEW_HASH

							const TreeNode* node = tree._sNodes.treeNodes[i];
							LocalDepth _depth = depth;
							int _slice = slice;
							while( tree._isValidSpaceNode( node->parent ) && SliceData::template HyperCubeTables< Dim , 1 , 0 >::Overlap[e.index][(unsigned int)(node-node->parent->children) ] )
							{
								node = node->parent , _depth-- , _slice >>= 1;
								_SliceValues& _pSliceValues = slabValues[_depth].sliceValues(_slice);
#ifdef NEW_HASH
								_pSliceValues.vertexPairKeyValues[ thread ].push_back( std::pair< _Key , _Key >( key1 , key2 ) );
#else // !NEW_HASH
								_pSliceValues.vertexPairKeyValues[ thread ].push_back( std::pair< long long , long long >( key1 , key2 ) );
#endif // NEW_HASH
							}
						}
					}
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
#ifdef NEW_THREADS
	static void _CopyFinerXSliceIsoEdgeKeys( ThreadPool &tp , const FEMTree< Dim , Real >& tree , LocalDepth depth , int slab , std::vector< _SlabValues>& slabValues )
#else // !NEW_THREADS
	static void _CopyFinerXSliceIsoEdgeKeys( const FEMTree< Dim , Real >& tree , LocalDepth depth , int slab , std::vector< _SlabValues>& slabValues )
#endif // NEW_THREADS
	{
		_XSliceValues& pSliceValues  = slabValues[depth  ].xSliceValues(slab);
		_XSliceValues& cSliceValues0 = slabValues[depth+1].xSliceValues( (slab<<1)|0 );
		_XSliceValues& cSliceValues1 = slabValues[depth+1].xSliceValues( (slab<<1)|1 );
		typename SliceData::XSliceTableData& pSliceData  = pSliceValues.xSliceData;
		typename SliceData::XSliceTableData& cSliceData0 = cSliceValues0.xSliceData;
		typename SliceData::XSliceTableData& cSliceData1 = cSliceValues1.xSliceData;
#ifdef NEW_THREADS
		tp.parallel_for( tree._sNodesBegin(depth,slab) , tree._sNodesEnd(depth,slab) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=tree._sNodesBegin(depth,slab) ; i<tree._sNodesEnd(depth,slab) ; i++ )
#else // !NEW_CODE
		for( int i=tree._sNodesBegin(depth,slab) ; i<tree._sNodesEnd(depth,slab) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i] ) ) if( IsActiveNode< Dim >( tree._sNodes.treeNodes[i]->children ) )
			{
#ifdef NEW_THREADS
#else // !NEW_THREADS
				int thread = omp_get_thread_num();
#endif // NEW_THREADS
				typename SliceData::SquareCornerIndices& pIndices = pSliceData.edgeIndices( i );
				for( typename HyperCube::Cube< Dim-1 >::template Element< 0 > _c ; _c<HyperCube::Cube< Dim-1 >::template ElementNum< 0 >() ; _c++ )
				{
					typename HyperCube::Cube< Dim >::template Element< 1 > e( HyperCube::CROSS , _c.index );
#ifdef NEW_CODE
					node_index_type pIndex = pIndices[ _c.index ];
#else // !NEW_CODE
					int pIndex = pIndices[ _c.index ];
#endif // NEW_CODE
					if( !pSliceValues.edgeSet[pIndex] )
					{
						typename HyperCube::Cube< Dim >::template Element< 0 > c0( HyperCube::BACK , _c.index ) , c1( HyperCube::FRONT , _c.index );

						// [SANITY CHECK]
						//					if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i]->children + c0 )!=tree._isValidSpaceNode( tree._sNodes.treeNodes[i]->children + c1 ) ) ERROR_OUT( "Finer edges should both be valid or invalid" );
						if( !tree._isValidSpaceNode( tree._sNodes.treeNodes[i]->children + c0.index ) || !tree._isValidSpaceNode( tree._sNodes.treeNodes[i]->children + c1.index ) ) continue;

#ifdef NEW_CODE
						node_index_type cIndex0 = cSliceData0.edgeIndices( tree._sNodes.treeNodes[i]->children + c0.index )[_c.index];
						node_index_type cIndex1 = cSliceData1.edgeIndices( tree._sNodes.treeNodes[i]->children + c1.index )[_c.index];
#else // !NEW_CODE
						int cIndex0 = cSliceData0.edgeIndices( tree._sNodes.treeNodes[i]->children + c0.index )[_c.index];
						int cIndex1 = cSliceData1.edgeIndices( tree._sNodes.treeNodes[i]->children + c1.index )[_c.index];
#endif // NEW_CODE
						// If there's one zero-crossing along the edge
						if( cSliceValues0.edgeSet[cIndex0] != cSliceValues1.edgeSet[cIndex1] )
						{
#ifdef NEW_HASH
							_Key key;
#else // !NEW_HASH
							long long key;
#endif // NEW_HASH
							if( cSliceValues0.edgeSet[cIndex0] ) key = cSliceValues0.edgeKeys[cIndex0]; //, vPair = cSliceValues0.edgeVertexMap.find( key )->second;
							else                                 key = cSliceValues1.edgeKeys[cIndex1]; //, vPair = cSliceValues1.edgeVertexMap.find( key )->second;
							pSliceValues.edgeKeys[ pIndex ] = key;
							pSliceValues.edgeSet[ pIndex ] = 1;
						}
						// If there's are two zero-crossings along the edge
						else if( cSliceValues0.edgeSet[cIndex0] && cSliceValues1.edgeSet[cIndex1] )
						{
#ifdef NEW_HASH
							_Key key0 = cSliceValues0.edgeKeys[cIndex0] , key1 = cSliceValues1.edgeKeys[cIndex1];
							pSliceValues.vertexPairKeyValues[ thread ].push_back( std::pair< _Key , _Key >( key0 , key1 ) );
#else // !NEW_HASH
							long long key0 = cSliceValues0.edgeKeys[cIndex0] , key1 = cSliceValues1.edgeKeys[cIndex1];
							pSliceValues.vertexPairKeyValues[ thread ].push_back( std::pair< long long , long long >( key0 , key1 ) );
#endif // NEW_HASH
							const TreeNode* node = tree._sNodes.treeNodes[i];
							LocalDepth _depth = depth;
							int _slab = slab;
							while( tree._isValidSpaceNode( node->parent ) && SliceData::template HyperCubeTables< Dim , 1 , 0 >::Overlap[e.index][(unsigned int)(node-node->parent->children) ] )
							{
								node = node->parent , _depth-- , _slab>>= 1;
								_SliceValues& _pSliceValues = slabValues[_depth].sliceValues(_slab);
#ifdef NEW_HASH
								_pSliceValues.vertexPairKeyValues[ thread ].push_back( std::pair< _Key , _Key >( key0 , key1 ) );
#else // !NEW_HASH
								_pSliceValues.vertexPairKeyValues[ thread ].push_back( std::pair< long long , long long >( key0 , key1 ) );
#endif // NEW_HASH
							}
						}
					}
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
#ifdef NEW_THREADS
	static void _SetSliceIsoEdges( ThreadPool &tp , const FEMTree< Dim , Real >& tree , LocalDepth depth , int slice , std::vector< _SlabValues >& slabValues )
	{
		if( slice>0          ) _SetSliceIsoEdges( tp , tree , depth , slice , HyperCube::FRONT , slabValues );
		if( slice<(1<<depth) ) _SetSliceIsoEdges( tp , tree , depth , slice , HyperCube::BACK  , slabValues );
	}
#else // !NEW_THREADS
	static void _SetSliceIsoEdges( const FEMTree< Dim , Real >& tree , LocalDepth depth , int slice , std::vector< _SlabValues >& slabValues )
	{
		if( slice>0          ) _SetSliceIsoEdges( tree , depth , slice , HyperCube::FRONT , slabValues );
		if( slice<(1<<depth) ) _SetSliceIsoEdges( tree , depth , slice , HyperCube::BACK  , slabValues );
	}
#endif // NEW_THREADS
#ifdef NEW_THREADS
	static void _SetSliceIsoEdges( ThreadPool &tp , const FEMTree< Dim , Real >& tree , LocalDepth depth , int slice , HyperCube::Direction zDir , std::vector< _SlabValues >& slabValues )
#else // !NEW_THREADS
	static void _SetSliceIsoEdges( const FEMTree< Dim , Real >& tree , LocalDepth depth , int slice , HyperCube::Direction zDir , std::vector< _SlabValues >& slabValues )
#endif // NEW_THREADS
	{
		_SliceValues& sValues = slabValues[depth].sliceValues( slice );
#ifdef NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( tree._localToGlobal( depth ) );
#ifdef NEW_THREADS
		tp.parallel_for( tree._sNodesBegin(depth, slice-(zDir==HyperCube::BACK ? 0 : 1)) , tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=tree._sNodesBegin(depth, slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i<tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i++ )
#else // !NEW_CODE
		for( int i=tree._sNodesBegin(depth, slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i<tree._sNodesEnd(depth,slice-(zDir==HyperCube::BACK ? 0 : 1)) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i] ) )
			{
				int isoEdges[ 2 * HyperCube::MarchingSquares::MAX_EDGES ];
#ifdef NEW_THREADS
				ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#else // !NEW_THREADS
				ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
				TreeNode* leaf = tree._sNodes.treeNodes[i];
				if( !IsActiveNode< Dim >( leaf->children ) )
				{
#ifdef NEW_CODE
					node_index_type idx = i - sValues.sliceData.nodeOffset;
#else // !NEW_CODE
					int idx = i - sValues.sliceData.nodeOffset;
#endif // NEW_CODE
					const typename SliceData::SquareEdgeIndices& eIndices = sValues.sliceData.edgeIndices( leaf );
					const typename SliceData::SquareFaceIndices& fIndices = sValues.sliceData.faceIndices( leaf );
					unsigned char mcIndex = sValues.mcIndices[idx];
					if( !sValues.faceSet[ fIndices[0] ] )
					{
						neighborKey.getNeighbors( leaf );
						unsigned int xx = WindowIndex< IsotropicUIntPack< Dim , 3 > , IsotropicUIntPack< Dim , 1 > >::Index + (zDir==HyperCube::BACK ? -1 : 1);
						if( !IsActiveNode< Dim >( neighborKey.neighbors[ tree._localToGlobal( depth ) ].neighbors.data[xx] ) || !IsActiveNode< Dim >( neighborKey.neighbors[ tree._localToGlobal( depth ) ].neighbors.data[xx]->children ) )
						{
							_FaceEdges fe;
							fe.count = HyperCube::MarchingSquares::AddEdgeIndices( mcIndex , isoEdges );
							for( int j=0 ; j<fe.count ; j++ ) for( int k=0 ; k<2 ; k++ )
							{
								if( !sValues.edgeSet[ eIndices[ isoEdges[2*j+k] ] ] ) ERROR_OUT( "Edge not set: " , slice , " / " , 1<<depth );
								fe.edges[j][k] = sValues.edgeKeys[ eIndices[ isoEdges[2*j+k] ] ];
							}
							sValues.faceSet[ fIndices[0] ] = 1;
							sValues.faceEdges[ fIndices[0] ] = fe;

							TreeNode* node = leaf;
							LocalDepth _depth = depth;
							int _slice = slice;
							typename HyperCube::Cube< Dim >::template Element< Dim-1 > f( zDir , 0 );
							std::vector< _IsoEdge > edges;
							edges.resize( fe.count );
							for( int j=0 ; j<fe.count ; j++ ) edges[j] = fe.edges[j];
							while( tree._isValidSpaceNode( node->parent ) && SliceData::template HyperCubeTables< Dim , 2 , 0 >::Overlap[f.index][(unsigned int)(node-node->parent->children) ] )
							{
								node = node->parent , _depth-- , _slice >>= 1;
								if( IsActiveNode< Dim >( neighborKey.neighbors[ tree._localToGlobal( _depth ) ].neighbors.data[xx] ) && IsActiveNode< Dim >( neighborKey.neighbors[ tree._localToGlobal( _depth ) ].neighbors.data[xx]->children ) ) break;
#ifdef NEW_HASH
								_Key key = _VertexData::FaceIndex( node , f , tree._localToGlobal( tree._maxDepth ) );
#else // !NEW_HASH
								long long key = _VertexData::FaceIndex( node , f , tree._localToGlobal( tree._maxDepth ) );
#endif // NEW_HASH
								_SliceValues& _sValues = slabValues[_depth].sliceValues( _slice );
#ifdef NEW_THREADS
#ifdef NEW_HASH
								_sValues.faceEdgeKeyValues[ thread ].push_back( std::pair< _Key , std::vector< _IsoEdge > >( key , edges ) );
#else // !NEW_HASH
								_sValues.faceEdgeKeyValues[ thread ].push_back( std::pair< long long , std::vector< _IsoEdge > >( key , edges ) );
#endif // NEW_HASH
#else // !NEW_THREADS
#ifdef NEW_HASH
								_sValues.faceEdgeKeyValues[ omp_get_thread_num() ].push_back( std::pair< _Key , std::vector< _IsoEdge > >( key , edges ) );
#else // !NEW_HASH
								_sValues.faceEdgeKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::vector< _IsoEdge > >( key , edges ) );
#endif // NEW_HASH
#endif // NEW_THREADS
							}
						}
					}
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}
#ifdef NEW_THREADS
	static void _SetXSliceIsoEdges( ThreadPool &tp , const FEMTree< Dim , Real >& tree , LocalDepth depth , int slab , std::vector< _SlabValues >& slabValues )
#else // !NEW_THREADS
	static void _SetXSliceIsoEdges( const FEMTree< Dim , Real >& tree , LocalDepth depth , int slab , std::vector< _SlabValues >& slabValues )
#endif // NEW_THREADS
	{
		_SliceValues& bValues = slabValues[depth].sliceValues ( slab   );
		_SliceValues& fValues = slabValues[depth].sliceValues ( slab+1 );
		_XSliceValues& xValues = slabValues[depth].xSliceValues( slab   );

#ifdef NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( tp.threadNum() );
#else // !NEW_THREADS
		std::vector< ConstOneRingNeighborKey > neighborKeys( omp_get_max_threads() );
#endif // NEW_THREADS
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( tree._localToGlobal( depth ) );
#ifdef NEW_THREADS
		tp.parallel_for( tree._sNodesBegin(depth,slab) , tree._sNodesEnd(depth,slab) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=tree._sNodesBegin(depth,slab) ; i<tree._sNodesEnd(depth,slab) ; i++ )
#else // !NEW_CODE
		for( int i=tree._sNodesBegin(depth,slab) ; i<tree._sNodesEnd(depth,slab) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i] ) )
			{
				int isoEdges[ 2 * HyperCube::MarchingSquares::MAX_EDGES ];
#ifdef NEW_THREADS
				ConstOneRingNeighborKey& neighborKey = neighborKeys[ thread ];
#else // !NEW_THREADS
				ConstOneRingNeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
#endif // NEW_THREADS
				TreeNode* leaf = tree._sNodes.treeNodes[i];
				if( !IsActiveNode< Dim >( leaf->children ) )
				{
					const typename SliceData::SquareCornerIndices& cIndices = xValues.xSliceData.edgeIndices( leaf );
					const typename SliceData::SquareEdgeIndices& eIndices = xValues.xSliceData.faceIndices( leaf );
					unsigned char mcIndex = ( bValues.mcIndices[ i - bValues.sliceData.nodeOffset ] ) | ( fValues.mcIndices[ i - fValues.sliceData.nodeOffset ]<<4 );
					{
						neighborKey.getNeighbors( leaf );
						// Iterate over the edges on the back
						for( typename HyperCube::Cube< Dim-1 >::template Element< 1 > _e ; _e<HyperCube::Cube< Dim-1 >::template ElementNum< 1 >() ; _e++ )
						{
							typename HyperCube::Cube< Dim >::template Element< 2 > f( HyperCube::CROSS , _e.index );
							unsigned char _mcIndex = HyperCube::Cube< Dim >::template ElementMCIndex< 2 >( f , mcIndex );

							unsigned int xx = SliceData::template HyperCubeTables< Dim , 2 >::CellOffsetAntipodal[f.index];
							if(	!xValues.faceSet[ eIndices[_e.index] ] && ( !IsActiveNode< Dim >( neighborKey.neighbors[ tree._localToGlobal( depth ) ].neighbors.data[xx] ) || !IsActiveNode< Dim >( neighborKey.neighbors[ tree._localToGlobal( depth ) ].neighbors.data[xx]->children ) ) )
							{
								_FaceEdges fe;
								fe.count = HyperCube::MarchingSquares::AddEdgeIndices( _mcIndex , isoEdges );
								for( int j=0 ; j<fe.count ; j++ ) for( int k=0 ; k<2 ; k++ )
								{
									typename HyperCube::Cube< Dim >::template Element< 1 > e( f , typename HyperCube::Cube< Dim-1 >::template Element< 1 >( isoEdges[2*j+k] ) );
									HyperCube::Direction dir ; unsigned int coIndex;
									e.factor( dir , coIndex );
									if( dir==HyperCube::CROSS ) // Cross-edge
									{
#ifdef NEW_CODE
										node_index_type idx = cIndices[ coIndex ];
#else // !NEW_CODE
										int idx = cIndices[ coIndex ];
#endif // NEW_CODE
										if( !xValues.edgeSet[ idx ] ) ERROR_OUT( "Edge not set: " , slab , " / " , 1<<depth );
										fe.edges[j][k] = xValues.edgeKeys[ idx ];
									}
									else
									{
										const _SliceValues& sValues = dir==HyperCube::BACK ? bValues : fValues;
#ifdef NEW_CODE
										node_index_type idx = sValues.sliceData.edgeIndices(i)[ coIndex ];
#else // !NEW_CODE
										int idx = sValues.sliceData.edgeIndices(i)[ coIndex ];
#endif // NEW_CODE
										if( !sValues.edgeSet[ idx ] ) ERROR_OUT( "Edge not set: " , slab , " / " , 1<<depth );
										fe.edges[j][k] = sValues.edgeKeys[ idx ];
									}
								}
								xValues.faceSet[ eIndices[_e.index] ] = 1;
								xValues.faceEdges[ eIndices[_e.index] ] = fe;

								TreeNode* node = leaf;
								LocalDepth _depth = depth;
								int _slab = slab;
								std::vector< _IsoEdge > edges;
								edges.resize( fe.count );
								for( int j=0 ; j<fe.count ; j++ ) edges[j] = fe.edges[j];
								while( tree._isValidSpaceNode( node->parent ) && SliceData::template HyperCubeTables< Dim , 2 , 0 >::Overlap[f.index][(unsigned int)(node-node->parent->children) ] )
								{
									node = node->parent , _depth-- , _slab >>= 1;
									if( IsActiveNode< Dim >( neighborKey.neighbors[ tree._localToGlobal( _depth ) ].neighbors.data[xx] ) && IsActiveNode< Dim >( neighborKey.neighbors[ tree._localToGlobal( _depth ) ].neighbors.data[xx]->children ) ) break;
#ifdef NEW_HASH
									_Key key = _VertexData::FaceIndex( node , f , tree._localToGlobal( tree._maxDepth ) );
#else // !NEW_HASH
									long long key = _VertexData::FaceIndex( node , f , tree._localToGlobal( tree._maxDepth ) );
#endif // NEW_HASH
									_XSliceValues& _xValues = slabValues[_depth].xSliceValues( _slab );
#ifdef NEW_THREADS
#ifdef NEW_HASH
									_xValues.faceEdgeKeyValues[ thread ].push_back( std::pair< _Key , std::vector< _IsoEdge > >( key , edges ) );
#else // !NEW_HASH
									_xValues.faceEdgeKeyValues[ thread ].push_back( std::pair< long long , std::vector< _IsoEdge > >( key , edges ) );
#endif // NEW_HASH
#else // !NEW_THREADS
#ifdef NEW_HASH
									_xValues.faceEdgeKeyValues[ omp_get_thread_num() ].push_back( std::pair< _Key , std::vector< _IsoEdge > >( key , edges ) );
#else // !NEW_HASH
									_xValues.faceEdgeKeyValues[ omp_get_thread_num() ].push_back( std::pair< long long , std::vector< _IsoEdge > >( key , edges ) );
#endif // NEW_HASH
#endif // NEW_THREADS
								}
							}
						}
					}
				}
			}		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}

#ifdef NEW_CODE
#ifdef NEW_THREADS
	static void _SetIsoSurface( ThreadPool &tp , const FEMTree< Dim , Real >& tree , LocalDepth depth , int offset , const _SliceValues& bValues , const _SliceValues& fValues , const _XSliceValues& xValues , CoredMeshData< Vertex , node_index_type >& mesh , bool polygonMesh , bool addBarycenter , node_index_type& vOffset , bool flipOrientation )
#else // !NEW_THREADS
	static void _SetIsoSurface( const FEMTree< Dim , Real >& tree , LocalDepth depth , int offset , const _SliceValues& bValues , const _SliceValues& fValues , const _XSliceValues& xValues , CoredMeshData< Vertex , node_index_type >& mesh , bool polygonMesh , bool addBarycenter , node_index_type& vOffset , bool flipOrientation )
#endif // NEW_THREADS
#else // !NEW_CODE
	static void _SetIsoSurface( const FEMTree< Dim , Real >& tree , LocalDepth depth , int offset , const _SliceValues& bValues , const _SliceValues& fValues , const _XSliceValues& xValues , CoredMeshData< Vertex >& mesh , bool polygonMesh , bool addBarycenter , int& vOffset , bool flipOrientation )
#endif // NEW_CODE
	{
#ifdef NEW_CODE
		std::vector< std::pair< node_index_type , Vertex > > polygon;
#else // !NEW_CODE
		std::vector< std::pair< int , Vertex > > polygon;
#endif // NEW_CODE
#ifdef NEW_THREADS
		std::vector< std::vector< _IsoEdge > > edgess( tp.threadNum() );
		tp.parallel_for( tree._sNodesBegin(depth,offset) , tree._sNodesEnd(depth,offset) , [&]( unsigned int thread , size_t i )
#else // !NEW_THREADS
		std::vector< std::vector< _IsoEdge > > edgess( omp_get_max_threads() );
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=tree._sNodesBegin(depth,offset) ; i<tree._sNodesEnd(depth,offset) ; i++ )
#else // !NEW_CODE
		for( int i=tree._sNodesBegin(depth,offset) ; i<tree._sNodesEnd(depth,offset) ; i++ )
#endif // NEW_CODE
#endif // NEW_THREADS
		{
			if( tree._isValidSpaceNode( tree._sNodes.treeNodes[i] ) )
			{
#ifdef NEW_THREADS
				std::vector< _IsoEdge >& edges = edgess[ thread ];
#else // !NEW_THREADS
				std::vector< _IsoEdge >& edges = edgess[ omp_get_thread_num() ];
#endif // NEW_THREADS
				TreeNode* leaf = tree._sNodes.treeNodes[i];
				int res = 1<<depth;
				LocalDepth d ; LocalOffset off;
				tree._localDepthAndOffset( leaf , d , off );
				bool inBounds = off[0]>=0 && off[0]<res && off[1]>=0 && off[1]<res && off[2]>=0 && off[2]<res;
				if( inBounds && !IsActiveNode< Dim >( leaf->children ) )
				{
					edges.clear();
					unsigned char mcIndex = ( bValues.mcIndices[ i - bValues.sliceData.nodeOffset ] ) | ( fValues.mcIndices[ i - fValues.sliceData.nodeOffset ]<<4 );
					// [WARNING] Just because the node looks empty doesn't mean it doesn't get eges from finer neighbors
					{
						// Gather the edges from the faces (with the correct orientation)
						for( typename HyperCube::Cube< Dim >::template Element< Dim-1 > f ; f<HyperCube::Cube< Dim >::template ElementNum< Dim-1 >() ; f++ )
						{
							int flip = HyperCube::Cube< Dim >::IsOriented( f ) ? 0 : 1;
							HyperCube::Direction fDir = f.direction();
							if( fDir==HyperCube::BACK || fDir==HyperCube::FRONT )
							{
								const _SliceValues& sValues = (fDir==HyperCube::BACK) ? bValues : fValues;
#ifdef NEW_CODE
								node_index_type fIdx = sValues.sliceData.faceIndices(i)[0];
#else // !NEW_CODE
								int fIdx = sValues.sliceData.faceIndices(i)[0];
#endif // NEW_CODE
								if( sValues.faceSet[fIdx] )
								{
									const _FaceEdges& fe = sValues.faceEdges[ fIdx ];
									for( int j=0 ; j<fe.count ; j++ ) edges.push_back( _IsoEdge( fe.edges[j][flip] , fe.edges[j][1-flip] ) );
								}
								else
								{
#ifdef NEW_HASH
									_Key key = _VertexData::FaceIndex( leaf , f , tree._localToGlobal( tree._maxDepth ) );
									typename std::unordered_map< _Key , std::vector< _IsoEdge > , typename _Key::Hasher >::const_iterator iter = sValues.faceEdgeMap.find(key);
#else // !NEW_HASH
									long long key = _VertexData::FaceIndex( leaf , f , tree._localToGlobal( tree._maxDepth ) );
									typename std::unordered_map< long long , std::vector< _IsoEdge > >::const_iterator iter = sValues.faceEdgeMap.find(key);
#endif // NEW_HASH
									if( iter!=sValues.faceEdgeMap.end() )
									{
										const std::vector< _IsoEdge >& _edges = iter->second;
										for( size_t j=0 ; j<_edges.size() ; j++ ) edges.push_back( _IsoEdge( _edges[j][flip] , _edges[j][1-flip] ) );
									}
									else ERROR_OUT( "Invalid faces: " , i , "  " , fDir==HyperCube::BACK ? "back" : ( fDir==HyperCube::FRONT ? "front" : ( fDir==HyperCube::CROSS ? "cross" : "unknown" ) ) );
								}
							}
							else
							{
#ifdef NEW_CODE
								node_index_type fIdx = xValues.xSliceData.faceIndices(i)[ f.coIndex() ];
#else // !NEW_CODE
								int fIdx = xValues.xSliceData.faceIndices(i)[ f.coIndex() ];
#endif // NEW_CODE
								if( xValues.faceSet[fIdx] )
								{
									const _FaceEdges& fe = xValues.faceEdges[ fIdx ];
									for( int j=0 ; j<fe.count ; j++ ) edges.push_back( _IsoEdge( fe.edges[j][flip] , fe.edges[j][1-flip] ) );
								}
								else
								{
#ifdef NEW_HASH
									_Key key = _VertexData::FaceIndex( leaf , f , tree._localToGlobal( tree._maxDepth ) );
									typename std::unordered_map< _Key , std::vector< _IsoEdge > , typename _Key::Hasher >::const_iterator iter = xValues.faceEdgeMap.find(key);
#else // !NEW_HASH
									long long key = _VertexData::FaceIndex( leaf , f , tree._localToGlobal( tree._maxDepth ) );
									typename std::unordered_map< long long , std::vector< _IsoEdge > >::const_iterator iter = xValues.faceEdgeMap.find(key);
#endif // NEW_HASH
									if( iter!=xValues.faceEdgeMap.end() )
									{
										const std::vector< _IsoEdge >& _edges = iter->second;
										for( size_t j=0 ; j<_edges.size() ; j++ ) edges.push_back( _IsoEdge( _edges[j][flip] , _edges[j][1-flip] ) );
									}
									else ERROR_OUT( "Invalid faces: " , i , "  " ,  fDir==HyperCube::BACK ? "back" : ( fDir==HyperCube::FRONT ? "front" : ( fDir==HyperCube::CROSS ? "cross" : "unknown" ) ) );
								}
							}
						}
						// Get the edge loops
#ifdef NEW_HASH
						std::vector< std::vector< _Key > > loops;
#else // !NEW_HASH
						std::vector< std::vector< long long  > > loops;
#endif // NEW_HASH
						while( edges.size() )
						{
							loops.resize( loops.size()+1 );
							_IsoEdge edge = edges.back();
							edges.pop_back();
#ifdef NEW_HASH
							_Key start = edge[0] , current = edge[1];
#else // !NEW_HASH
							long long start = edge[0] , current = edge[1];
#endif // NEW_HASH
							while( current!=start )
							{
								int idx;
								for( idx=0 ; idx<(int)edges.size() ; idx++ ) if( edges[idx][0]==current ) break;
								if( idx==edges.size() )
								{
#ifdef NEW_HASH
									typename std::unordered_map< _Key , _Key , typename _Key::Hasher >::const_iterator iter;
#else // !NEW_HASH
									typename std::unordered_map< long long, long long >::const_iterator iter;
#endif // NEW_HASH
									if     ( (iter=bValues.vertexPairMap.find(current))!=bValues.vertexPairMap.end() ) loops.back().push_back( current ) , current = iter->second;
									else if( (iter=fValues.vertexPairMap.find(current))!=fValues.vertexPairMap.end() ) loops.back().push_back( current ) , current = iter->second;
									else if( (iter=xValues.vertexPairMap.find(current))!=xValues.vertexPairMap.end() ) loops.back().push_back( current ) , current = iter->second;
									else
									{
										LocalDepth d ; LocalOffset off;
										tree._localDepthAndOffset( leaf , d , off );
#ifdef NEW_HASH
										ERROR_OUT( "Failed to close loop [" , d-1 , ": " , off[0] , " " , off[1] , " " , off[2] , "] | (" , i , "): " , current.to_string() );
#else // !NEW_HASH
										ERROR_OUT( "Failed to close loop [" , d-1 , ": " , off[0] , " " , off[1] , " " , off[2] , "] | (" , i , "): " , current );
#endif // NEW_HASH
									}
								}
								else
								{
									loops.back().push_back( current );
									current = edges[idx][1];
									edges[idx] = edges.back() , edges.pop_back();
								}
							}
							loops.back().push_back( start );
						}
						// Add the loops to the mesh
						for( size_t j=0 ; j<loops.size() ; j++ )
						{
#ifdef NEW_CODE
							std::vector< std::pair< node_index_type , Vertex > > polygon( loops[j].size() );
#else // !NEW_CODE
							std::vector< std::pair< int , Vertex > > polygon( loops[j].size() );
#endif // NEW_CODE
							for( size_t k=0 ; k<loops[j].size() ; k++ )
							{
#ifdef NEW_HASH
								_Key key = loops[j][k];
								typename std::unordered_map< _Key , std::pair< node_index_type , Vertex > , typename _Key::Hasher >::const_iterator iter;
#else // !NEW_HASH
								long long key = loops[j][k];
#ifdef NEW_CODE
								typename std::unordered_map< long long, std::pair< node_index_type , Vertex > , typename _Key::Hasher >::const_iterator iter;
#else // !NEW_CODE
								typename std::unordered_map< long long, std::pair< int, Vertex > >::const_iterator iter;
#endif // NEW_CODE
#endif // NEW_HASH
								size_t kk = flipOrientation ? loops[j].size()-1-k : k;
								if     ( ( iter=bValues.edgeVertexMap.find( key ) )!=bValues.edgeVertexMap.end() ) polygon[kk] = iter->second;
								else if( ( iter=fValues.edgeVertexMap.find( key ) )!=fValues.edgeVertexMap.end() ) polygon[kk] = iter->second;
								else if( ( iter=xValues.edgeVertexMap.find( key ) )!=xValues.edgeVertexMap.end() ) polygon[kk] = iter->second;
								else ERROR_OUT( "Couldn't find vertex in edge map" );
							}
#ifdef NEW_THREADS
							_AddIsoPolygons( thread , mesh , polygon , polygonMesh , addBarycenter , vOffset );
#else // !NEW_THREADS
							_AddIsoPolygons( mesh , polygon , polygonMesh , addBarycenter , vOffset );
#endif // NEW_THREADS
						}
					}
				}
			}
		}
#ifdef NEW_THREADS
		);
#endif // NEW_THREADS
	}

	template< unsigned int WeightDegree , typename Data , unsigned int DataSig >
	static bool _GetIsoVertex( const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey , ConstPointSupportKey< IsotropicUIntPack< Dim , FEMSignature< DataSig >::Degree > >& dataKey , const TreeNode* node , typename HyperCube::template Cube< Dim-1 >::template Element< 1 > _e , HyperCube::Direction zDir , const _SliceValues& sValues , Vertex& vertex , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
	{
		static const unsigned int DataDegree = FEMSignature< DataSig >::Degree;
		Point< Real , Dim > position;
		int c0 , c1;
		const typename HyperCube::Cube< Dim-1 >::template Element< 0 > *_c = SliceData::template HyperCubeTables< Dim-1 , 1 , 0 >::OverlapElements[_e.index];
		c0 = _c[0].index , c1 = _c[1].index;

		bool nonLinearFit = sValues.cornerGradients!=NullPointer( Point< Real , Dim > );
		const typename SliceData::SquareCornerIndices& idx = sValues.sliceData.cornerIndices( node );
		Real x0 = sValues.cornerValues[idx[c0]] , x1 = sValues.cornerValues[idx[c1]];
		Point< Real , Dim > s;
		Real start , width;
		tree._startAndWidth( node , s , width );
		int o;
		{
			const HyperCube::Direction* dirs = SliceData::template HyperCubeTables< Dim-1 , 1 >::Directions[ _e.index ];
			for( int d=0 ; d<Dim-1 ; d++ ) if( dirs[d]==HyperCube::CROSS )
			{
				o = d;
				start = s[d];
				for( int dd=1 ; dd<Dim-1 ; dd++ ) position[(d+dd)%(Dim-1)] = s[(d+dd)%(Dim-1)] + width * ( dirs[(d+dd)%(Dim-1)]==HyperCube::BACK ? 0 : 1 );
			}
		}
		position[ Dim-1 ] = s[Dim-1] + width * ( zDir==HyperCube::BACK ? 0 : 1 );

		double averageRoot;
		bool rootFound = false;
		if( nonLinearFit )
		{
			double dx0 = sValues.cornerGradients[idx[c0]][o] * width , dx1 = sValues.cornerGradients[idx[c1]][o] * width;

			// The scaling will turn the Hermite Spline into a quadratic
			double scl = (x1-x0) / ( (dx1+dx0 ) / 2 );
			dx0 *= scl , dx1 *= scl;

			// Hermite Spline
			Polynomial< 2 > P;
			P.coefficients[0] = x0;
			P.coefficients[1] = dx0;
			P.coefficients[2] = 3*(x1-x0)-dx1-2*dx0;

			double roots[2];
			int rCount = 0 , rootCount = P.getSolutions( isoValue , roots , 0 );
			averageRoot = 0;
			for( int i=0 ; i<rootCount ; i++ ) if( roots[i]>=0 && roots[i]<=1 ) averageRoot += roots[i] , rCount++;
			if( rCount ) rootFound = true;
			averageRoot /= rCount;
		}
		if( !rootFound )
		{
			// We have a linear function L, with L(0) = x0 and L(1) = x1
			// => L(t) = x0 + t * (x1-x0)
			// => L(t) = isoValue <=> t = ( isoValue - x0 ) / ( x1 - x0 )
			if( x0==x1 ) ERROR_OUT( "Not a zero-crossing root: " , x0 , " " , x1 );
			averageRoot = ( isoValue - x0 ) / ( x1 - x0 );
		}
		if( averageRoot<=0 || averageRoot>=1 )
		{
			WARN( "Bad average root: " , averageRoot , "\t(" , x0 ," " , x1 ,") (" , isoValue , ")" );
			if( averageRoot<0 ) averageRoot = 0;
			if( averageRoot>1 ) averageRoot = 1;
		}
		position[o] = Real( start + width*averageRoot );
		Real depth = (Real)1.;
		Data dataValue;
		if( densityWeights )
		{
			Real weight;
			tree._getSampleDepthAndWeight( *densityWeights , node , position , weightKey , depth , weight );
		}
		if( data )
		{
			if( DataDegree==0 ) 
			{
				Point< Real , 3 > center( s[0] + width/2 , s[1] + width/2 , s[2] + width/2 );
				dataValue = tree.template _evaluate< ProjectiveData< Data , Real > , SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > , 0 >( *data , center , *pointEvaluator , dataKey ).value();
			}
			else dataValue = tree.template _evaluate< ProjectiveData< Data , Real > , SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > , 0 >( *data , position , *pointEvaluator , dataKey ).value();
		}
		SetVertex( vertex , position , depth , dataValue );
		return true;
	}
	template< unsigned int WeightDegree , typename Data , unsigned int DataSig >
	static bool _GetIsoVertex( const FEMTree< Dim , Real >& tree , typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , Real isoValue , ConstPointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey , ConstPointSupportKey< IsotropicUIntPack< Dim , FEMSignature< DataSig >::Degree > >& dataKey , const TreeNode* node , typename HyperCube::template Cube< Dim-1 >::template Element< 0 > _c , const _SliceValues& bValues , const _SliceValues& fValues , Vertex& vertex , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex )
	{
		static const unsigned int DataDegree = FEMSignature< DataSig >::Degree;
		Point< Real , Dim > position;

		bool nonLinearFit = bValues.cornerGradients!=NullPointer( Point< Real , Dim > ) && fValues.cornerGradients!=NullPointer( Point< Real , Dim > );
		const typename SliceData::SquareCornerIndices& idx0 = bValues.sliceData.cornerIndices( node );
		const typename SliceData::SquareCornerIndices& idx1 = fValues.sliceData.cornerIndices( node );
		Real x0 = bValues.cornerValues[ idx0[_c.index] ] , x1 = fValues.cornerValues[ idx1[_c.index] ];
		Point< Real , Dim > s;
		Real start , width;
		tree._startAndWidth( node , s , width );
		start = s[2];
		int x , y;
		{
			const HyperCube::Direction* xx = SliceData::template HyperCubeTables< Dim-1 , 0 >::Directions[ _c.index ];
			x = xx[0]==HyperCube::BACK ? 0 : 1 , y = xx[1]==HyperCube::BACK ? 0 : 1;
		}

		position[0] = s[0] + width*x;
		position[1] = s[1] + width*y;

		double averageRoot;
		bool rootFound = false;

		if( nonLinearFit )
		{
			double dx0 = bValues.cornerGradients[ idx0[_c.index] ][2] * width , dx1 = fValues.cornerGradients[ idx1[_c.index] ][2] * width;
			// The scaling will turn the Hermite Spline into a quadratic
			double scl = (x1-x0) / ( (dx1+dx0 ) / 2 );
			dx0 *= scl , dx1 *= scl;

			// Hermite Spline
			Polynomial< 2 > P;
			P.coefficients[0] = x0;
			P.coefficients[1] = dx0;
			P.coefficients[2] = 3*(x1-x0)-dx1-2*dx0;

			double roots[2];
			int rCount = 0 , rootCount = P.getSolutions( isoValue , roots , 0 );
			averageRoot = 0;
			for( int i=0 ; i<rootCount ; i++ ) if( roots[i]>=0 && roots[i]<=1 ) averageRoot += roots[i] , rCount++;
			if( rCount ) rootFound = true;
			averageRoot /= rCount;
		}
		if( !rootFound )
		{
			// We have a linear function L, with L(0) = x0 and L(1) = x1
			// => L(t) = x0 + t * (x1-x0)
			// => L(t) = isoValue <=> t = ( isoValue - x0 ) / ( x1 - x0 )
			if( x0==x1 ) ERROR_OUT( "Not a zero-crossing root: " , x0 , " " , x1 );
			averageRoot = ( isoValue - x0 ) / ( x1 - x0 );
		}
		if( averageRoot<=0 || averageRoot>=1 )
		{
			WARN( "Bad average root: " , averageRoot , "\t(" , x0 ," " , x1 ,") (" , isoValue , ")" );
			if( averageRoot<0 ) averageRoot = 0;
			if( averageRoot>1 ) averageRoot = 1;
		}
		position[2] = Real( start + width*averageRoot );
		Real depth = (Real)1.;
		Data dataValue;
		if( densityWeights )
		{
			Real weight;
			tree._getSampleDepthAndWeight( *densityWeights , node , position , weightKey , depth , weight );
		}
		if( data )
		{
			if( DataDegree==0 ) 
			{
				Point< Real , 3 > center( s[0] + width/2 , s[1] + width/2 , s[2] + width/2 );
				dataValue = tree.template _evaluate< ProjectiveData< Data , Real > , SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > , 0 >( *data , center , *pointEvaluator , dataKey ).value();
			}
			else dataValue = tree.template _evaluate< ProjectiveData< Data , Real > , SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > , 0 >( *data , position , *pointEvaluator , dataKey ).value();
		}
		SetVertex( vertex , position , depth , dataValue );
		return true;
	}

#ifdef NEW_CODE
#ifdef NEW_THREADS
	static unsigned int _AddIsoPolygons( unsigned int thread , CoredMeshData< Vertex , node_index_type >& mesh , std::vector< std::pair< node_index_type , Vertex > >& polygon , bool polygonMesh , bool addBarycenter , node_index_type &vOffset )
#else // !NEW_THREADS
	static unsigned int _AddIsoPolygons( CoredMeshData< Vertex , node_index_type >& mesh , std::vector< std::pair< node_index_type , Vertex > >& polygon , bool polygonMesh , bool addBarycenter , node_index_type &vOffset )
#endif // NEW_THREADS
#else // !NEW_CODE
	static int _AddIsoPolygons( CoredMeshData< Vertex >& mesh , std::vector< std::pair< int , Vertex > >& polygon , bool polygonMesh , bool addBarycenter , int& vOffset )
#endif // NEW_CODE
	{
		if( polygonMesh )
		{
#ifdef NEW_CODE
			std::vector< node_index_type > vertices( polygon.size() );
			for( unsigned int i=0 ; i<polygon.size() ; i++ ) vertices[i] = polygon[polygon.size()-1-i].first;
#else // !NEW_CODE
			std::vector< int > vertices( polygon.size() );
			for( int i=0 ; i<(int)polygon.size() ; i++ ) vertices[i] = polygon[polygon.size()-1-i].first;
#endif // NEW_CODE
#ifdef NEW_THREADS
			mesh.addPolygon_s( thread , vertices );
#else // !NEW_THREADS
			mesh.addPolygon_s( vertices );
#endif // NEW_THREADS
			return 1;
		}
		if( polygon.size()>3 )
		{
			bool isCoplanar = false;
#ifdef NEW_CODE
			std::vector< node_index_type > triangle( 3 );
#else // !NEW_CODE
			std::vector< int > triangle( 3 );
#endif // NEW_CODE

			if( addBarycenter )
#ifdef NEW_CODE
				for( unsigned int i=0 ; i<polygon.size() ; i++ ) for( unsigned int j=0 ; j<i ; j++ )
					if( (i+1)%polygon.size()!=j && (j+1)%polygon.size()!=i )
					{
						Vertex v1 = polygon[i].second , v2 = polygon[j].second;
						for( int k=0 ; k<3 ; k++ ) if( v1.point[k]==v2.point[k] ) isCoplanar = true;
					}
#else // !NEW_CODE
				for( int i=0 ; i<(int)polygon.size() ; i++ )
					for( int j=0 ; j<i ; j++ )
						if( (i+1)%polygon.size()!=j && (j+1)%polygon.size()!=i )
						{
							Vertex v1 = polygon[i].second , v2 = polygon[j].second;
							for( int k=0 ; k<3 ; k++ ) if( v1.point[k]==v2.point[k] ) isCoplanar = true;
						}
#endif // NEW_CODE
			if( isCoplanar )
			{
				Vertex c;
				c *= 0;
#ifdef NEW_CODE
				for( unsigned int i=0 ; i<polygon.size() ; i++ ) c += polygon[i].second;
#else // !NEW_CODE
				for( int i=0 ; i<(int)polygon.size() ; i++ ) c += polygon[i].second;
#endif // NEW_CODE
				c /= ( typename Vertex::Real )polygon.size();
#ifdef NEW_CODE
				node_index_type cIdx;
#else // !NEW_CODE
				int cIdx;
#endif // NEW_CODE
#ifdef NEW_THREADS
				{
					std::lock_guard< std::mutex > lock( _pointInsertionMutex );
					cIdx = mesh.addOutOfCorePoint( c );
					vOffset++;
				}
#else // !NEW_THREADS
#pragma omp critical (add_barycenter_point_access)
				{
					cIdx = mesh.addOutOfCorePoint( c );
					vOffset++;
				}
#endif // NEW_THREADS
#ifdef NEW_CODE
				for( unsigned i=0 ; i<polygon.size() ; i++ )
#else // !NEW_CODE
				for( int i=0 ; i<(int)polygon.size() ; i++ )
#endif // NEW_CODE
				{
					triangle[0] = polygon[ i                  ].first;
					triangle[1] = cIdx;
					triangle[2] = polygon[(i+1)%polygon.size()].first;
#ifdef NEW_THREADS
					mesh.addPolygon_s( thread , triangle );
#else // !NEW_THREADS
					mesh.addPolygon_s( triangle );
#endif // NEW_THREADS
				}
#ifdef NEW_CODE
				return (unsigned int)polygon.size();
#else // !NEW_CODE
				return (int)polygon.size();
#endif // NEW_CODE
			}
			else
			{
				std::vector< Point< Real , Dim > > vertices( polygon.size() );
#ifdef NEW_CODE
				for( unsigned int i=0 ; i<polygon.size() ; i++ ) vertices[i] = polygon[i].second.point;
				std::vector< TriangleIndex< node_index_type > > triangles = MinimalAreaTriangulation< node_index_type , Real , Dim >( ( ConstPointer( Point< Real , Dim > ) )GetPointer( vertices ) , (node_index_type)vertices.size() );
#else // !NEW_CODE
				for( int i=0 ; i<(int)polygon.size() ; i++ ) vertices[i] = polygon[i].second.point;
				std::vector< TriangleIndex > triangles = MinimalAreaTriangulation< Real , Dim >( ( ConstPointer( Point< Real , Dim > ) )GetPointer( vertices ) , vertices.size() );
#endif // NEW_CODE
				if( triangles.size()!=polygon.size()-2 ) ERROR_OUT( "Minimal area triangulation failed:" , triangles.size() , " != " , polygon.size()-2 );
#ifdef NEW_CODE
				for( unsigned int i=0 ; i<triangles.size() ; i++ )
#else // !NEW_CODE
				for( int i=0 ; i<(int)triangles.size() ; i++ )
#endif // NEW_CODE
				{
					for( int j=0 ; j<3 ; j++ ) triangle[2-j] = polygon[ triangles[i].idx[j] ].first;
#ifdef NEW_THREADS
					mesh.addPolygon_s( thread , triangle );
#else // !NEW_THREADS
					mesh.addPolygon_s( triangle );
#endif // NEW_THREADS
				}
			}
		}
		else if( polygon.size()==3 )
		{
#ifdef NEW_CODE
			std::vector< node_index_type > vertices( 3 );
#else // !NEW_CODE
			std::vector< int > vertices( 3 );
#endif // NEW_CODE
			for( int i=0 ; i<3 ; i++ ) vertices[2-i] = polygon[i].first;
#ifdef NEW_THREADS
			mesh.addPolygon_s( thread , vertices );
#else // !NEW_THREADS
			mesh.addPolygon_s( vertices );
#endif // NEW_THREADS
		}
#ifdef NEW_CODE
		return (unsigned int)polygon.size()-2;
#else // !NEW_CODE
		return (int)polygon.size()-2;
#endif // NEW_CODE
	}
public:
	struct IsoStats
	{
		double cornersTime , verticesTime , edgesTime , surfaceTime;
		double copyFinerTime , setTableTime;
		IsoStats( void ) : cornersTime(0) , verticesTime(0) , edgesTime(0) , surfaceTime(0) , copyFinerTime(0) , setTableTime(0) {;}
		std::string toString( void ) const
		{
			std::stringstream stream;
			stream << "Corners / Vertices / Edges / Surface / Set Table / Copy Finer: ";
			stream << std::fixed << std::setprecision(1) << cornersTime << " / " << verticesTime << " / " << edgesTime << " / " << surfaceTime << " / " << setTableTime << " / " << copyFinerTime;
			stream << " (s)";
			return stream.str();
		}
	};
#ifdef NEW_THREADS
	template< typename Data , unsigned int ... FEMSigs , unsigned int WeightDegree , unsigned int DataSig >
	static IsoStats Extract( UIntPack< FEMSigs ... > , UIntPack< WeightDegree > , UIntPack< DataSig > , const FEMTree< Dim , Real >& tree , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , const DenseNodeData< Real , UIntPack< FEMSigs ... > >& coefficients , Real isoValue , CoredMeshData< Vertex , node_index_type >& mesh , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex , bool nonLinearFit , bool addBarycenter , bool polygonMesh , bool flipOrientation )
	{
		ThreadPool tp;
		return Extract( UIntPack< FEMSigs ... >() , UIntPack< WeightDegree >() , UIntPack< DataSig >() , tp , tree , densityWeights , data , coefficients , isoValue , mesh , SetVertex , nonLinearFit , addBarycenter , polygonMesh , flipOrientation );
	}
#endif // NEW_THREADS
	template< typename Data , unsigned int ... FEMSigs , unsigned int WeightDegree , unsigned int DataSig >
#ifdef NEW_CODE
#ifdef NEW_THREADS
	static IsoStats Extract( UIntPack< FEMSigs ... > , UIntPack< WeightDegree > , UIntPack< DataSig > , ThreadPool &tp , const FEMTree< Dim , Real >& tree , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , const DenseNodeData< Real , UIntPack< FEMSigs ... > >& coefficients , Real isoValue , CoredMeshData< Vertex , node_index_type >& mesh , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex , bool nonLinearFit , bool addBarycenter , bool polygonMesh , bool flipOrientation )
#else // !NEW_THREADS
	static IsoStats Extract( UIntPack< FEMSigs ... > , UIntPack< WeightDegree > , UIntPack< DataSig > , const FEMTree< Dim , Real >& tree , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , const DenseNodeData< Real , UIntPack< FEMSigs ... > >& coefficients , Real isoValue , CoredMeshData< Vertex , node_index_type >& mesh , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex , bool nonLinearFit , bool addBarycenter , bool polygonMesh , bool flipOrientation )
#endif // NEW_THREADS
#else // !NEW_CODE
	static IsoStats Extract( UIntPack< FEMSigs ... > , UIntPack< WeightDegree > , UIntPack< DataSig > , const FEMTree< Dim , Real >& tree , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > >* data , const DenseNodeData< Real , UIntPack< FEMSigs ... > >& coefficients , Real isoValue , CoredMeshData< Vertex >& mesh , std::function< void ( Vertex& , Point< Real , Dim > , Real , Data ) > SetVertex , bool nonLinearFit , bool addBarycenter , bool polygonMesh , bool flipOrientation )
#endif // NEW_CODE
	{
		IsoStats isoStats;
		static_assert( sizeof...(FEMSigs)==Dim , "[ERROR] Number of signatures should match dimension" );
		tree._setFEM1ValidityFlags( UIntPack< FEMSigs ... >() );
		static const unsigned int DataDegree = FEMSignature< DataSig >::Degree;
		static const int FEMDegrees[] = { FEMSignature< FEMSigs >::Degree ... };
		for( int d=0 ; d<Dim ; d++ ) if( FEMDegrees[d]==0 && nonLinearFit ) WARN( "Constant B-Splines do not support non-linear interpolation" ) , nonLinearFit = false;

		SliceData::SetHyperCubeTables();

		typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >* pointEvaluator = NULL;
		if( data ) pointEvaluator = new typename FEMIntegrator::template PointEvaluator< IsotropicUIntPack< Dim , DataSig > , ZeroUIntPack< Dim > >( tree._maxDepth );
		DenseNodeData< Real , UIntPack< FEMSigs ... > > coarseCoefficients( tree._sNodesEnd( tree._maxDepth-1 ) );
		memset( coarseCoefficients() , 0 , sizeof(Real)*tree._sNodesEnd( tree._maxDepth-1 ) );
#ifdef NEW_THREADS
		tp.parallel_for( tree._sNodesBegin(0) , tree._sNodesEnd( tree._maxDepth-1 ) , [&]( unsigned int , size_t i ){ coarseCoefficients[i] = coefficients[i]; } );
#else // !NEW_THREADS
#pragma omp parallel for
#ifdef NEW_CODE
		for( node_index_type i=tree._sNodesBegin(0) ; i<tree._sNodesEnd( tree._maxDepth-1 ) ; i++ ) coarseCoefficients[i] = coefficients[i];
#else // !NEW_CODE
		for( int i=tree._sNodesBegin(0) ; i<tree._sNodesEnd( tree._maxDepth-1 ) ; i++ ) coarseCoefficients[i] = coefficients[i];
#endif // NEW_CODE
#endif // NEW_THREADS
		typename FEMIntegrator::template RestrictionProlongation< UIntPack< FEMSigs ... > > rp;
#ifdef NEW_THREADS
		for( LocalDepth d=1 ; d<tree._maxDepth ; d++ ) tree._upSample( UIntPack< FEMSigs ... >() , tp , rp , d , coarseCoefficients() );
#else // !NEW_THREADS
		for( LocalDepth d=1 ; d<tree._maxDepth ; d++ ) tree._upSample( UIntPack< FEMSigs ... >() , rp , d , coarseCoefficients() );
#endif // NEW_THREADS
		FEMTree< Dim , Real >::MemoryUsage();

		std::vector< _Evaluator< UIntPack< FEMSigs ... > , 1 > > evaluators( tree._maxDepth+1 );
		for( LocalDepth d=0 ; d<=tree._maxDepth ; d++ ) evaluators[d].set( tree._maxDepth );

#ifdef NEW_CODE
		node_index_type vertexOffset = 0;
#else // !NEW_CODE
		int vertexOffset = 0;
#endif // NEW_CODE

#ifdef NEW_THREADS
		std::vector< _SlabValues > slabValues;
		slabValues.reserve( tree._maxDepth+1 );
		for( unsigned int d=0 ; d<(unsigned int)tree._maxDepth+1 ;  d++ ) slabValues.emplace_back( tp );
#else // !NEW_THREADS
		std::vector< _SlabValues > slabValues( tree._maxDepth+1 );
#endif //NEW_THREADS

		// Initialize the back slice
		for( LocalDepth d=tree._maxDepth ; d>=0 ; d-- )
		{
			double t = Time();
#ifdef NEW_THREADS
			SliceData::SetSliceTableData( tp , tree._sNodes , &slabValues[d].sliceValues(0).sliceData , &slabValues[d].xSliceValues(0).xSliceData , &slabValues[d].sliceValues(1).sliceData , tree._localToGlobal( d ) , tree._localInset( d ) );
#else // !NEW_THREADS
			SliceData::SetSliceTableData( tree._sNodes , &slabValues[d].sliceValues(0).sliceData , &slabValues[d].xSliceValues(0).xSliceData , &slabValues[d].sliceValues(1).sliceData , tree._localToGlobal( d ) , tree._localInset( d ) );
#endif // NEW_THREADS
			isoStats.setTableTime += Time()-t;
			slabValues[d].sliceValues (0).reset( nonLinearFit );
			slabValues[d].sliceValues (1).reset( nonLinearFit );
			slabValues[d].xSliceValues(0).reset( );
		}
		for( LocalDepth d=tree._maxDepth ; d>=0 ; d-- )
		{
			// Copy edges from finer
			double t = Time();
#ifdef NEW_THREADS
			if( d<tree._maxDepth ) _CopyFinerSliceIsoEdgeKeys( tp , tree , d , 0 , slabValues );
#else // !NEW_THREADS
			if( d<tree._maxDepth ) _CopyFinerSliceIsoEdgeKeys( tree , d , 0 , slabValues );
#endif // NEW_THREADS
			isoStats.copyFinerTime += Time()-t , t = Time();
#ifdef NEW_THREADS
			_SetSliceIsoCorners< FEMSigs ... >( tp , tree , coefficients() , coarseCoefficients() , isoValue , d , 0 , slabValues , evaluators[d] );
#else // !NEW_THREADS
			_SetSliceIsoCorners< FEMSigs ... >( tree , coefficients() , coarseCoefficients() , isoValue , d , 0 , slabValues , evaluators[d] );
#endif // NEW_THREADS
			isoStats.cornersTime += Time()-t , t = Time();
#ifdef NEW_THREADS
			_SetSliceIsoVertices< WeightDegree , Data , DataSig >( tp , tree , pointEvaluator , densityWeights , data , isoValue , d , 0 , vertexOffset , mesh , slabValues , SetVertex );
#else // !NEW_THREADS
			_SetSliceIsoVertices< WeightDegree , Data , DataSig >( tree , pointEvaluator , densityWeights , data , isoValue , d , 0 , vertexOffset , mesh , slabValues , SetVertex );
#endif // NEW_THREADS
			isoStats.verticesTime += Time()-t , t = Time();
#ifdef NEW_THREADS
			_SetSliceIsoEdges( tp , tree , d , 0 , slabValues );
#else // !NEW_THREADS
			_SetSliceIsoEdges( tree , d , 0 , slabValues );
#endif // NEW_THREADS
			isoStats.edgesTime += Time()-t , t = Time();
		}

		// Iterate over the slices at the finest level
		for( int slice=0 ; slice<( 1<<tree._maxDepth ) ; slice++ )
		{
			// Process at all depths that contain this slice
			LocalDepth d ; int o;
			for( d=tree._maxDepth , o=slice+1 ; d>=0 ; d-- , o>>=1 )
			{
				// Copy edges from finer (required to ensure we correctly track edge cancellations)
				double t = Time();
				if( d<tree._maxDepth )
				{
#ifdef NEW_THREADS
					_CopyFinerSliceIsoEdgeKeys( tp , tree , d , o , slabValues );
					_CopyFinerXSliceIsoEdgeKeys( tp , tree , d , o-1 , slabValues );
#else // !NEW_THREADS
					_CopyFinerSliceIsoEdgeKeys( tree , d , o , slabValues );
					_CopyFinerXSliceIsoEdgeKeys( tree , d , o-1 , slabValues );
#endif // NEW_THREADS
				}
				isoStats.copyFinerTime += Time()-t , t = Time();
				// Set the slice values/vertices
#ifdef NEW_THREADS
				_SetSliceIsoCorners< FEMSigs ... >( tp , tree , coefficients() , coarseCoefficients() , isoValue , d , o , slabValues , evaluators[d] );
#else // !NEW_THREADS
				_SetSliceIsoCorners< FEMSigs ... >( tree , coefficients() , coarseCoefficients() , isoValue , d , o , slabValues , evaluators[d] );
#endif // NEW_THREADS
				isoStats.cornersTime += Time()-t , t = Time();
#ifdef NEW_THREADS
				_SetSliceIsoVertices< WeightDegree , Data , DataSig >( tp , tree , pointEvaluator , densityWeights , data , isoValue , d , o , vertexOffset , mesh , slabValues , SetVertex );
#else // !NEW_THREADS
				_SetSliceIsoVertices< WeightDegree , Data , DataSig >( tree , pointEvaluator , densityWeights , data , isoValue , d , o , vertexOffset , mesh , slabValues , SetVertex );
#endif // NEW_THREADS
				isoStats.verticesTime += Time()-t , t = Time();
#ifdef NEW_THREADS
				_SetSliceIsoEdges( tp , tree , d , o , slabValues );
#else // !NEW_THREADS
				_SetSliceIsoEdges( tree , d , o , slabValues );
#endif // NEW_THREADS
				isoStats.edgesTime += Time()-t , t = Time();

				// Set the cross-slice edges
#ifdef NEW_THREADS
				_SetXSliceIsoVertices< WeightDegree , Data , DataSig >( tp , tree , pointEvaluator , densityWeights , data , isoValue , d , o-1 , vertexOffset , mesh , slabValues , SetVertex );
#else // !NEW_THREADS
				_SetXSliceIsoVertices< WeightDegree , Data , DataSig >( tree , pointEvaluator , densityWeights , data , isoValue , d , o-1 , vertexOffset , mesh , slabValues , SetVertex );
#endif // NEW_THREADS
				isoStats.verticesTime += Time()-t , t = Time();
#ifdef NEW_THREADS
				_SetXSliceIsoEdges( tp , tree , d , o-1 , slabValues );
#else // !NEW_THREADS
				_SetXSliceIsoEdges( tree , d , o-1 , slabValues );
#endif // NEW_THREADS
				isoStats.edgesTime += Time()-t , t = Time();

#ifdef NEW_THREADS
				auto SlabSet = [&]( unsigned int , size_t i )
				{
					switch( i )
					{
					case 0: slabValues[d]. sliceValues(o-1).setEdgeVertexMap() ; return;
					case 1: slabValues[d]. sliceValues(o  ).setEdgeVertexMap() ; return;
					case 2: slabValues[d].xSliceValues(o-1).setEdgeVertexMap() ; return;
					case 3: slabValues[d]. sliceValues(o-1).setVertexPairMap() ; return;
					case 4: slabValues[d]. sliceValues(o  ).setVertexPairMap() ; return;
					case 5: slabValues[d].xSliceValues(o-1).setVertexPairMap() ; return;
					case 6: slabValues[d]. sliceValues(o-1).setFaceEdgeMap() ; return;
					case 7: slabValues[d]. sliceValues(o  ).setFaceEdgeMap() ; return;
					case 8: slabValues[d].xSliceValues(o-1).setFaceEdgeMap() ; return;
					}
				};
				tp.parallel_for( 0 , 9 , SlabSet , 1 );
#else // !NEW_THREADS
#pragma omp parallel sections
				{
#pragma omp section
					slabValues[d]. sliceValues(o-1).setEdgeVertexMap();
#pragma omp section
					slabValues[d]. sliceValues(o  ).setEdgeVertexMap();
#pragma omp section
					slabValues[d].xSliceValues(o-1).setEdgeVertexMap();
#pragma omp section
					slabValues[d]. sliceValues(o-1).setVertexPairMap();
#pragma omp section
					slabValues[d]. sliceValues(o  ).setVertexPairMap();
#pragma omp section
					slabValues[d].xSliceValues(o-1).setVertexPairMap();
#pragma omp section
					slabValues[d]. sliceValues(o-1).setFaceEdgeMap();
#pragma omp section
					slabValues[d]. sliceValues(o  ).setFaceEdgeMap();
#pragma omp section
					slabValues[d].xSliceValues(o-1).setFaceEdgeMap();
				}
#endif // NEW_THREADS
				// Add the triangles
				t = Time();
#ifdef NEW_THREADS
				_SetIsoSurface( tp , tree , d , o-1 , slabValues[d].sliceValues(o-1) , slabValues[d].sliceValues(o) , slabValues[d].xSliceValues(o-1) , mesh , polygonMesh , addBarycenter , vertexOffset , flipOrientation );
#else // !NEW_THREADS
				_SetIsoSurface( tree , d , o-1 , slabValues[d].sliceValues(o-1) , slabValues[d].sliceValues(o) , slabValues[d].xSliceValues(o-1) , mesh , polygonMesh , addBarycenter , vertexOffset , flipOrientation );
#endif // NEW_THREADS
				isoStats.surfaceTime += Time()-t;

				if( o&1 ) break;
			}

			for( d=tree._maxDepth , o=slice+1 ; d>=0 ; d-- , o>>=1 )
			{
				// Initialize for the next pass
				if( o<(1<<(d+1)) )
				{
					double t = Time();
#ifdef NEW_THREADS
					SliceData::SetSliceTableData( tp , tree._sNodes , NULL , &slabValues[d].xSliceValues(o).xSliceData , &slabValues[d].sliceValues(o+1).sliceData , tree._localToGlobal( d ) , o + tree._localInset( d ) );
#else // !NEW_THREADS
					SliceData::SetSliceTableData( tree._sNodes , NULL , &slabValues[d].xSliceValues(o).xSliceData , &slabValues[d].sliceValues(o+1).sliceData , tree._localToGlobal( d ) , o + tree._localInset( d ) );
#endif // NEW_THREADS
					isoStats.setTableTime += Time()-t;
					slabValues[d].sliceValues(o+1).reset( nonLinearFit );
					slabValues[d].xSliceValues(o).reset();
				}
				if( o&1 ) break;
			}
		}
		FEMTree< Dim , Real >::MemoryUsage();
		if( pointEvaluator ) delete pointEvaluator;
		return isoStats;
	}
};
#ifdef NEW_THREADS
template< class Real , class Vertex > std::mutex IsoSurfaceExtractor< 3 , Real , Vertex >::_pointInsertionMutex;
#endif // NEW_THREADS

template< class Real , class Vertex > template< unsigned int D , unsigned int K >
unsigned int IsoSurfaceExtractor< 3 , Real , Vertex >::SliceData::HyperCubeTables< D , K >::CellOffset[ HyperCube::Cube< D >::template ElementNum< K >() ][ HyperCube::Cube< D >::template IncidentCubeNum< K >() ];
template< class Real , class Vertex > template< unsigned int D , unsigned int K >
unsigned int IsoSurfaceExtractor< 3 , Real , Vertex >::SliceData::HyperCubeTables< D , K >::IncidentElementCoIndex[ HyperCube::Cube< D >::template ElementNum< K >() ][ HyperCube::Cube< D >::template IncidentCubeNum< K >() ];
template< class Real , class Vertex > template< unsigned int D , unsigned int K >
unsigned int IsoSurfaceExtractor< 3 , Real , Vertex >::SliceData::HyperCubeTables< D , K >::CellOffsetAntipodal[ HyperCube::Cube< D >::template ElementNum< K >() ];
template< class Real , class Vertex > template< unsigned int D , unsigned int K >
typename HyperCube::Cube< D >::template IncidentCubeIndex < K > IsoSurfaceExtractor< 3 , Real , Vertex >::SliceData::HyperCubeTables< D , K >::IncidentCube[ HyperCube::Cube< D >::template ElementNum< K >() ];
template< class Real , class Vertex > template< unsigned int D , unsigned int K >
typename HyperCube::Direction IsoSurfaceExtractor< 3 , Real , Vertex >::SliceData::HyperCubeTables< D , K >::Directions[ HyperCube::Cube< D >::template ElementNum< K >() ][ D ];
template< class Real , class Vertex > template< unsigned int D , unsigned int K1 , unsigned int K2 >
typename HyperCube::Cube< D >::template Element< K2 > IsoSurfaceExtractor< 3 , Real , Vertex >::SliceData::HyperCubeTables< D , K1 , K2 >::OverlapElements[ HyperCube::Cube< D >::template ElementNum< K1 >() ][ HyperCube::Cube< D >::template OverlapElementNum< K1 , K2 >() ];
template< class Real , class Vertex > template< unsigned int D , unsigned int K1 , unsigned int K2 >
bool IsoSurfaceExtractor< 3 , Real , Vertex >::SliceData::HyperCubeTables< D , K1 , K2 >::Overlap[ HyperCube::Cube< D >::template ElementNum< K1 >() ][ HyperCube::Cube< D >::template ElementNum< K2 >() ];
