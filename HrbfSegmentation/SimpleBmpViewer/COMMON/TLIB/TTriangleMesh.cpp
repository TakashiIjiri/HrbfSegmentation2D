//////////////////////////////////////////////////////////////////////
// TPolygonMesh.cpp: TPolygonMesh クラスのインプリメンテーション
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "TTriangleMesh.h"
#include <algorithm>
#include <map>

//#include <fstream>
//#include <cfloat>
//#include <cmath>
//#include <deque>

TTriangleMesh::~TTriangleMesh(){ Clear();}
TTriangleMesh::TTriangleMesh(){ 
	m_verts   = 0;	 m_polys  = 0;  m_edges  = 0;
	m_v_norm  = 0;   m_p_edges= 0;
	m_v_edge  = 0;   m_p_norm = 0;
	m_v_uv    = 0;
	m_vSize = m_eSize = m_pSize = 0;
	m_WingEdgeIsAvailable = false;
}
TTriangleMesh::TTriangleMesh(const TTriangleMesh &src)
{
	m_verts   = 0;	 m_polys  = 0;  m_edges  = 0;
	m_v_norm  = 0;   m_p_edges= 0;
	m_v_edge  = 0;   m_p_norm = 0;
	m_v_uv    = 0;
	m_vSize = m_eSize = m_pSize = 0;
	init( src );
}


void TTriangleMesh::Clear()
{
	if( m_verts   != 0) delete[] m_verts  ; m_verts   = 0;
	if( m_polys   != 0) delete[] m_polys  ; m_polys   = 0;
	if( m_edges   != 0) delete[] m_edges  ; m_edges   = 0;
	
	if( m_v_norm  != 0) delete[] m_v_norm ; m_v_norm  = 0;
	if( m_v_edge  != 0) delete[] m_v_edge ; m_v_edge  = 0;
	if( m_v_uv    != 0) delete[] m_v_uv   ; m_v_uv    = 0;

	if( m_p_edges != 0) delete[] m_p_edges; m_p_edges = 0;
	if( m_p_norm  != 0) delete[] m_p_norm ; m_p_norm  = 0;
	m_v_1ring.clear();
	m_vSize = m_eSize = m_pSize = 0;
	m_WingEdgeIsAvailable = false;
}

void TTriangleMesh::init(const TTriangleMesh &pm)
{
	Clear();
	m_vSize = pm.m_vSize;
	m_pSize = pm.m_pSize;
	m_eSize = pm.m_eSize;
	if( pm.m_verts  != 0) m_verts  = new TVector3  [m_vSize];  
	if( pm.m_v_norm != 0) m_v_norm = new TVector3  [m_vSize];   
	if( pm.m_v_edge != 0) m_v_edge = new int       [m_vSize];  
	if( pm.m_v_uv   != 0) m_v_uv   = new TVector3  [m_vSize];   

	if( pm.m_polys  != 0) m_polys  = new TTriangle [m_pSize];
	if( pm.m_p_edges!= 0) m_p_edges= new TTriangle [m_pSize];
	if( pm.m_p_norm != 0) m_p_norm = new TVector3  [m_pSize];
	
	if( pm.m_edges  != 0) m_edges  = new TWingEdge [m_eSize];

	for( int i=0; i<m_vSize; ++i){ m_verts  [i] = pm.m_verts  [i];
	                               m_v_norm [i] = pm.m_v_norm [i];
								   m_v_edge [i] = pm.m_v_edge [i];
								   m_v_uv   [i] = pm.m_v_uv   [i];}

	for( int i=0; i<m_pSize; ++i){ m_polys  [i] = pm.m_polys  [i];
	                               m_p_edges[i] = pm.m_p_edges[i];
								   m_p_norm [i] = pm.m_p_norm [i];}
	for( int i=0; i<m_eSize; ++i)  m_edges  [i] = pm.m_edges  [i];

	m_v_1ring             = pm.m_v_1ring            ;
	m_WingEdgeIsAvailable = pm.m_WingEdgeIsAvailable;

}

bool TTriangleMesh::initFromVertsPolys(const vector<TVector3> &verts, const vector<TTriangle > &polys){
	Clear();
	m_vSize = (int)verts.size();   
	m_verts  = new TVector3[m_vSize]; 
	m_v_norm = new TVector3[m_vSize];
	m_v_uv   = new TVector3[m_vSize];
	
	m_pSize = (int)polys.size(); 
	m_polys  = new TTriangle[m_pSize]; 
	m_p_norm = new TVector3 [m_pSize];


	for( int i=0; i<m_vSize; ++i) m_verts[i].Set( verts[i] );
	for( int i=0; i<m_pSize; ++i) m_polys[i].Set( polys[i] );
	updateNormal  ();
	updateWingEdge();
	m_WingEdgeIsAvailable = false ;

	if( checkVtxMoreThan4BoundEdge() ) {
		fprintf( stderr, "this mesh has one or more viertices that connected to more than 2 boundary edges\n"); 
	}else if( !checkWingedEdge() ){
		fprintf( stderr, "this mesh has strange polygon structure \n"); 
	}
	else {
		m_WingEdgeIsAvailable = true ;
		updateOneRing ();
	}
	return m_WingEdgeIsAvailable;
}





//see http://www.hiramine.com/programming/3dmodelfileformat/objfileformat.html
bool TTriangleMesh::initFromObjFile( const char* fname )
{
	Clear();
	FILE* fp = fopen(fname,"r") ; if( !fp ) return false;

	list<TVector3>  vs_list, uvs_list;
	list<TTriangle >  polys_list, uvpolys_list;
	bool uContainMinus = false, 
		 vContainMinus = false;

	char buf[256] ;		
	while(fgets(buf,255,fp)){//一行読む

		char* bkup = _strdup(buf) ;        //文字列をコピー
		char* token = strtok( buf, " \t" );//最初のトークンを取得

		if( stricmp( token,"vt" ) == 0 ){ // Texture coordinate
			TVector3 vt ;
			sscanf( bkup,"vt %lf %lf",&vt.data[0],&vt.data[1] ) ;
			if(vt[0] < 0) uContainMinus = true;
			if(vt[1] < 0) vContainMinus = true;
			uvs_list.push_back( vt ) ;
		} 
		else if( stricmp( token,"v" ) == 0 ){ // Vertex location
			TVector3 v;
			sscanf( bkup,"v %lf %lf %lf",&v.data[0],&v.data[1],&v.data[2] ) ;
			vs_list.push_back( v ) ;
		} 
		else if( stricmp( token,"f" ) == 0 ){ // Polygon 
			TTriangle p, pUV;
			int tmp;
			int vtxnum = sscanf( bkup,"f %d %d %d %d", &p.idx[0], &p.idx[1], &p.idx[2], &tmp) ;//sscanfの返り値は正常に読めたフィールド数 (/が入ったら2フィールドしか読めない)
			
			if( vtxnum < 3 ) vtxnum = sscanf( bkup,"f %d/%d %d/%d %d/%d %d/%d" ,            &p.idx[0], &pUV.idx[0],
																			                &p.idx[1], &pUV.idx[1],
																			                &p.idx[2], &pUV.idx[2], &tmp, &tmp)/2 ;
			if( vtxnum < 3 ) vtxnum = sscanf( bkup,"f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", &p.idx[0], &pUV.idx[0], &tmp,
																						    &p.idx[1], &pUV.idx[1] ,&tmp,
																						    &p.idx[2], &pUV.idx[2] ,&tmp, &tmp, &tmp, &tmp)/3 ;
			if( vtxnum < 3 ) vtxnum = sscanf( bkup,"f %d//%d %d//%d %d//%d %d//%d" ,        &p.idx[0], &tmp,  
																						    &p.idx[1], &tmp,
																						    &p.idx[2], &tmp,  &tmp, &tmp )/2 ;
			--p  .idx[0]; --p  .idx[1]; --p  .idx[2];
			--pUV.idx[0]; --pUV.idx[1]; --pUV.idx[2];

			polys_list  .push_back( p ) ;
			uvpolys_list.push_back(pUV) ;
		}
		free(bkup) ;
	}
	fclose(fp) ;

	vector< TVector3 > verts; t_copy<TVector3 >( vs_list   , verts ) ;
	vector< TTriangle> polys; t_copy<TTriangle>( polys_list, polys ) ;
	return initFromVertsPolys( verts, polys );
}

bool TTriangleMesh::SaveToObjFile(const char* fname)
{
	fprintf( stderr, "export obj file %s\n\n", fname);
	FILE* fp = fopen(fname,"w") ;
	if( !fp )	return false ;
	fprintf(fp,"#Polygon model by ILPolygonModel\n") ;
	for( int i=0;i<m_vSize; ++i ) fprintf(fp,"v %f %f %f\n",m_verts[i].data[0]  , m_verts[i].data[1]  , m_verts[i].data[2]  ) ;
	for( int i=0;i<m_pSize; ++i ) fprintf(fp,"f %d %d %d\n",m_polys[i].idx [0]+1, m_polys[i].idx [1]+1, m_polys[i].idx [2]+1) ;
	fclose(fp) ;
	return true ;
}
//m_v_norm and m_p_norm are already allocated 
void TTriangleMesh::updateNormal  ()
{
	//if( m_vSize == 0 || m_pSize == 0 ) return;
	for(int i=0; i<m_vSize; ++i) m_v_norm[i].Set(0,0,0);
	for(int i=0; i<m_pSize; ++i)
	{
		int *idx = m_polys[i].idx;
		t_V1subtV2_cros_V3subtV4( m_verts[ idx[1] ], m_verts[ idx[0] ], 
			                      m_verts[ idx[2] ], m_verts[ idx[0] ], m_p_norm[i] );
		if( m_p_norm[i].Normalize_Self() ){
			m_v_norm[ idx[0] ]  += m_p_norm[i];
			m_v_norm[ idx[1] ]  += m_p_norm[i];
			m_v_norm[ idx[2] ]  += m_p_norm[i];
		}
	}
	for(int i=0; i<m_vSize; ++i) m_v_norm[i].Normalize_Self();
}

//Winged Edge Construction///////////////////////////////////////////////////////////////////////////
// consider boundaries// error may occur when the vertex has two boudary(more than three boundary edges
// this construct "m_edges / m_p_edges / m_v_edge from "m_verts / m_polys" 
//If more than three triangle share the same edge (v0-v1), this method return strange results (This is the fault of the input)
static const int edg_mat[3][2]            = {{0,1},{1,2},{2,0}} ;
static const int edgn_prev_next_mat[3][2] = {{2,1},{0,2},{1,0}} ;
void TTriangleMesh::updateWingEdge()
{
	if( m_p_edges != 0 ) delete[] m_p_edges; m_p_edges = new TTriangle[m_pSize];
	if( m_v_edge  != 0 ) delete[] m_v_edge ; m_v_edge  = new int      [m_vSize]; 
	
	vector<TWingEdge> Es; //m_e_*は まだサイズが不明なので new出来ない
	
	for( int i=0;i<m_vSize; ++i) m_v_edge[i] = -1;

	vector< list<int> > emanatingEdges( m_vSize ) ;//list for Eminating edge IDs

	for( int polygon_id = 0; polygon_id <m_pSize; ++polygon_id )
	{
		int *pVtx = m_polys  [polygon_id].idx ;
		int *eIdx = m_p_edges[polygon_id].idx ;
		
		//edgeを生成: v[0,1], p[0,1], oneEdge, polygon.edge登録
		for( int i = 0; i < 3; i++ )
		{
			TWingEdge _we ;
			_we.v[0] = pVtx[ edg_mat[i][0] ] ;
			_we.v[1] = pVtx[ edg_mat[i][1] ] ;
 
			bool bInverted = false ;
			if( _we.v[0] > _we.v[1] ){ bInverted = true ; std::swap( _we.v[0], _we.v[1] );}
			
			list<int>& emanEs = emanatingEdges[ _we.v[0] ] ;
			
			list<int>::iterator we_it;
			for( we_it = emanEs.begin() ; we_it != emanEs.end() ; we_it++ ) if( Es[*we_it].v[1] == _we.v[1] ) break ;
			
			if( we_it == emanEs.end() ){
				Es.push_back( _we ) ;
				emanatingEdges[ _we.v[0] ].push_back( (int)Es.size()-1 ) ;
				emanatingEdges[ _we.v[1] ].push_back( (int)Es.size()-1 ) ;
			}
			
			int edg_id = (we_it == emanEs.end() ? (int)Es.size()-1 : *we_it ) ;
			
			eIdx[i] = m_v_edge[ Es[edg_id].v[0] ] = m_v_edge[ Es[edg_id].v[1] ] = edg_id ;
			if( bInverted ) Es[edg_id].p[1] = polygon_id ;
			else		    Es[edg_id].p[0] = polygon_id ;
		}
	}
	//set e[0,1,2,3]
	for( int polygon_id = 0; polygon_id <m_pSize; ++polygon_id )
	{
		int *pVtx = m_polys  [polygon_id].idx ;// polygonの vertex の [0], [1], [2]
		int *eIdx = m_p_edges[polygon_id].idx ;// polygonの edge   の [0], [1], [2] (これは v0-v1, v1-v2, v2-v0に対応)

		for( int i = 0; i < 3; i++ )
		{
			TWingEdge& we = Es[ eIdx[i] ] ;
			if( we.v[0] == pVtx[ edg_mat[i][0]] ){	// Not inverted
				we.e[0] = eIdx[ edgn_prev_next_mat[i][0] ] ;
				we.e[1] = eIdx[ edgn_prev_next_mat[i][1] ] ;
			} else {
				we.e[3] = eIdx[ edgn_prev_next_mat[i][0] ] ;
				we.e[2] = eIdx[ edgn_prev_next_mat[i][1] ] ;
			}
		}
	}
	//copy
	m_eSize = (int)Es.size();
	if( m_edges != 0) delete[] m_edges; m_edges = new TWingEdge[ Es.size() ];
	for( int i=0, s=(int) Es.size();i<s; ++i) m_edges[i].Set( Es[i] );


}

//update 1 ring neighborhood////////////////////////////
void TTriangleMesh::updateOneRing ()
{
	m_v_1ring.clear();  m_v_1ring.resize( m_vSize );
	for( int i=0; i < m_vSize; ++i) GetVsAroundVertex( i, m_v_1ring[i] );
}

void TTriangleMesh::drawPolygons() const
{
	glBegin(GL_TRIANGLES );
	for( int i = 0; i < m_pSize; ++i){
		int *pIds = m_polys[i].idx;
		glNormal3dv( m_v_norm[ pIds[0] ].data ); glVertex3dv( m_verts[ pIds[0] ].data );
		glNormal3dv( m_v_norm[ pIds[1] ].data ); glVertex3dv( m_verts[ pIds[1] ].data );
		glNormal3dv( m_v_norm[ pIds[2] ].data ); glVertex3dv( m_verts[ pIds[2] ].data );
	}
	glEnd();
}

void TTriangleMesh::drawPolygons_flatShading() const
{
	glBegin(GL_TRIANGLES );
	for( int i = 0; i < m_pSize; ++i){
		
		int *pIds = m_polys[i].idx;
		glNormal3dv( m_p_norm[ i ].data ); 
		glVertex3dv( m_verts[ pIds[0] ].data );
		glVertex3dv( m_verts[ pIds[1] ].data );
		glVertex3dv( m_verts[ pIds[2] ].data );
	}
	glEnd();
}


void TTriangleMesh::drawEdges()
{
	glBegin( GL_LINES );
	for( int i = 0; i < m_eSize; ++i){
		glVertex3dv( m_verts[ m_edges[i].v[0] ].data );
		glVertex3dv( m_verts[ m_edges[i].v[1] ].data );
	}
	glEnd();
}
void TTriangleMesh::drawBoundEdges()
{
	glBegin( GL_LINES );
	for( int i = 0; i < m_eSize; ++i)if( m_edges[i].isBound() ){
		glVertex3dv( m_verts[ m_edges[i].v[0] ].data );
		glVertex3dv( m_verts[ m_edges[i].v[1] ].data );
	}
	glEnd();
}

void TTriangleMesh::drawBoundPoints()
{
	fprintf( stderr, "#");
	glBegin( GL_POINTS );
	for( int i = 0; i < m_eSize; ++i)if( m_edges[i].isBound() )
	{
		glVertex3dv( m_verts[ m_edges[i].v[0] ].data );
		glVertex3dv( m_verts[ m_edges[i].v[1] ].data );
	}
	glEnd();
}

void   TTriangleMesh::trace() 
{
	fprintf( stderr, "vSize = %d  eSize = %d  pSize = %d\n", m_vSize, m_eSize, m_pSize );

	for( int i=0;i<m_vSize;++i ){ TVector3 &v = m_verts[i] ; fprintf(stderr, "Vertex%d:(%f,%f,%f)\n",i, v[0],v[1],v[2]) ;}
	fprintf(stderr, "-----\n") ;
	for( int i=0;i<m_eSize;++i ){ TWingEdge  &e = m_edges[i] ; fprintf(stderr, "Edge%d:%d-%d:%d,%d,%d,%d\n",i,e.v[0],e.v[1] ,e.e[0],e.e[1],e.e[2],e.e[3]) ;}
	fprintf(stderr, "-----\n") ;
	for( int i=0;i<m_pSize;++i ){

		fprintf(stderr, "Polygon%d\n\tVertices: ",i) ;
		for( int j=0;j<3;++j ) fprintf(stderr, " ,%d",m_polys[i].idx[j]) ;		
		fprintf(stderr, "\n\tEdges: ") ;
		if( m_eSize!=0) for( int j=0;j<3;++j ) fprintf(stderr, " ,%d",m_p_edges[i].idx[j]) ;
		fprintf(stderr, "\n") ;
	}
}


//check the existance of strange vertex / more than four boundary edges
bool TTriangleMesh::checkVtxMoreThan4BoundEdge(){
	vector< int > v_bEdgeNum( m_vSize, 0);
	for( int i=0; i<m_eSize; ++i) if( m_edges[i].isBound() ) {v_bEdgeNum[ m_edges[i].v[0] ]++;
	                                                          v_bEdgeNum[ m_edges[i].v[1] ]++;}
	for( int i=0; i<m_vSize; ++i) if( v_bEdgeNum[i] >= 3 ) return true;
	return false;

}

//WingedEdge operation////////////////
bool TTriangleMesh::checkWingedEdge()
{
	vector<int> e_numOfReferred( m_eSize );
	for( int p =0; p<m_pSize; ++p){
		e_numOfReferred[ m_p_edges[p].idx[0] ]++;
		e_numOfReferred[ m_p_edges[p].idx[1] ]++;
		e_numOfReferred[ m_p_edges[p].idx[2] ]++;
	}

	for( int ei = 0; ei < m_eSize; ++ei)
	{
		if( e_numOfReferred[ei] > 2 ) {
			fprintf( stderr, "edge[%d] is referred from %d polygons!\n", ei, e_numOfReferred[ei ] );
			return false;
		}
		const TWingEdge &e = m_edges[ ei ];
		const int v0_i     = e.v[0];
		const int v1_i     = e.v[1];

		if( e.p[0] != -1 )
		{
			const int p0_i = e.p[0];  const int* pVid = m_polys  [p0_i].idx;
			                          const int* pEid = m_p_edges[p0_i].idx;
			const int e0_i = e.e[0];  const TWingEdge &e0 = m_edges[e0_i];
			const int e1_i = e.e[1];  const TWingEdge &e1 = m_edges[e1_i];

			//頂点が一致するかcheck!
			int vId_fromE0 = (e0.v[0] != v0_i ) ? e0.v[0] : e0.v[1];
			int vId_fromE1 = (e1.v[0] != v1_i ) ? e1.v[0] : e1.v[1];

			if( vId_fromE0 != vId_fromE1 ) { fprintf( stderr, "wrong #1 edge == %d\n", ei); return false;}

			bool polygonVtx = false;
			if( pVid[0] == v0_i && pVid[1] == v1_i && pVid[2] == vId_fromE0 ) polygonVtx = true;
			if( pVid[1] == v0_i && pVid[2] == v1_i && pVid[0] == vId_fromE0 ) polygonVtx = true;
			if( pVid[2] == v0_i && pVid[0] == v1_i && pVid[1] == vId_fromE0 ) polygonVtx = true;
			if( pVid[0] == v0_i && pVid[2] == v1_i && pVid[1] == vId_fromE0 ) polygonVtx = true;
			if( pVid[1] == v0_i && pVid[0] == v1_i && pVid[2] == vId_fromE0 ) polygonVtx = true;
			if( pVid[2] == v0_i && pVid[1] == v1_i && pVid[0] == vId_fromE0 ) polygonVtx = true;
			if( !polygonVtx ){ fprintf( stderr, "wrong #2 edge == %d\n", ei);  return false; }

			bool polygonEdge = false;
			if( pEid[0] == ei && pEid[1] == e0_i && pEid[2] == e1_i ) polygonEdge = true;
			if( pEid[1] == ei && pEid[2] == e0_i && pEid[0] == e1_i ) polygonEdge = true;
			if( pEid[2] == ei && pEid[0] == e0_i && pEid[1] == e1_i ) polygonEdge = true;
			if( pEid[0] == ei && pEid[2] == e0_i && pEid[1] == e1_i ) polygonEdge = true;
			if( pEid[1] == ei && pEid[0] == e0_i && pEid[2] == e1_i ) polygonEdge = true;
			if( pEid[2] == ei && pEid[1] == e0_i && pEid[0] == e1_i ) polygonEdge = true;
			if( !polygonEdge ){ fprintf( stderr, "wrong #3 edge == %d\n", ei); return false;}  
		}

		if( e.p[1] != -1)
		{
			const int p1_i = e.p[1];   const int     *pVid = m_polys  [p1_i].idx;
			                           const int     *pEid = m_p_edges[p1_i].idx;
			const int e2_i = e.e[2];   const TWingEdge &e2 = m_edges  [e2_i];
			const int e3_i = e.e[3];   const TWingEdge &e3 = m_edges  [e3_i];

			//頂点が一致するかcheck!
			int vId_fromE2 = (e2.v[0] != v0_i ) ? e2.v[0] : e2.v[1];
			int vId_fromE3 = (e3.v[0] != v1_i ) ? e3.v[0] : e3.v[1];

			if( vId_fromE2 != vId_fromE3 ) { fprintf( stderr, "wrong #4 edge == %d\n", ei); return false; }

			bool polygonVtx = false;
			if( pVid[0] == v0_i && pVid[1] == v1_i && pVid[2] == vId_fromE2 ) polygonVtx = true;
			if( pVid[1] == v0_i && pVid[2] == v1_i && pVid[0] == vId_fromE2 ) polygonVtx = true;
			if( pVid[2] == v0_i && pVid[0] == v1_i && pVid[1] == vId_fromE2 ) polygonVtx = true;
			if( pVid[0] == v0_i && pVid[2] == v1_i && pVid[1] == vId_fromE2 ) polygonVtx = true;
			if( pVid[1] == v0_i && pVid[0] == v1_i && pVid[2] == vId_fromE2 ) polygonVtx = true;
			if( pVid[2] == v0_i && pVid[1] == v1_i && pVid[0] == vId_fromE2 ) polygonVtx = true;
			if( !polygonVtx ){ fprintf( stderr, "wrong #2 edge == %d\n", ei);  return false;}

			bool polygonEdge = false;
			if( pEid[0] == ei && pEid[1] == e2_i && pEid[2] == e3_i ) polygonEdge = true;
			if( pEid[1] == ei && pEid[2] == e2_i && pEid[0] == e3_i ) polygonEdge = true;
			if( pEid[2] == ei && pEid[0] == e2_i && pEid[1] == e3_i ) polygonEdge = true;
			if( pEid[0] == ei && pEid[2] == e2_i && pEid[1] == e3_i ) polygonEdge = true;
			if( pEid[1] == ei && pEid[0] == e2_i && pEid[2] == e3_i ) polygonEdge = true;
			if( pEid[2] == ei && pEid[1] == e2_i && pEid[0] == e3_i ) polygonEdge = true;
			if( !polygonEdge ){ fprintf( stderr, "wrong #2 edge == %d\n", ei);  return false;}
		}
	}

	for( int ei = 0; ei < m_eSize; ++ei){
		TWingEdge& we = m_edges[ ei ];
		if( we.p[0]==-1||we.p[1]==-1 ) continue ;
		int *p1vtx = m_polys[ we.p[0] ].idx ;
		int *p2vtx = m_polys[ we.p[1] ].idx ;
		if     ( p1vtx[0] == p2vtx[0] ){ if( p1vtx[1] != p2vtx[2] && p1vtx[2] != p2vtx[1] ){ fprintf( stderr, "wrong11\n");return false;} }
		else if( p1vtx[0] == p2vtx[1] ){ if( p1vtx[1] != p2vtx[0] && p1vtx[2] != p2vtx[2] ){ fprintf( stderr, "wrong11\n");return false;} }
		else if( p1vtx[0] == p2vtx[2] ){ if( p1vtx[1] != p2vtx[1] && p1vtx[2] != p2vtx[0] ){ fprintf( stderr, "wrong11\n");return false;} }
		else if( p1vtx[1] == p2vtx[0] ){ if( p1vtx[2] !=p2vtx[2] ){ fprintf( stderr, "wrong122\n");return false;} }
		else if( p1vtx[1] == p2vtx[1] ){ if( p1vtx[2] !=p2vtx[0] ){ fprintf( stderr, "wrong12\n");return false;} }
		else                           { if( p1vtx[1] !=p2vtx[2] || p1vtx[2] != p2vtx[1] ){ fprintf( stderr, "wrong1123143\n");return false;} }
	}
	return true ;

	fprintf( stderr, "no wing edge error is found\n");
	return true;
}

void TTriangleMesh::GetVsPsEsAroundVertex(const int vID, vector<int> &ringVs, vector<int> &ringPs, vector<int> &ringEs) const
{
	ringVs.clear(); ringPs.clear(); ringEs.clear();
	if( m_eSize == 0 || m_edges == 0 ) return;

	bool hasBoundary = false  ;
	const int& e0    = m_v_edge[ vID ];
	int        ei    = e0;

	if( ei == -1 ) return;

	list<int> vs, ps, es;
	do {
		const TWingEdge& wee = m_edges[ei] ;	
		bool bFore = (wee.v[0] == vID);

		es.push_back( ei                    );
		vs.push_back( wee.v[ bFore ? 1 : 0] ); 
		ei = wee.e[ bFore ? 0:3];
		if( ei == -1){ hasBoundary =true; break;}

		ps.push_back( wee.p[ bFore ? 0 : 1] );
	} while( ei != e0 ) ;

	//逆に検索して行く
	if(hasBoundary)	
	{
		ei = e0;
		while(true)
		{
			const TWingEdge& wee = m_edges[ei];
			bool bFore = (wee.v[0] == vID);

			if( ei != e0 ){
				es.push_back( ei                    );
				vs.push_back( wee.v[ bFore ? 1 : 0] ); 
			}
			ei = wee.e[bFore ? 2:1] ;
			if(ei == -1) break; //boundaryに到達
			else         ps.push_front(wee.p[ bFore ? 1:0] );
		}
	}
	t_copy<int>(vs, ringVs) ;
	t_copy<int>(ps, ringPs) ;
	t_copy<int>(es, ringEs) ;
}

void TTriangleMesh::GetPsAroundVertex(const int vID, vector<int> &Polys) const
{
	Polys.clear();
	if( m_eSize == 0 || m_edges == 0 ) return;
	bool hasBoundary = false  ;
	const int& e0    = m_v_edge[vID];
	int        ei    = e0;

	if( ei == -1 ) return;

	list<int> ps;
	do {
		const TWingEdge& wee = m_edges[ei] ;	
		bool bFore = (wee.v[0] == vID);

		ei = wee.e[ bFore ? 0:3];
		if( ei == -1){ hasBoundary =true; break;}

		ps.push_back( wee.p[ bFore ? 0 : 1] );
	} while( ei != e0 ) ;

	//逆に検索して行く
	if(hasBoundary)	
	{
		ei = e0;
		while(true)
		{
			const TWingEdge& wee = m_edges[ei];
			bool bFore = (wee.v[0] == vID);

			ei = wee.e[bFore ? 2:1] ;
			if(ei == -1) break; //boundaryに到達
			else         ps.push_front(wee.p[ bFore ? 1:0] );
		}
	}
	t_copy<int>(ps, Polys) ;
}

void TTriangleMesh::GetVsAroundVertex(const int vID, vector<int> &Verts) const
{
	Verts.clear();
	if( m_eSize == 0 || m_edges == 0 ) return;
	if( m_eSize == 0 || m_edges == 0 ) return;
	bool hasBoundary = false  ;
	const int& e0    = m_v_edge[vID];
	int        ei    = e0;

	if( ei == -1 ) return;

	list<int> vs;
	do {
		const TWingEdge& wee = m_edges[ei] ;	
		bool bFore = (wee.v[0] == vID);

		vs.push_back( wee.v[ bFore ? 1 : 0] ); 
		ei = wee.e[ bFore ? 0:3];
		if( ei == -1){ hasBoundary =true; break;}

	} while( ei != e0 ) ;

	//逆に検索して行く
	if(hasBoundary)	
	{
		ei = e0;
		while(true)
		{
			const TWingEdge& wee = m_edges[ei];
			bool bFore = (wee.v[0] == vID);

			if( ei != e0 ){
				vs.push_back( wee.v[ bFore ? 1 : 0] ); 
			}
			ei = wee.e[bFore ? 2:1] ;
			if(ei == -1) break; //boundaryに到達
		}
	}
	t_copy<int>(vs, Verts) ;
}


void TTriangleMesh::GetEsAroundVertex(const int vID, vector<int> &Edges) const
{
	Edges.clear();
	if( m_eSize == 0 || m_edges == 0 ) return;
	bool hasBoundary = false  ;
	const int& e0    = m_v_edge[vID];
	int        ei    = e0;

	if( ei == -1 ) return;

	list<int> es;
	do {
		const TWingEdge& wee = m_edges[ei] ;	
		bool bFore = (wee.v[0] == vID);

		es.push_back( ei                    );
		ei = wee.e[ bFore ? 0:3];
		if( ei == -1){ hasBoundary =true; break;}

	} while( ei != e0 ) ;

	//逆に検索して行く
	if(hasBoundary)	
	{
		ei = e0;
		while(true)
		{
			const TWingEdge& wee = m_edges[ei];
			bool bFore = (wee.v[0] == vID);

			if( ei != e0 ){
				es.push_back( ei                    );
			}
			ei = wee.e[bFore ? 2:1] ;
			if(ei == -1) break; //boundaryに到達
		}
	}
	t_copy<int>(es, Edges) ;
}

int  TTriangleMesh::GetDegreeOfVertex(const int vID){
	vector<int> Es; GetEsAroundVertex( vID, Es ); return (int) Es.size();
}

bool TTriangleMesh::bBoundaryVertex  (const int vID)const{
	vector<int> Es;	GetEsAroundVertex( vID, Es);
	for( int i=0; i<(int)Es.size(); ++i) if( m_edges[ Es[i] ].isBound() ) return true;
	return false;
}


bool t_intersect_RayToPolygons(const TVector3 &rayP, const TVector3 &rayD, const TTriangleMesh &mesh,
													   int &polyID, TVector3 &pos, int bothFrontBack)
{
	const TTriangle *polys = mesh.m_polys;
	const TVector3  *verts = mesh.m_verts;

	TMatrix9 M;
	TVector3 normal, b, stk;

	double distMin = DBL_MAX;
	polyID    = -1;
	bool info = false;
	for( int i = 0; i < mesh.getPnum() ; ++i) 
	{
		const int *pVid = mesh.m_polys[i].idx;
		const TVector3 &v0 = verts[pVid[0]];
		const TVector3 &v1 = verts[pVid[1]];
		const TVector3 &v2 = verts[pVid[2]];

		normal.Set_V1subtV2_crsprd_V3subtV4( v1,v0,v2,v0);
		normal.Normalize_Self();

		if     (bothFrontBack == 1 && normal * rayD >= 0) continue;
		else if(bothFrontBack == 2 && normal * rayD <= 0) continue;
		
		M.Set( v1.data[0]-v0.data[0], v2.data[0]-v0.data[0], - rayD.data[0], 
			   v1.data[1]-v0.data[1], v2.data[1]-v0.data[1], - rayD.data[1], 
			   v1.data[2]-v0.data[2], v2.data[2]-v0.data[2], - rayD.data[2] );
		if( !M.getInvertSelf() ) continue;
		b.SetSubtract( rayP, v0);
		t_MatMultVec( M, b, stk);
		double s = stk.data[0];
		double t = stk.data[1];
		double k = stk.data[2];
		if (s < 0 || 1 < s || t < 0 || 1 < t || 1 < s + t) continue;
		if (k < distMin) {
			distMin = k;
			pos.Set_V1_Add_CoefMultV2( rayP, k, rayD);
			polyID = i;
			info = true;
		}
	}
	return info;
}



bool t_intersect_LineSegmToPolygons(const TVector3 &p0, const TVector3 &p1, const TTriangleMesh &mesh, int &polyID, TVector3 &pos, double &distInP0P1)
{
	const TTriangle *polys = mesh.m_polys;
	const TVector3 *verts  = mesh.m_verts;

	TMatrix9 M;
	TVector3 b, stk;

	for( int i = 0; i < mesh.getPnum() ; ++i) 
	{
		const int *pVid = mesh.m_polys[i].idx;
		const TVector3 &v0 = verts[pVid[0]];
		const TVector3 &v1 = verts[pVid[1]];
		const TVector3 &v2 = verts[pVid[2]];

		M.Set( v1[0]-v0[0], v2[0]-v0[0], - (p1[0]-p0[0]), 
			   v1[1]-v0[1], v2[1]-v0[1], - (p1[1]-p0[1]), 
			   v1[2]-v0[2], v2[2]-v0[2], - (p1[2]-p0[2]) );
		if( !M.getInvertSelf() ) continue;

		b.SetSubtract( p0, v0);
		t_MatMultVec( M, b, stk);

		if (stk[0] < 0 || 1 < stk[0] || stk[1] < 0 || 1 < stk[1] || 1 < stk[0] + stk[1] || stk[2] < 0 || 1 < stk[2]) continue;
		distInP0P1 = stk[2];
		pos.Set_V1_Add_CoefMultV2( p0, stk[2], p1-p0);
		polyID = i;
		return true;
	}
	return false;
}

	
bool   TTriangleMesh::pickSurface          (const TVector3 &rayP, const TVector3 &rayD, TVector3 &pos, int &polyID, int bothFrontBack ) const{
	return t_intersect_RayToPolygons( rayP, rayD, *this, polyID, pos, bothFrontBack );
}
bool TTriangleMesh::getCrossingPolygonToRay(const TVector3 &rayP, const TVector3 &rayD, int &polyID, TVector3 &pos, int bothFrontBack) const{
	return t_intersect_RayToPolygons( rayP, rayD, *this, polyID, pos, bothFrontBack );
}
bool TTriangleMesh::getCrossingPolygonToLineSegm(const TVector3 &p0  , const TVector3 &p1  , int &polyID, TVector3 &pos, double &distInP0P1) const{
	return t_intersect_LineSegmToPolygons( p0, p1, *this, polyID, pos, distInP0P1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//Root 3 Subdivision boundaryにも対応 //trgtPs[i] = trueのpolygonのみsubdivision////////////////////
void TTriangleMesh::subdivisionRoot3( const byte *trgtPs )
{
	updateWingEdge();

	int newVertNum = m_vSize;
	for( int i= 0; i< m_pSize; ++i) if( trgtPs[i] ) newVertNum++;

	vector< TVector3 > newVerts( newVertNum );
	vector< TTriangle  > newPolys;

	//vertexを作る//
	for(int i = 0; i < m_vSize; ++i) newVerts[i] = m_verts[i];

	vector<int> idmap_oldP2newV( m_pSize, -1 );
	int vIdx = m_vSize;
	for(int i = 0; i < m_pSize; ++i) if( trgtPs[i] )
	{
		const int *pVid  = m_polys[i].idx;
		newVerts[ vIdx ].SetGravityCenter( m_verts[ pVid[0] ], m_verts[ pVid[1] ], m_verts[ pVid[2] ]);
		idmap_oldP2newV[i] = vIdx;
		vIdx++;
	}

	//A polygonを作る 1 非分割polygonを挿入し 2 各エッジ上にpolygonを二個ずつ生成していく
	for(int i = 0; i < m_pSize; ++i) if( !trgtPs[i] ) newPolys.push_back( m_polys[i] );

	for(int i = 0; i < m_eSize; ++i)
	{
		TWingEdge &we = m_edges[i];
		if     ( (we.p[0] == -1 || !trgtPs[we.p[0]] ) && (we.p[1] == -1 || !trgtPs[we.p[1]])){ continue;}
		else if( (we.p[0] != -1 &&  trgtPs[we.p[0]] ) && (we.p[1] == -1 || !trgtPs[we.p[1]])) newPolys.push_back( TTriangle( we.v[1], idmap_oldP2newV[we.p[0]], we.v[0] ) );
		else if( (we.p[0] == -1 || !trgtPs[we.p[0]] ) && (we.p[1] != -1 &&  trgtPs[we.p[1]])) newPolys.push_back( TTriangle( we.v[0], idmap_oldP2newV[we.p[1]], we.v[1] ) );
		else if( (we.p[0] != -1 &&  trgtPs[we.p[0]] ) && (we.p[1] != -1 &&  trgtPs[we.p[1]])){newPolys.push_back( TTriangle( we.v[0], idmap_oldP2newV[we.p[1]], idmap_oldP2newV[we.p[0]] ) );
																							  newPolys.push_back( TTriangle( we.v[1], idmap_oldP2newV[we.p[0]], idmap_oldP2newV[we.p[1]] ) );}
	}
	this->initFromVertsPolys(newVerts, newPolys);

	fprintf( stderr, "subdiv %d %d\n", newVerts.size(), newPolys.size());
	fprintf( stderr, "subdiv %d %d %d\n", m_vSize, m_pSize, m_eSize);
}

//boundaryにも対応 //trgtPs[i] = trueのpolygonのみsubdivision
void TTriangleMesh::subdivisionRoot3( const  vector<byte> &trgtPs )
{
	byte *trgt = new byte[ trgtPs.size() ];
	for( int i=0, s=(int)trgtPs.size(); i<s;++i) trgt[i] = trgtPs[i];
	subdivisionRoot3( trgt );
	delete[] trgt;
}


void TTriangleMesh::subdivisionRoot3()
{
	byte *trgt = new byte[ m_pSize ];
	for( int i=0; i<m_pSize; ++i) trgt[i] =1;
	subdivisionRoot3( trgt );
	delete[] trgt;
}




/*

   v0
    | \
    |   \             e0:subdivするedge
  e0|     \
    |       \
    |         \
  v1|___________\v2

*/
static void subDiv_1( const int *vIds, const int *eIds, const TWingEdge *Es, const vector<byte> e_bSubdiv, vector<int > &map_e2v,  vector<TTriangle> &Ps )
{
	int v0,v1,v2;

	const int  eId = (e_bSubdiv[ eIds[0] ] ) ? eIds[0]: 
		             (e_bSubdiv[ eIds[1] ] ) ? eIds[1]: eIds[2];
	const TWingEdge &e = Es[ eId ];

	if     ( e.hasV0V1( vIds[0], vIds[1]) ) { v0=vIds[0]; v1=vIds[1]; v2=vIds[2]; }
	else if( e.hasV0V1( vIds[1], vIds[2]) ) { v0=vIds[1]; v1=vIds[2]; v2=vIds[0]; }
	else                                    { v0=vIds[2]; v1=vIds[0]; v2=vIds[1]; }

	Ps.push_back( TTriangle( v0,map_e2v[eId ],v2 ) );
	Ps.push_back( TTriangle( v1,v2,map_e2v[eId ] ) );
}


/*

   v0
    | \
    |   \             e0:subdivしない
  e0|     \e2　　　   e1:subdivする
    |       \         e2:subdivする
    |         \
  v1|___________\v2
         e1
*/
static void subDiv_2( const int *vIds, const int *eIds, const TWingEdge *Es, const TVector3 *Vs, const vector<byte> e_bSubdiv, vector<int > &map_e2v,  vector<TTriangle> &Ps )
{
	int v0,v1,v2, e0,e1,e2;

	e0  = (!e_bSubdiv[ eIds[0] ] ) ? eIds[0]: 
		  (!e_bSubdiv[ eIds[1] ] ) ? eIds[1]: eIds[2];
	const TWingEdge &e = Es[ e0 ];

	if     ( e.hasV0V1( vIds[0], vIds[1]) ) { v0=vIds[0]; v1=vIds[1]; v2=vIds[2]; }
	else if( e.hasV0V1( vIds[1], vIds[2]) ) { v0=vIds[1]; v1=vIds[2]; v2=vIds[0]; }
	else                                    { v0=vIds[2]; v1=vIds[0]; v2=vIds[1]; }


	//e1, e2のidを探す
	if     (  eIds[0] != e0 && Es[eIds[0]].hasV0V1( v1, v2 ) )   e1 = eIds[0];
	else if(  eIds[1] != e0 && Es[eIds[1]].hasV0V1( v1, v2 ) )   e1 = eIds[1];
	else/*if( eIds[2] != e0 && Es[eIds[2]].hasV0V1( v1, v2 ) )*/ e1 = eIds[2];
	if     (  eIds[0] != e0 && eIds[0] != e1)   e2 = eIds[0];
	else if(  eIds[1] != e0 && eIds[1] != e1)   e2 = eIds[1];
	else/*if( eIds[2] != e0 && eIds[2] != e1) */e2 = eIds[2];


	double d_v0e1 = t_distance_sq( Vs[v0], Vs[ map_e2v[e1] ] );
	double d_v1e2 = t_distance_sq( Vs[v1], Vs[ map_e2v[e2] ] );

	Ps.push_back( TTriangle(v2, map_e2v[e2], map_e2v[e1]) );
	if( d_v0e1 < d_v1e2 ){
		Ps.push_back( TTriangle(v0, v1, map_e2v[e1] ) );
		Ps.push_back( TTriangle(v0, map_e2v[e1], map_e2v[e2]) );
	}else{
		Ps.push_back( TTriangle(v0, v1, map_e2v[e2]) );
		Ps.push_back( TTriangle(v1, map_e2v[e1], map_e2v[e2]) );
	}
}








//新しい頂点は後ろに追加
//面はばらばらになるかも
bool TTriangleMesh::subdivision_divideEdges( const vector<byte> e_bSubdiv, const vector<TVector3> &e_vPos )
{
	vector<TTriangle >  Ps;
	vector<TVector3  >  Vs;

	//まず頂点を作る
	Vs.reserve( m_vSize );
	for( int i=0; i<m_vSize; ++i ) Vs.push_back( m_verts[i] );

	//分割できるedgeを計算
	int DivideEdgeNum = 0;
	for( int i=0; i<m_eSize; ++i) if( e_bSubdiv[i] ) DivideEdgeNum++;

	vector<int > map_e2v( m_eSize, -1);//edge 分割後に生成されるvertexId
	for( int i=0; i<m_eSize; ++i)if( e_bSubdiv[i] )
	{
		map_e2v[i] = (int)Vs.size();
		Vs.push_back( e_vPos[i] );
	}

	//triangle作成
	for( int i=0; i<m_pSize; ++i)
	{
		int num = 0;
		if( e_bSubdiv[ m_p_edges[i].idx[0] ] ) num++;
		if( e_bSubdiv[ m_p_edges[i].idx[1] ] ) num++;
		if( e_bSubdiv[ m_p_edges[i].idx[2] ] ) num++;


		if(      num == 0) Ps.push_back( m_polys[i] );
		else if( num == 1) subDiv_1( m_polys[i].idx, m_p_edges[i].idx, m_edges         , e_bSubdiv, map_e2v, Ps );
		else if( num == 2) subDiv_2( m_polys[i].idx, m_p_edges[i].idx, m_edges, m_verts, e_bSubdiv, map_e2v, Ps );
		else {
			fprintf( stderr, "error , 今はこの場合分けはサポートしない。なぜなら来ないはずだから。このコメントが表示されたら、絶対なんかおかしい");
			return false;
		}
	}

	return initFromVertsPolys( Vs, Ps );
}








void TTriangleMesh::initAsIcosahedron( double r )
{
	double a = 1.0/ sqrt(5.) ; 
	double b = (1 - a) / 2   ; 
	double c = (1 + a) / 2   ; 
	double d = sqrt(b)       ; 
	double e = sqrt(c)       ; 
	vector<TVector3> Vs(12);
	vector<TTriangle > Ps(20);
	Vs[ 0] = TVector3( 0,  1,  0    );
	Vs[ 1] = TVector3( 0,  a,  2 * a);
	Vs[ 2] = TVector3( e,  a,  b    );
	Vs[ 3] = TVector3( d,  a, -c    );
	Vs[ 4] = TVector3(-d,  a, -c    );
	Vs[ 5] = TVector3(-e,  a,  b    );
	Vs[ 6] = TVector3( d, -a,  c    );
	Vs[ 7] = TVector3( e, -a, -b    );
	Vs[ 8] = TVector3( 0, -a, -2 * a);
	Vs[ 9] = TVector3(-e, -a, -b    );
	Vs[10] = TVector3(-d, -a,  c    );
	Vs[11] = TVector3( 0, -1,  0    ); 
	for( int i=0; i<12; ++i) Vs[i] *= r;
	Ps[ 0].Set( 0, 1, 2); Ps[ 1].Set( 0, 2, 3); Ps[ 2].Set( 0, 3, 4); Ps[ 3].Set( 0, 4, 5);
	Ps[ 4].Set( 0, 5, 1); Ps[ 5].Set( 1, 6, 2); Ps[ 6].Set( 2, 7, 3); Ps[ 7].Set( 3, 8, 4);
	Ps[ 8].Set( 4, 9, 5); Ps[ 9].Set( 5,10, 1); Ps[10].Set( 1,10, 6); Ps[11].Set( 2, 6, 7);
	Ps[12].Set( 3, 7, 8); Ps[13].Set( 4, 8, 9); Ps[14].Set( 5, 9,10); Ps[15].Set(11, 6,10);
	Ps[16].Set(11, 7, 6); Ps[17].Set(11, 8, 7); Ps[18].Set(11, 9, 8); Ps[19].Set(11,10, 9);
	initFromVertsPolys( Vs, Ps );
}
void TTriangleMesh::initAsSphere     ( double r, int resH, int resV)
{
	vector<TVector3> Vs( 2 + resH * (resV - 2)            );
	vector<TTriangle > Ps( 2 * resH + (resV - 3) * resH * 2 );
	Vs.front().Set(0, r, 0); 
	Vs.back ().Set(0,-r, 0);
	int vIdx = 1;
	for( int i = 1; i < resV-1; ++i)
	{
		double phai = M_PI / (double)(resV-1) * i ;
		double tmpR = sin( phai );
		double y    = cos( phai );
		for( int j = 0; j < resH  ; ++j, ++vIdx){
			double theta = 2 * M_PI / (double)resH * j;
			Vs[ vIdx ].Set(  r * tmpR * cos( theta ), 
				             r *  y                 , 
				           - r * tmpR * sin( theta ) );
		}
	}
	
	//polygonを作る
	for( int j = 0; j < resH; ++j) Ps[j].Set( 0, 1 + j, 1 + ( (j==resH-1)? 0 : j+1) );

	int pIdx = resH;
	for( int i = 1; i < resV -2 ; ++i)
	for( int j = 0; j < resH    ; ++j)
	{
		int v0 = 1 + (i-1    ) * resH + j                       ;//i,j   v0------v3
		int v1 = 1 + (i-1 + 1) * resH + j                       ;//i,j    |      |
		int v2 = 1 + (i-1 + 1) * resH + ((j==resH-1)? 0 : j + 1);//i,j    |      |
		int v3 = 1 + (i-1    ) * resH + ((j==resH-1)? 0 : j + 1);//i,j   v1------v2
		Ps[pIdx  ].Set( v0, v1, v2);
		Ps[pIdx+1].Set( v0, v2, v3);
		pIdx += 2;
	}
	for( int j = 0; j < resH; ++j)
	{
		Ps[pIdx].Set( (int) Vs.size() -1, 
			          1 + (resV-3) * resH + ( (j==resH-1)? 0 : j+1),
					  1 + (resV-3) * resH + j );
		pIdx++;
	}
	initFromVertsPolys( Vs, Ps );
}
void TTriangleMesh::initAsCube( double sizeX, double sizeY, double sizeZ)
{	
	vector<TVector3> Vs( 8);
	vector<TTriangle > Ps(20);
	Vs[0].Set( 0    ,  0   , 0	   );
	Vs[1].Set( sizeX,  0   , 0	   );
	Vs[2].Set( sizeX, sizeY, 0	   );
	Vs[3].Set( 0    , sizeY, 0     );
	Vs[4].Set( 0    ,  0   , sizeZ );
	Vs[5].Set( sizeX,  0   , sizeZ );
	Vs[6].Set( sizeX, sizeY, sizeZ );
	Vs[7].Set( 0    , sizeY, sizeZ );
	Ps[ 0].Set( 0,2,1) ;
	Ps[ 1].Set( 0,3,2) ;
	Ps[ 2].Set( 4,5,6) ;
	Ps[ 3].Set( 4,6,7) ;
	Ps[ 4].Set( 0,1,5) ;
	Ps[ 5].Set( 0,5,4) ;
	Ps[ 6].Set( 3,6,2) ;
	Ps[ 7].Set( 3,7,6) ;
	Ps[ 8].Set( 1,2,6) ;
	Ps[ 9].Set( 1,6,5) ;
	Ps[10].Set( 0,7,3) ;
	Ps[11].Set( 0,4,7) ;
	initFromVertsPolys( Vs, Ps );
}

double TTriangleMesh::CalcAverageEdgeLength()
{
	if( m_eSize == 0 || m_edges == 0 ) return 0;
	double d = 0;
	for( int i=0; i< m_eSize;++i) d += t_distance( m_verts[ m_edges[i].v[0] ], m_verts[ m_edges[i].v[1] ]);
	return d / (double)m_eSize;

}

TVector3 TTriangleMesh::CalcAverageNormal    ()
{
	TVector3 n;
	for( int i=0; i<m_vSize; ++i) n += m_v_norm[i];
	return n.Normalize();
}


struct comp_sortVertsPolys_byXvalue{
	bool operator()(const pair<double,int>& a, const pair<double,int>&b) const { return b.first > a.first; }
};

void   TTriangleMesh::sortVerts_byXvalue()//ランダムアクセスが多少改善される
{
	vector< TVector3 > Vs( m_vSize );
	vector< TTriangle  > Ps( m_pSize );

	//sort//
	vector< pair<double,int> > sortedVid( m_vSize );
	for( int i=0; i<m_vSize; ++i) sortedVid[i] = make_pair( m_verts[i].data[0], i);
	sort( sortedVid.begin(), sortedVid.end(), comp_sortVertsPolys_byXvalue()); 

	//map //
	vector< int > map_oldV2sortV( m_vSize, -1);
	for( int i=0; i<m_vSize; ++i) {
		map_oldV2sortV [ sortedVid[i].second ] = i;
		Vs[i] = m_verts[ sortedVid[i].second ];
	}

	for( int i=0; i<m_pSize; ++i){
		Ps[i].idx[0] = map_oldV2sortV[ m_polys[i].idx[0] ];
		Ps[i].idx[1] = map_oldV2sortV[ m_polys[i].idx[1] ];
		Ps[i].idx[2] = map_oldV2sortV[ m_polys[i].idx[2] ];
	}
	initFromVertsPolys( Vs, Ps );
}

void   TTriangleMesh::sortPolys_byXvalue()//ランダムアクセスが多少改善され
{
	//sort//
	vector< pair<double,int> > sortedPid( m_pSize );
	for( int i=0; i<m_pSize ; ++i) 
		sortedPid[i] = make_pair( m_verts[ m_polys[i].idx[0] ].data[0] + 
			                      m_verts[ m_polys[i].idx[1] ].data[0] + 
							      m_verts[ m_polys[i].idx[2] ].data[0] , i);
	sort( sortedPid.begin(), sortedPid.end(), comp_sortVertsPolys_byXvalue()); 

	vector< TVector3 > Vs( m_vSize );
	vector< TTriangle  > Ps( m_pSize );
	for( int i=0; i<m_pSize; ++i) Ps[i] = m_polys[ sortedPid[i].second ];
	for( int i=0; i<m_vSize; ++i) Vs[i] = m_verts[i];
	initFromVertsPolys(Vs,Ps);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//SMOOTHING///////////////////////////////////////////////////////////////////////////////////////////////
void TTriangleMesh::SmoothingByLaplacianNormal(int time)
{
	TVector3 norm, gCenter, lapracian;	
	vector<TVector3> newVertPos( m_vSize );

	for( int count = 0; count < time; ++count)
	{
		updateNormal();
		for( int i = 0; i < m_vSize; ++i)
		{
			const vector<int> &Vs = m_v_1ring[i];

			norm   .Set( m_v_norm[i] );

			gCenter.Set(0,0,0);
			for(int k = 0,s=(int)Vs.size(); k < s; ++k) gCenter += m_verts[ Vs[k] ];
			gCenter  /= (double) Vs.size();
			
			lapracian.SetSubtract( m_verts[i], gCenter);
			newVertPos[i].SetAdditionWithCoef( lapracian * norm, norm, 1, gCenter);
		}
		for(int i = 0 ; i < m_vSize; ++i) m_verts[i] = newVertPos[i];
	}
	updateNormal();
}
void TTriangleMesh::SmoothingByLaplacianZero(int time)
{
	TVector3 *newVerts = new TVector3[ m_vSize ];

	for( int count = 0; count < time; ++count)
	{
		for( int i = 0; i < m_vSize; ++i)
		{			
			vector<int> &Vs = m_v_1ring[i];
			newVerts[i].Set(0,0,0);
			for(int k = 0, s = (int)Vs.size(); k < s; ++k) newVerts[i] += m_verts[ Vs[k] ];
			newVerts[i] /= (double) Vs.size();
		}

		for(int i = 0 ; i < m_vSize; ++i) m_verts[i] = newVerts[i];
	}
	updateNormal();
	delete[] newVerts;
}

void TTriangleMesh::SmoothingByLaplacianZero_fixBound(int time)
{
	TVector3 *newVerts = new TVector3[ m_vSize ];

	for( int count = 0; count < time; ++count)
	{
		for( int i = 0; i < m_vSize; ++i) 
		{			
			vector<int> &Vs = m_v_1ring[i];
			newVerts[i].Set(0,0,0);
			for(int k = 0, s = (int)Vs.size(); k < s; ++k) newVerts[i] += m_verts[ Vs[k] ];
			newVerts[i] /= (double) Vs.size();
		}

		for(int i = 0 ; i < m_vSize; ++i) if( !bBoundaryVertex(i) ) m_verts[i] = newVerts[i];
	}
	updateNormal();
	delete[] newVerts;
}

void TTriangleMesh::SmoothingByLaplaceMagDiff(int time)
{
	fprintf( stderr, "smoothing LM diff verts=%d edges = %d polys = %d", m_vSize, m_eSize, m_pSize );

	TVector3 *newVerts = new TVector3[m_vSize];
	double     *LM_curr  = new double    [m_vSize];
	double     *LM_new   = new double    [m_vSize];

	vector<short> v_bBound(m_vSize, false);
	for( int i=0; i<m_eSize; ++i) if( m_edges[i].isBound() ) v_bBound[ m_edges[i].v[0] ] = v_bBound[ m_edges[i].v[1] ] = true;


	TVector3 lap, gc;
	for( int count = 0; count < time; ++count)
	{
		//calc LM
		memset( LM_new , 0, sizeof( double ) * m_vSize);
		memset( LM_curr, 0, sizeof( double ) * m_vSize);
		for( int i = 0; i < m_vSize; ++i) if( !v_bBound[i] )
		{
			gc.Set(0,0,0);
			for(int k = 0, s=(int)m_v_1ring[i].size(); k < s; ++k) gc += m_verts[ m_v_1ring[i][k] ];
			gc /= (double)m_v_1ring[i].size();
			lap.SetSubtract( m_verts[i], gc );
			LM_curr[i] = lap.Length();
			if( lap * m_v_norm[i] < 0 ) LM_curr[0] *= -1;
		}
		//diffuse LM (Laplacian == 0)
		for( int i = 0; i < m_vSize; ++i) if( !v_bBound[i] )
		{
			for(int k = 0, s=(int)m_v_1ring[i].size(); k < s; ++k) LM_new[i] += LM_curr[ m_v_1ring[i][k] ];
			LM_new[i] /= (double)m_v_1ring[i].size();
			LM_new[i] = 0.5 * ( LM_new[i] + LM_curr[i]);
		}

		updateNormal();
		for( int i = 0; i < m_vSize; ++i) 
		{
			if( v_bBound[i] ) { newVerts[i] = m_verts[i]; continue;}
			gc.Set(0,0,0);
			for(int k = 0, s=(int)m_v_1ring[i].size(); k < s; ++k) gc += m_verts[ m_v_1ring[i][k] ];
			gc /= (double)m_v_1ring[i].size();

			newVerts[i] = gc + LM_new[i] * m_v_norm[i];
		}
		for(int i = 0 ; i < m_vSize; ++i) m_verts[i] = newVerts[i];
	}
	updateNormal();

	delete[] newVerts;
	delete[] LM_curr ;
	delete[] LM_new  ;
}


void TTriangleMesh::SmoothingByLaplacianZero_OnlyShapeEdge(int time,  double dotThresh)
{
	fprintf( stderr, "smoothing LM diff verts=%d edges = %d polys = %d", m_vSize, m_eSize, m_pSize );

	TVector3 *newVerts = new TVector3[ m_vSize ];

	for( int count = 0; count < time; ++count)
	{
		for( int i = 0; i < m_vSize; ++i)
		{			
			vector<int> Es, Vs, Ps; GetVsPsEsAroundVertex(i, Vs, Ps, Es);

			bool isSharpEdge = false;
			for( int kk=0; kk<(int)Es.size(); ++kk) if( !m_edges[Es[kk]].isBound() )
			{
				const TVector3 &n0 = m_p_norm[ m_edges[ Es[kk] ].p[0] ];
				const TVector3 &n1 = m_p_norm[ m_edges[ Es[kk] ].p[1] ];
				if( n0*n1 < dotThresh) {isSharpEdge = true; break;}
			}
			if( !isSharpEdge ) { newVerts[i] = m_verts[i]; continue;}

			newVerts[i].Set(0,0,0);
			for(int k = 0; k < (int)Vs.size(); ++k) newVerts[i] += m_verts[ Vs[k] ];
			newVerts[i]  /= (double) Vs.size();

		}
		for(int i = 0 ; i < m_vSize; ++i) m_verts[i] = newVerts[i];
	}
	updateNormal();
	delete[] newVerts;
}

void TTriangleMesh::SmoothingAtOneVertex(int vIdx)
{
	if( vIdx < 0 || m_vSize-1 < vIdx) return;
	vector<int> Vs; GetVsAroundVertex( vIdx, Vs);
	m_verts[vIdx].Set(0,0,0);
	for( int i=0;i<(int)Vs.size();++i) m_verts[vIdx] += m_verts[Vs[i]]; 
	m_verts[vIdx] /= (double)Vs.size();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//Topological operations/////////////////////////////////////////////////////////////////////////////
//m_edges       の末尾に初期化された edgeをnum個追加する //wingEdge should be constructed 
void TTriangleMesh::pushbackNewEdge( int num )
{

	TWingEdge *newEdges = new TWingEdge[m_eSize + num];
	for( int i=0; i<m_eSize; ++i) newEdges[i].Set( m_edges[i] );
	swap( m_edges, newEdges );
	
	m_eSize += num;
	if( newEdges != 0 ) delete[] newEdges;

}
//m_polysとm_p_*の末尾に初期化された edgeをnum個追加する//wingEdge should be constructed 
void TTriangleMesh::pushbackNewPoly( int num )
{
	TTriangle  *newPolys  = new TTriangle [m_pSize + num];
	TTriangle  *newPEdges = new TTriangle [m_pSize + num];
	TVector3 *newPnorm  = new TVector3[m_pSize + num];
	for( int i=0; i< m_pSize; ++i){
		newPolys [i].Set( m_polys  [i] );
		newPEdges[i].Set( m_p_edges[i] );
		newPnorm [i].Set( m_p_norm [i] );
	}
	swap( m_polys  , newPolys );
	swap( m_p_edges, newPEdges);
	swap( m_p_norm , newPnorm );
	if( newPolys  != 0 ) delete[] newPolys;
	if( newPEdges != 0 ) delete[] newPEdges;
	if( newPnorm  != 0 ) delete[] newPnorm;
	m_pSize += num;
}

//m_vertsとm_v_*の末尾に初期化された edgeをnum個追加する
void TTriangleMesh::pushbackNewVert( int num)
{
	TVector3 *newVerts = new TVector3[m_vSize + num];
	TVector3 *newVnorm = new TVector3[m_vSize + num];
	int        *newVEdge = new int       [m_vSize + num];
	for( int i=0; i<m_vSize; ++i)
	{
		newVerts[i].Set( m_verts[i] );
		newVnorm[i].Set( m_v_norm[i] );
		newVEdge[i] = m_v_edge[i];
	}
	for( int i=0; i<num; ++i) newVEdge[m_vSize+i] = -1;
	
	swap( m_verts , newVerts );
	swap( m_v_norm, newVnorm );
	swap( m_v_edge, newVEdge );
	if( newVerts  != 0 ) delete[] newVerts;
	if( newVnorm  != 0 ) delete[] newVnorm;
	if( newVEdge  != 0 ) delete[] newVEdge;
	m_vSize += num;
}


/*----------------------------------------------------------------------------
新しいvertexを配置し、trgtEdgeを分割する
topology周りの処理を行う
push backはすべて最初に行う. 参照のリンク切れが起こらないように注意
---------------------------------------------------------------------------------*/
void TTriangleMesh::splitEdge( int trgtEdgeId, const TVector3 &newVertexPos, int &newVertexId)
{
	if( trgtEdgeId < 0 || m_eSize -1 < trgtEdgeId ) return;

	//add 1vertex / 2 or 3 edges / 1 or 2 polygons///////////////////////////////////
	newVertexId   = m_vSize;  pushbackNewVert(1);
	int newEpivId = m_eSize;  pushbackNewEdge(1);
	m_verts[m_vSize-1] = newVertexPos;

	int newE0_Id = -1, newP0_Id = -1; 
    int newE1_Id = -1, newP1_Id = -1;

	if( m_edges[trgtEdgeId].p[0] != -1 ){
		newE0_Id = m_eSize;  pushbackNewEdge(1);
		newP0_Id = m_pSize;  pushbackNewPoly(1);
	}
	if( m_edges[trgtEdgeId].p[1] != -1 ){
		newE1_Id = m_eSize;  pushbackNewEdge(1);
		newP1_Id = m_pSize;  pushbackNewPoly(1);
	}

	TWingEdge  &e     = m_edges[trgtEdgeId ];
	TVector3 &v0    = m_verts[ e.v[0]    ]; int &v0_oneEdge = m_v_edge[ e.v[0] ];
	TVector3 &v1    = m_verts[ e.v[1]    ]; int &v1_oneEdge = m_v_edge[ e.v[1] ];
	TVector3 &v_n   = m_verts[newVertexId]; int &vn_oneEdge = m_v_edge[newVertexId];
	TWingEdge  &e_nPiv= m_edges[newEpivId  ];

	const int pre_v0 = e.v[0], pre_v1 = e.v[1];
	const int pre_p0 = e.p[0], pre_p1 = e.p[1];
	const int pre_e0 = e.e[0], pre_e1 = e.e[1];
	const int pre_e2 = e.e[2], pre_e3 = e.e[3];

	e     .v[1] = newVertexId;
	e_nPiv.v[0] = pre_v1     ;
	e_nPiv.v[1] = newVertexId;
	v0_oneEdge  = trgtEdgeId;
	v1_oneEdge  = newEpivId;
	vn_oneEdge  = newEpivId;

	if( pre_p0 != -1 )
	{
		TWingEdge &e_n0 = m_edges[newE0_Id];
		TWingEdge &e0   = m_edges[ pre_e0 ];
		TWingEdge &e1   = m_edges[ pre_e1 ];
		int *p0vtx      = m_polys  [ pre_p0 ].idx;

		int vRightIdx = (p0vtx[0] != pre_v0 && p0vtx[0] != pre_v1) ? p0vtx[0] : 
                        (p0vtx[1] != pre_v0 && p0vtx[1] != pre_v1) ? p0vtx[1] : p0vtx[2];

		m_polys  [ pre_p0 ].ReplaceID( pre_v1, newVertexId);
		m_p_edges[ pre_p0 ].ReplaceID( pre_e1, newE0_Id   );

		m_polys  [newP0_Id].Set( newVertexId, pre_v1 , vRightIdx);
		m_p_edges[newP0_Id].Set(   newEpivId, pre_e1 ,  newE0_Id);

		m_v_edge[vRightIdx] = pre_e0;

		e.e[0] = pre_e0;
		e.e[1] = newE0_Id;

		e_nPiv.e[2] = pre_e1;
		e_nPiv.e[3] = newE0_Id;
		e_nPiv.p[1] = newP0_Id;
		
		e_n0.v[0] = vRightIdx  ; 
		e_n0.v[1] = newVertexId; 

		e_n0.e[0] = pre_e1    ; 
		e_n0.e[1] = newEpivId ; 
		e_n0.e[2] = pre_e0    ; 
		e_n0.e[3] = trgtEdgeId; 

		e_n0.p[0] = newP0_Id  ; 
		e_n0.p[1] = pre_p0    ;

		e0.ReplaceEdgeID(     pre_e1, newE0_Id);
		e1.ReplaceEdgeID(     pre_e0, newE0_Id);
		e1.ReplaceEdgeID( trgtEdgeId, newEpivId);
		e1.ReplacePolyID(     pre_p0, newP0_Id);
	}

	if( pre_p1 != -1 )
	{
		TWingEdge &e_n1 = m_edges[newE1_Id] ;
		TWingEdge &e2   = m_edges[ pre_e2 ] ;
		TWingEdge &e3   = m_edges[ pre_e3 ] ;
		int *p1vtx = m_polys[ pre_p1 ].idx; 

		int vLeftIdx = (p1vtx[0] != pre_v0 && p1vtx[0] != pre_v1) ? p1vtx[0] : 
                       (p1vtx[1] != pre_v0 && p1vtx[1] != pre_v1) ? p1vtx[1] : p1vtx[2];

		m_polys  [ pre_p1 ].ReplaceID( pre_v1, newVertexId);
		m_p_edges[ pre_p1 ].ReplaceID( pre_e3, newE1_Id   );

		m_polys  [newP1_Id].Set( newVertexId,  vLeftIdx,  pre_v1);
		m_p_edges[newP1_Id].Set(   newEpivId,  newE1_Id,  pre_e3);

		m_v_edge[ vLeftIdx ] = pre_e2;

		e.e[2] = pre_e2;
		e.e[3] = newE1_Id;

		e_nPiv.e[0] = pre_e3;
		e_nPiv.e[1] = newE1_Id;
		e_nPiv.p[0] = newP1_Id;
		
		e_n1.v[0] = vLeftIdx   ;
		e_n1.v[1] = newVertexId; 

		e_n1.e[0] = pre_e2    ; 
		e_n1.e[1] = trgtEdgeId; 
		e_n1.e[2] = pre_e3    ; 
		e_n1.e[3] = newEpivId ; 

		e_n1.p[0] = pre_p1    ; 
		e_n1.p[1] = newP1_Id  ;
		e2.ReplaceEdgeID(     pre_e3, newE1_Id);
		e3.ReplaceEdgeID(     pre_e2, newE1_Id);
		e3.ReplaceEdgeID( trgtEdgeId, newEpivId);
		e3.ReplacePolyID(     pre_p1, newP1_Id );
	}
	updateNormal();
	updateOneRing();
}

//TTriangle に使用されていない不必要なvertexを削除する
void TTriangleMesh::RemoveUnUsedVertices()
{
	if( m_eSize == 0 || m_edges == 0) updateWingEdge();

	vector< TVector3 > Vs         ; 
	vector< TTriangle  > Ps(m_pSize);

	Vs.reserve( m_vSize );
	int *vMap = new int[m_vSize];
	for( int i=0; i<m_vSize; ++i) if( m_v_edge[i] >= 0 ){ 
		Vs.push_back( m_verts[i]); 
		vMap[i] = (int)Vs.size()-1; 
	}
	for( int i=0; i<m_pSize; ++i) Ps[i].Set(vMap[ m_polys[i].idx[0] ], 
											vMap[ m_polys[i].idx[1] ], 
											vMap[ m_polys[i].idx[2] ] );
	initFromVertsPolys( Vs, Ps );
	delete[] vMap;
}


bool TTriangleMesh::RemoveEdgePossible(int edgeid) const{
	const TWingEdge& we = m_edges[edgeid] ;
	const TWingEdge& e0 = m_edges[we.e[0]] ; 
	const TWingEdge& e1 = m_edges[we.e[1]] ;
	const TWingEdge& e2 = m_edges[we.e[2]] ;
	const TWingEdge& e3 = m_edges[we.e[3]] ;

	const int v0_id = we.v[0] ;
	const int v1_id = we.v[1] ;
	if( bBoundaryVertex( v0_id ) || bBoundaryVertex(v1_id) ) return false;

	const int v0edge = m_v_edge[v0_id];
	const int v1edge = m_v_edge[v1_id];
	
	///// Check if this removal causes nonmanifold mesh.
	{
		int v_neigh_num[2]={0,0} ;

		int start_we_id = v0edge;
		int we_id = start_we_id ;
		do{
			TWingEdge *pwe = &m_edges[we_id] ;

			v_neigh_num[0]++ ;
			if( pwe->v[0]==v0_id ) we_id = pwe->e[0] ; 
			else                   we_id = pwe->e[3] ;
			
		} while ( we_id != start_we_id ) ;
		
		start_we_id = v1edge;
		we_id = start_we_id ;
		do{
			TWingEdge* pwe = &m_edges[we_id] ;
			
			v_neigh_num[1]++ ;
			if( pwe->v[0]==v1_id ) we_id = pwe->e[0] ;
			else                   we_id = pwe->e[3] ;
		} while ( we_id != start_we_id ) ;
		
		if( v_neigh_num[0] <= 3 || v_neigh_num[1] <= 3 ) return false ;

		
		int skip_vid[2] = { e0.v[0]==v0_id?e0.v[1]:e0.v[0] , e2.v[0]==v0_id?e2.v[1]:e2.v[0] } ;
		list<int> v0_neighbor_list ;
		// find v0 neighbors
		int start_id , id ;
		start_id = id = v0edge;
		do {
			TWingEdge* we = &m_edges[id] ;
			int neighborid = we->v[0]==v0_id?we->v[1]:we->v[0] ;
			if( neighborid != skip_vid[0] && neighborid != skip_vid[1] )
				v0_neighbor_list.push_back( neighborid ) ;

			id = (we->v[0]==v0_id?we->e[0]:we->e[3]) ;
		} while( id != start_id && id != -1 ) ;

		// check v1 neighbors
		start_id = id = v1edge;
		do {
			TWingEdge* we = &m_edges[id] ;
			int neighborid = we->v[0]==v1_id?we->v[1]:we->v[0] ;

			for(list<int>::iterator iit = v0_neighbor_list.begin() ; iit != v0_neighbor_list.end() ; iit++ )
				if( neighborid == *iit ) return false ;

			id = (we->v[0]==v1_id?we->e[0]:we->e[3]) ;
		} while( id != start_id && id != -1) ;
	}

	return true ;
}

// e = m_edges[edge_id], とe.e[1], e.e[3],   e.p[0], e.p[1],  e.v[1] を同時に削除
//
// edges[edgeid],edges[edges[edgeid].e[1]],edges[edges[edgeid].e[3]] and vertices[edges[edgeid].v[1]],
// polygons[edges[edgeid].p[0]],polygons[edges[edgeid].p[1]] are simultaneously removed. However, whether the actual instances are
// removed from array edges,vertices,polygons depends on the last argument.
// And also note position of vertices[edges[edgeid].v[0]] is not changed.
// if topologyOnly == true, the system does NOT remove targetVs/Ps/Es from array
bool TTriangleMesh::RemoveEdge(int edgeid, bool topologyOnly)
{
	if( !RemoveEdgePossible(edgeid) ) return false ;
	if( m_edges[edgeid].isBound ()  ) return false;
	const TWingEdge &we = m_edges[edgeid] ;
	TWingEdge &e0 = m_edges[we.e[0]], &e1 = m_edges[we.e[1]];//will remove
	TWingEdge &e2 = m_edges[we.e[2]], &e3 = m_edges[we.e[3]];//will remove

	const int v0_id    = we.v[0] ;
	const int v1_id    = we.v[1] ;
	m_verts[v0_id] = 0.5 * (m_verts[v0_id] + m_verts[v1_id] );
	vector< int > v1_Vs, v1_Ps, v1_Es; GetVsPsEsAroundVertex( v1_id, v1_Vs, v1_Ps, v1_Es);

	if(e0.v[0] == v0_id ){//e0.e[] e0.p[]の付け替え
		if( e1.v[0]==v1_id ){ e0.p[1] = e1.p[1];   e0.e[2] = (e1.e[2]==we.e[3] ? we.e[2] : e1.e[2]);   e0.e[3] = e1.e[3] ;} 
		else                { e0.p[1] = e1.p[0];   e0.e[2] = (e1.e[1]==we.e[3] ? we.e[2] : e1.e[1]);   e0.e[3] = e1.e[0] ;}
		m_polys[e0.p[1]].ReplaceID(v1_id,v0_id) ;
	} else               { 
		if( e1.v[0]==v1_id ){ e0.p[0] = e1.p[1];   e0.e[1] = (e1.e[2]==we.e[3] ? we.e[2] : e1.e[2]);   e0.e[0] = e1.e[3] ;}
		else                { e0.p[0] = e1.p[0];   e0.e[1] = (e1.e[1]==we.e[3] ? we.e[2] : e1.e[1]);   e0.e[0] = e1.e[0] ;}
		m_polys[e0.p[0]].ReplaceID(v1_id,v0_id);
	}
	if(e2.v[0] == v0_id ){//e2.e[] e2.p[]の付け替え
		if( e3.v[0]==v1_id ){ e2.p[0] = e3.p[0];   e2.e[0] = (e3.e[0]==we.e[1] ? we.e[0] : e3.e[0]) ;  e2.e[1] = e3.e[1] ;} 
		else                { e2.p[0] = e3.p[1];   e2.e[0] = (e3.e[3]==we.e[1] ? we.e[0] : e3.e[3]) ;  e2.e[1] = e3.e[2] ;}
		m_polys[e2.p[0]].ReplaceID(v1_id,v0_id);
	} else {
		if( e3.v[0]==v1_id ){ e2.p[1] = e3.p[0];   e2.e[3] = (e3.e[0]==we.e[1] ? we.e[0] : e3.e[0]) ;  e2.e[2] = e3.e[1] ;}
		else {                e2.p[1] = e3.p[1];   e2.e[3] = (e3.e[3]==we.e[1] ? we.e[0] : e3.e[3]) ;  e2.e[2] = e3.e[2] ;}
		m_polys[e2.p[1]].ReplaceID(v1_id,v0_id) ;
	}

	m_v_edge[v0_id] = we.e[0] ;
	m_v_edge[ e0.v[0]==v0_id? e0.v[1] : e0.v[0] ] = we.e[0] ; //左 wing vertexのone edge
	m_v_edge[ e2.v[0]==v0_id? e2.v[1] : e2.v[0] ] = we.e[2] ; //右 wing vertexのone edge

	//影響を受ける2 polygonsのエッジそのものと　そのエッジのwinged edgeの付け替え
	const int e1_p1_id = e0.p[e0.v[0]==v0_id?1:0]; //e1のwing polygonで，we.p[0]ではない方  *e0.pは書き換え済みでこれを利用
	const int e3_p0_id = e2.p[e2.v[0]==v0_id?0:1]; //e3のwing polygonで，we.p[1]ではない方  *e2.pは書き換え済みでこれを利用
	
	m_p_edges[ e1_p1_id ].ReplaceID( we.e[1], we.e[0] );
	m_p_edges[ e3_p0_id ].ReplaceID( we.e[3], we.e[2] );
	for(int i=0;i<3;i++ ) m_edges[ m_p_edges[ e1_p1_id ].idx[i] ].ReplaceEdgeID( we.e[1], we.e[0] );
	for(int i=0;i<3;i++ ) m_edges[ m_p_edges[ e3_p0_id ].idx[i] ].ReplaceEdgeID( we.e[3], we.e[2] );

	//v1の周りにあった edge polygonのv1への参照をv0につなぎかえる
	for( int i=0, s=(int)v1_Es.size(); i<s; ++i){
		if( v1_Ps[i] != we.p[0] && v1_Ps[i] != we.p[1]                       ) m_polys[v1_Ps[i]].ReplaceID   ( v1_id, v0_id);
		if( v1_Es[i] != edgeid  && v1_Es[i] != we.e[1] && v1_Es[i] != we.e[3]){
			m_edges[v1_Es[i]].ReplaceVtxID( v1_id, v0_id);
			if( m_edges[v1_Es[i]].v[0] > m_edges[v1_Es[i]].v[1] ) m_edges[v1_Es[i]].invert();
		}
	}
	//Remove  deleted instances from the arraies
	if(!topologyOnly)
	{
		set<int> r_vid; r_vid.insert( v1_id) ;
		set<int> r_pid; r_pid.insert( we.p[0] ); r_pid.insert( we.p[1] ); 
		set<int> r_eid; r_eid.insert( edgeid  ); r_eid.insert( we.e[1] ); r_eid.insert( we.e[3] );
		RemoveVsPsEs( r_vid, r_pid, r_eid);
	}
	return true ;
}

//RemoveShortEdges 平均より短い(*edge thresh rate)エッジを全て削除
//v_removable [ edge.v[0/1] = falseならそのエッジは消さない
//削除前verticesIDに対する削除後vertices id のmapをv_mapPreToNewに入れる
void TTriangleMesh::RemoveShortEdges( const vector<short> &v_removable, vector< int >  &v_mapPreToNew, double edgeThreshRate)
{
	double edgeLengthAve = CalcAverageEdgeLength();
	set<int> v_removed, e_removed, p_removed;

	while( true )
	{
		bool bRemoved = false;
		for( int i=0; i< m_eSize; ++i)
		{
			int v0 = m_edges[i].v[0];
			int v1 = m_edges[i].v[1];

			if( e_removed.find( i ) != e_removed.end() ) continue;
			if( !v_removable[ v0 ] || !v_removable[v1] ) continue;
		    if( t_distance( m_verts[ v0 ], m_verts[ v1 ] ) > edgeLengthAve * edgeThreshRate ) continue;

			//remove edge
			TWingEdge trgtEdgeCopy = m_edges[i];
			if( RemoveEdge( i, true ) )
			{
				v_removed.insert( trgtEdgeCopy.v[1] );
				p_removed.insert( trgtEdgeCopy.p[0] ); p_removed.insert( trgtEdgeCopy.p[1] );
				e_removed.insert( trgtEdgeCopy.e[1] ); e_removed.insert( trgtEdgeCopy.e[3] ); e_removed.insert( i );
				bRemoved = true;
			}
		}
		if( !bRemoved ) break;
	}
	
	v_mapPreToNew.clear();
	v_mapPreToNew.resize( m_vSize, -1 );
	int idx = 0;
	for( int i=0; i< getVnum(); ++i) if( v_removed.find( i ) == v_removed.end() ){
		v_mapPreToNew[i] = idx;
		idx++;
	}
	RemoveVsPsEs( v_removed, p_removed, e_removed );
}

//IDにより指定されたvertex/edge/polygonをそのまま削除する
//winged edgeが構築されている必要がある
//参照されている(使われている)vertex/edge/polygonの削除を行うと，後にエラーが起こる
void TTriangleMesh::RemoveVsPsEs(const set<int> &Vs, const set<int> &Ps, const set<int> &Es) 
{
	//if( m_eSize == 0 || m_edges == 0 ) return;
	int* vid_map = new int[ m_vSize ];
	int* eid_map = new int[ m_eSize ];
	int* pid_map = new int[ m_pSize ];
	TVector3 *newVerts  = new TVector3[m_vSize-Vs.size()];
	TVector3 *newVnorm  = new TVector3[m_vSize-Vs.size()];
	int        *newVedge  = new int       [m_vSize-Vs.size()];
	
	TTriangle  *newPolys  = new TTriangle [m_pSize-Ps.size()];
	TTriangle  *newPedge  = new TTriangle [m_pSize-Ps.size()];
	TVector3 *newPnorm  = new TVector3[m_pSize-Ps.size()];
	
	TWingEdge  *newEdges  = new TWingEdge [m_eSize-Es.size()];

	for( int i=0, id=0;i< m_vSize; ++i ) if( Vs.find(i) == Vs.end()){ newVerts[id] = m_verts [i] ;
																	  newVnorm[id] = m_v_norm[i] ;
																	  newVedge[id] = m_v_edge[i] ; vid_map[i] = id; ++id ;   }
	for( int i=0, id=0;i< m_pSize; ++i ) if( Ps.find(i) == Ps.end()){ newPolys[id] = m_polys  [i];
																	  newPedge[id] = m_p_edges[i];
																	  newPnorm[id] = m_p_norm [i]; pid_map[i] = id; ++id;      }
	for( int i=0, id=0;i< m_eSize; ++i ) if( Es.find(i) == Es.end()){ newEdges[id] = m_edges  [i]; eid_map[i] = id; ++id;      }

	m_vSize -= (int)Vs.size();
	m_pSize -= (int)Ps.size();
	m_eSize -= (int)Es.size();
	for( int i=0;i<m_vSize; ++i)  newVedge[i] = newVedge[i]==-1? -1 : eid_map[ newVedge[i] ];
	
	for( int i=0;i<m_pSize; ++i){ newPolys[i].idx[0]  = vid_map[newPolys[i].idx[0]] ; 
								  newPolys[i].idx[1]  = vid_map[newPolys[i].idx[1]] ;
								  newPolys[i].idx[2]  = vid_map[newPolys[i].idx[2]] ;
								  newPedge[i].idx[0]  = eid_map[newPedge[i].idx[0]] ;
								  newPedge[i].idx[1]  = eid_map[newPedge[i].idx[1]] ;
								  newPedge[i].idx[2]  = eid_map[newPedge[i].idx[2]] ; }
	for( int i=0;i<m_eSize; ++i){
								  TWingEdge& we = newEdges[i] ;
								  we.v[0] =                vid_map[we.v[0]];  we.v[1] =                vid_map[we.v[1]] ;
							 	  we.p[0] = we.p[0]==-1?-1:pid_map[we.p[0]];  we.p[1] = we.p[1]==-1?-1:pid_map[we.p[1]] ;
								  we.e[0] = we.e[0]==-1?-1:eid_map[we.e[0]];  we.e[1] = we.e[1]==-1?-1:eid_map[we.e[1]] ;
								  we.e[2] = we.e[2]==-1?-1:eid_map[we.e[2]];  we.e[3] = we.e[3]==-1?-1:eid_map[we.e[3]] ;}


	swap( newVerts, m_verts  ); if( newVerts != 0) delete[] newVerts;
	swap( newVnorm, m_v_norm ); if( newVnorm != 0) delete[] newVnorm;
	swap( newVedge, m_v_edge ); if( newVedge != 0) delete[] newVedge;

	swap( newPolys, m_polys  ); if( newPolys != 0) delete[] newPolys;
	swap( newPedge, m_p_edges); if( newPedge != 0) delete[] newPedge;
	swap( newPnorm, m_p_norm ); if( newPnorm != 0) delete[] newPnorm;

	swap( newEdges, m_edges  ); if( newEdges != 0) delete[] newEdges;

	delete[] vid_map ;
	delete[] pid_map ;
	delete[] eid_map ;

	updateNormal();
	updateOneRing();
}


//remove e = m_edges[ edgeId ]  e.v[1], e.p[0or1], e.e[1or3]
bool TTriangleMesh::RemoveEdge_onBound( int edgeId)
{
	const TWingEdge &e = m_edges[ edgeId ];
	if( e.p[0] != -1 && e.p[1] != -1 ) return false; 
	if( e.p[0] == -1 && e.p[1] == -1 ) return false; 

	//remove vertex * 1,  polygon * 1
	const int removedVid = e.v[1];
	const int newVid     = e.v[0];
	const int removedPid = ( e.p[0] != -1) ? e.p[0] : e.p[1];

	vector< int > v1_Ps; GetPsAroundVertex( e.v[1], v1_Ps );//v1がとられて、v0にマージされる

	TTriangle &po = m_polys[removedPid];
	if( GetDegreeOfVertex( po.idx[0] ) <= (bBoundaryVertex( po.idx[0]) ? 2 : 3)) return false;
	if( GetDegreeOfVertex( po.idx[1] ) <= (bBoundaryVertex( po.idx[1]) ? 2 : 3)) return false;
	if( GetDegreeOfVertex( po.idx[2] ) <= (bBoundaryVertex( po.idx[2]) ? 2 : 3)) return false;

	m_verts[e.v[0]] = 0.5 * (m_verts[e.v[0]] + m_verts[e.v[1]]);

	//removed vertexの付け替え//
	for( int i=0, s=(int) v1_Ps.size(); i<s; ++i) if( v1_Ps[i] != removedPid )
		m_polys[ v1_Ps[i] ].ReplaceID( removedVid, newVid );

	vector< TVector3 > Vs; Vs.reserve( m_vSize-1 );
	vector< TTriangle  > Ps; Ps.reserve( m_pSize-1 );
	for( int i=0; i<m_vSize; ++i) if( i != removedVid ) Vs.push_back( m_verts[i] );
	for( int i=0; i<m_pSize; ++i) if( i != removedPid ) Ps.push_back( m_polys[i] );


	//vertexが取り除かれるので,Vidをずらす
	for( int i=0; i<m_pSize-1; ++i)
	{
		TTriangle &p = Ps[i];
		if( p.idx[0] > removedVid ) --p.idx[0];
		if( p.idx[1] > removedVid ) --p.idx[1];
		if( p.idx[2] > removedVid ) --p.idx[2];
	}
	initFromVertsPolys( Vs, Ps );
	return true;
}

//////Region growing!////////////////////////////////////////////////////////////////////////
//lavel値は 0 ～ : 連続領域にlabelingする
//戻り値はlabel数
int TTriangleMesh::calcLinkedLabel( vector<int> &v_labels )
{
	if( m_vSize != m_v_1ring.size() ) updateOneRing();

	v_labels.clear(); v_labels.resize( m_vSize, -1 );

	int labelId = 0;
	for( ; true; ++labelId )
	{
		//search pivot id//
		int pivVertex = -1;
		for( int i=0; i<m_vSize; ++i) if( v_labels[i] == -1 ) { pivVertex = i; break;}
		if( pivVertex == -1 ) break;

		//growth//
		list< int > growthFront;
		growthFront.push_back( pivVertex );
		v_labels[pivVertex] = labelId;

		while( !growthFront.empty() )
		{
			int idx = growthFront.front(); growthFront.pop_front();
			for( int i=0;i< (int)m_v_1ring[idx].size(); ++i) if( v_labels[ m_v_1ring[idx][i] ] == -1 )
			{
				v_labels[ m_v_1ring[idx][i] ] = labelId;
				growthFront.push_back( m_v_1ring[idx][i] );
			}
		}
	}
	return labelId;
}
//lavel値は 0 ～ : 連続領域にlabelingする
int TTriangleMesh::calcLinkedLabel_poly( vector<int> &p_labels )
{
	if( m_vSize != m_v_1ring.size() ) updateOneRing();
	p_labels.clear(); p_labels.resize( m_pSize, -1 );

	int labelId = 0;
	for( ; true; ++labelId )
	{
		//search pivot id//
		int pivPolygon = -1;
		for( int i=0; i< m_pSize; ++i) if( p_labels[i] == -1 ) { pivPolygon = i; break;}
		if( pivPolygon == -1 ) break;

		//growth
		list< int > growthFront;
		growthFront.push_back( pivPolygon );
		p_labels[pivPolygon] = labelId;

		while( !growthFront.empty() )
		{
			int idx = growthFront.front(); growthFront.pop_front();
			
			vector<int> Ps; GetPsAroundVertex( m_polys[idx].idx[0], Ps );
			for( int i=0;i<(int)Ps.size(); ++i) if( p_labels[ Ps[i] ] == -1 ){ p_labels[ Ps[i] ] = labelId; growthFront.push_back( Ps[i] );}
			Ps.clear();     GetPsAroundVertex( m_polys[idx].idx[1], Ps );
			for( int i=0;i<(int)Ps.size(); ++i) if( p_labels[ Ps[i] ] == -1 ){ p_labels[ Ps[i] ] = labelId; growthFront.push_back( Ps[i] );}
			Ps.clear();     GetPsAroundVertex( m_polys[idx].idx[2], Ps );
			for( int i=0;i<(int)Ps.size(); ++i) if( p_labels[ Ps[i] ] == -1 ){ p_labels[ Ps[i] ] = labelId; growthFront.push_back( Ps[i] );}
		}
	}
	return labelId;
}

//////////////calcRoiMesh////////////////////////////////////////////////////////////////////////////
//Region growing from SeedVertices( number of growing iteration is defined by roiSize
//SeedVertsからroiSize分だけregion growingさせた Region of interest meshを生成する
//roiVs, roiPsに頂点情報とtriangle情報が格納される
//vMap_ROItoFUL, vMap_FULtoROIには ROI mesh  vertex <--> このmeshのvertexの頂点IDの変換mapが挿入される
void TTriangleMesh::calcRoiMesh( const	vector<int>		  &seedVertIdx, 
								        int				   roiSize,
										vector<TVector3 > &roiVs, 
										vector<TTriangle> &roiPs,
										vector<int>       &vmap_ROItoFUL,
										vector<int>       &vmap_FULtoROI) 
{
	if( m_vSize != m_v_1ring.size() ) updateOneRing();

	vector<int> bRoiVerts( m_vSize, 0 );//0:init, 2:visiting  ループで訪れてる、1:visited 

	for( int i=0, s=(int)seedVertIdx.size(); i<s; ++i) bRoiVerts[ seedVertIdx[i] ] = 1;

	for( int count = 0; count < roiSize; ++count)
	{
		for( int i = 0; i<m_vSize;                        ++i) if( bRoiVerts[       i       ] == 1 )
		for( int k = 0,ki=(int)m_v_1ring[i].size(); k<ki; ++k) if( bRoiVerts[m_v_1ring[i][k]] != 1 ) bRoiVerts[ m_v_1ring[i][k] ] = 2;

		for( int i = 0; i<m_vSize; ++i) if( bRoiVerts[i]==2 ) bRoiVerts[i] = 1;
	}

	//ROI meshを構築する////////////////////////////////////////////////////////////
	vmap_FULtoROI.clear();  vmap_FULtoROI.resize ( m_vSize, -1 );
	vmap_ROItoFUL.clear();  vmap_ROItoFUL.reserve( m_vSize     );
	roiVs.clear();
	roiPs.clear();

	for( int i = 0; i < m_vSize; ++i) if( bRoiVerts[i] )
	{
		vmap_FULtoROI[i] = (int)roiVs.size();
		roiVs        .push_back( m_verts[i] );
		vmap_ROItoFUL.push_back( i          );
	}

	for( int i=0; i<m_pSize; ++i)
		if( bRoiVerts[ m_polys[i].idx[0] ] && bRoiVerts[ m_polys[i].idx[1] ] && bRoiVerts[ m_polys[i].idx[2] ])
			roiPs.push_back( TTriangle(	vmap_FULtoROI[ m_polys[i].idx[0] ], 
			                            vmap_FULtoROI[ m_polys[i].idx[1] ], 
										vmap_FULtoROI[ m_polys[i].idx[2] ]) );
}

//getRingNeighborhood/////////////////////////////////////////////////////////////////
//m_verts[ vIdx ] から ringSize分region growingを行った N-ring neighborhood を生成する
//ringVerts 内に対象となる Vertex idxが挿入される(vIdx自信を含む)
void TTriangleMesh::getRingNeighborhood( int vIdx, int ringSize, set<int> &ringVerts)
{
	if( m_vSize != m_v_1ring.size() ) updateOneRing();
	ringVerts.clear();
	set<int> borderVerts; 
	borderVerts.insert( vIdx );

	for( int count = -1; count < ringSize; ++count)
	{
		set<int> newBorder;
		for( set<int>::iterator it = borderVerts.begin(); it != borderVerts.end(); ++it)
		{
			ringVerts.insert( *it );
			for( int i=0, s = (int)m_v_1ring[*it].size(); i<s; ++i)
			{
				int trgtI = m_v_1ring[*it][i];
				if( ringVerts  .find( trgtI ) == ringVerts  .end() && borderVerts.find( trgtI ) == borderVerts.end() ) newBorder.insert( trgtI );
			}
		}
		borderVerts = newBorder;
	}
}

//calcNearestPolygon_insideCheck////////////////////////////////////////////////
//点pから一番近いtriangleを検索し、そのtriangleの法線方向を利用して 点ｐをtriangleに射影
//射影した点がtriangleの内部に入るときのみtrueを返す
bool TTriangleMesh::calcNearestPolygon_insideCheck(const TVector3 &p, TVector3 &nearPos, int &polyIdx,  double &dist) const
{
	int vIdx; getNearestVertexIdx(p, vIdx );
	vector<int > Ps; GetPsAroundVertex( vIdx, Ps );
	
	TVector3 pos, bary, P0,P1,P2, zero(0,0,0), one(1,1,1);
	for( int i=0; i<(int)Ps.size(); ++i)
	{
		P0  = m_verts[ m_polys[ Ps[i] ].idx[0] ];
		P1  = m_verts[ m_polys[ Ps[i] ].idx[1] ];
		P2  = m_verts[ m_polys[ Ps[i] ].idx[2] ];
		const TVector3 &N = m_p_norm [ Ps[i] ];

		pos = p;
		pos.Add_CoefMultVec( -t_V1subtV2_multV3( pos, P0, N), N);
		pos *= 10.0; P0  *= 10.0; P1  *= 10.0; P2  *= 10.0;
		bool tf = t_calcBaryCoordOfPositionOnPoly( pos, P0, P1, P2, bary) ;
		pos *= 0.1 ; P0  *= 0.1 ; P1  *= 0.1 ; P2  *= 0.1 ;

		if( !tf ) continue;

		if( t_isInWindow3D_MinMax( zero, one, bary )  ){
			nearPos = pos;
			polyIdx = Ps[i];
			return true;
		}
	}
	return false;
}

//Hole filling///////////////////////////////////////////////////////////////////////////////////////

///////////////////GetSequentialBoundVsEs/////////////////////////////////////////
// Seed Vtxからスタートし、boundary上にある一連のvertex + edgesを取得する
//boundaryが枝分かれするようなときには, falseを返す
//方向は,隣接triangleのnormal方向を、上向きにとった時に、
//右側に左側にbouondary、右側にtriangleをおく
//
//つまり、穴にたいして左回りになる
//       x1
//       |\
//       |  \
//こっち↑   \
//bound  |    \   normalは紙面手前向き
//       |     \
//	     |------ x2
//       x0

bool TTriangleMesh::GetSequentialBoundVsEs( int seedVtx, vector<int> &seqVs, vector<int> &seqEs)
{
	seqVs.clear(); 
	seqEs.clear();
	if( !bBoundaryVertex( seedVtx ) ) return false;
	seqVs.push_back( seedVtx );

	{   //進行方向確定
		vector<int> Es; GetEsAroundVertex( seedVtx, Es );
		vector<int> bEs;
		for( int i=0, s=(int) Es.size(); i<s; ++i) if( m_edges[Es[i]].isBound() ) bEs.push_back( Es[i] );
		if( (int) bEs.size() != 2 ) {fprintf(stderr,"strange ! 789asfd21\n"); return false;}
		
		const TWingEdge &we = m_edges[ bEs.front() ];
		
		int polyId = (we.p[0] !=-1) ? we.p[0] : we.p[1]; 
		const TTriangle &p = m_polys[ polyId      ]; 
		TVector3 norm; CalcNormalOfOnePolygon( polyId, norm);
		
		//polygonのx0-x2とは無関係
		TVector3 &e_v0 = m_verts[  seedVtx ]; //edgeのv[ 0 ]
		TVector3 &e_v1 = m_verts[ (seedVtx  != we.v[0]                        ) ? we.v[0]  : we.v[1] ]; //edgeのv[ 1 ]
		TVector3 &p_v2 = m_verts[ (p.idx[0] != we.v[0] && p.idx[0] != we.v[1] ) ? p.idx[0] :
			                        (p.idx[1] != we.v[0] && p.idx[1] != we.v[1] ) ? p.idx[1] : p.idx[2]]; //edgeに隣接するpolygonの頂点
		
		bool bFored =  (p_v2 - e_v0) * ( (e_v1 - e_v0) ^ norm) >= 0 ;

		TWingEdge &e = m_edges[  bFored ? bEs.front():bEs.back()  ]; 
	
		seqEs.push_back( bFored  ?  bEs.front() : bEs.back()  );
		seqVs.push_back( (seqVs.back() != e.v[0])?e.v[0]:e.v[1] );
	}

	//trace boundary
	while( true )
	{
		if( (int) seqVs.size() > m_vSize ) return false;//1点に４本以上boundary edgeがある場合，無限ループに落ちる場合がある 
		vector<int> Es; GetEsAroundVertex( seqVs.back(), Es );

		//check//
		int bEdgeCounter = 0;
		for( int i=0; i<(int) Es.size(); ++i) if( m_edges[Es[i]].isBound() ) bEdgeCounter++;
		if( bEdgeCounter != 2 ) { fprintf(stderr,"strange ! 789asfd21\n");return false;}

		for( int i=0; i<(int) Es.size(); ++i) if( m_edges[Es[i]].isBound() && Es[i] != seqEs.back()) {
			seqEs.push_back( Es[i]); break; 
		}

		TWingEdge &e = m_edges[ seqEs.back() ];

		int newVid = (e.v[0] != seqVs.back() ) ? e.v[0] : e.v[1]; 
		if( newVid == seedVtx ) return true;
	
		seqVs.push_back( newVid );
	}

	return false;
}


//////fillHoles_noExtraVertex/////////////////////////////////////////////////////
//頂点を追加しないsimpleなhole filling
//seedIdxを含むholeを, 頂点を新たに追加せえず，triangle fanにより埋める
//凸でない穴についてはtopology的には埋まるが，重なりなどが起きる可能性は十分ある
bool TTriangleMesh::fillHole_noExtraVertex(int seedIdx)
{
	vector< int > seqVs, seqEs;
	if( !GetSequentialBoundVsEs( seedIdx, seqVs, seqEs ) ) return false; 

	//重心に一番近い頂点を利用(一応)
	TVector3 gc; for( int i=0;i<(int)seqVs.size(); ++i) gc += m_verts[seqVs[i]]; gc /= (double) seqVs.size();

	int pivIdx = -1;
	double dist = 100000000;
	for( int i=0;i<(int)seqVs.size(); ++i) {
		double d = t_distance_sq( gc, m_verts[ seqVs[i] ]);
		if( d < dist ){dist=d; pivIdx=seqVs[i]; }
	}

	vector< TVector3 > Vs( m_vSize ); for( int i=0; i<m_vSize; ++i) Vs[i]=m_verts[i] ;
	vector< TTriangle  > Ps( m_pSize ); for( int i=0; i<m_pSize; ++i) Ps[i]=m_polys[i] ;

	for( int i=0;i<(int)seqVs.size(); ++i) 
	{
		int id1 = seqVs[ (i==0)?seqVs.size()-1 : i-1];
		int id2 = seqVs[        i                   ];
		if( id1 != pivIdx && id2!= pivIdx ) Ps.push_back( TTriangle(pivIdx, id1, id2) );
	}
	initFromVertsPolys( Vs, Ps );
	return true;
}
//全てのholeを埋める
bool TTriangleMesh::fillHoles_noExtraVertex()
{
	while( true )
	{
		int seedIdx = -1;		
		for( int i=0; i<m_eSize; ++i) if( m_edges[i].isBound() ) { seedIdx = m_edges[i].v[0]; break;}
		if( seedIdx == -1 ) break;
		if( !fillHole_noExtraVertex( seedIdx ) ) return false;
		if( !checkWingedEdge()                 ) return false;
	}
	return true;
}
bool TTriangleMesh::fillHoles( )
{
	while( true )
	{
		int seedIdx = -1;		
		for( int i=0; i<m_eSize; ++i) if( m_edges[i].isBound() ) { seedIdx = m_edges[i].v[0]; break;}
		if( seedIdx == -1 ) break;
		if( !fillHole( seedIdx ) ) return false;
	}
	return true;
}


///////fillHole//////////////////////////////////////////////////////
//seedVtxを含むHoleをちょうどいいgrid mesh patchを用いることで埋める
//vertexのindexは変化しない(新たなvertexは、後ろに追加される/////////
bool TTriangleMesh::fillHole( int seedVtx)
{
	fprintf( stderr, "startFilling Hole\n");
	vector< TVector3 > Vs( m_vSize ); for( int i=0; i<m_vSize; ++i) Vs[i]=m_verts[i] ;
	vector< TTriangle  > Ps( m_pSize ); for( int i=0; i<m_pSize; ++i) Ps[i]=m_polys[i] ;

	//Get boundary sequens(バウンダリは左回り(holeを左手にみる)) 
	vector< int > seqVs, seqEs;
	if( !GetSequentialBoundVsEs( seedVtx, seqVs, seqEs ) ) return false; 

	vector<int> newPolygons;
	if(     seqVs.size() <= 2 ) return false;
	else if(seqVs.size() <= 3 )//新たな頂点必要ない
	{
		fprintf( stderr, " fill hole by TRIANGLE !!!!!!!!!!!!!!!!!\n");
		Ps.push_back( TTriangle(seqVs[0], seqVs[1], seqVs[2] )); 
	}
	else if(seqVs.size() <= 9 )//single vertexでつなぐ
	{
		fprintf( stderr, " fill hole by SINGLE POINT !!!!!!!!!!!!!!!!!\n");

		int newVid = (int)Vs.size();
		Vs.push_back( TVector3() );

		for( int i=0; i<(int)seqVs.size();++i)
		{
			Ps.push_back( TTriangle( newVid, seqVs[ i==0?seqVs.size()-1:i-1], seqVs[ i ]));
			Vs.back() += m_verts[ seqVs[i] ];
		}
		Vs.back() /= (double)seqVs.size();
	}
	else //pachを作成
	{
		fprintf( stderr, " fill hole by  PATCH !!!!!!!!!!!!!!!!!\n");

		vector<TVector3> verts2D(seqVs.size()); 
		for( int i=0; i<(int)seqVs.size();++i) verts2D[i] = Vs[seqVs[i]];

		TVector3 tmp, tmp1, dir, gCenter;
		TMatrix16 rotTo2D, rotFrom2D;

		t_PCA_3D( verts2D, tmp, tmp1, dir );//dirがboundaryを平面近似した時の平面のnormal
		rotTo2D.RotateFromVectorToVector( dir, TVector3(0,0,1) );
		rotFrom2D = rotTo2D; 
		rotFrom2D.GetInvMatrix_Self();
		
		t_stroke_getGravityCenter( verts2D, gCenter );

		for( int i=0; i<(int)verts2D.size();++i){
			verts2D[i] -= gCenter;
			verts2D[i] = rotTo2D * verts2D[i]; 
			verts2D[i].data[2]  = 0;
		}

		double gridR = t_stroke_Length( verts2D ) / ((double)verts2D.size()-1.0);

		TTriangleMesh patch;
		vector<int>    mapVerts2D_to_pachBound;//verts2Dのidx --> pachboundaryのidx　というmap
		patch.initAs2DGridPatch_boundVCorrespond( verts2D, gridR, mapVerts2D_to_pachBound);//pach作成

		//縫い合わせ boundaryの方向と数を合わせて貼り合わせる
		patch.TransformAllVertices( rotFrom2D );
		patch.TranslateAllVertices( gCenter   );

		vector< int > mapPachVtxToNewVtx( patch.getVnum(), -1);
		for( int i=0; i<(int) mapVerts2D_to_pachBound.size(); ++i) mapPachVtxToNewVtx[ mapVerts2D_to_pachBound[i] ] = seqVs[i];

		for( int i=0; i<(int) patch.getVnum() ; ++i) if( mapPachVtxToNewVtx[i] == -1 ){
			mapPachVtxToNewVtx[i] = (int)Vs.size();
			Vs.push_back( patch.m_verts[i] );
		}

		for( int i=0; i<(int) patch.getPnum(); ++i){
			int *vtx = patch.m_polys[i].idx;
			Ps.push_back( TTriangle( mapPachVtxToNewVtx[ vtx[0] ], mapPachVtxToNewVtx[ vtx[1] ], mapPachVtxToNewVtx[ vtx[2] ] ) );
		}
	}
	initFromVertsPolys( Vs, Ps );

	fprintf( stderr, "finish Hole filling \n");

	return true;
}

///////////initAsPatch////////////////////////////////////////////////////////////////////
//2Dの輪郭 boundPointsより一回り小さいgrid状のpatchを生成する
//boundaryの位置がちょうど会う物ではない
//     0   1   2   3   4   5   6   7   8   
//   0   1   2   3   4   5   6   7   8   9    
//     0   1   2   3   4   5   6   7   8
//   0   1   2   3   4   5   6   7   8   9    
//     0   1   2   3   4   5   6   7   8
//   0   1   2   3   4   5   6   7   8   9    
//
//   const vector<TVector3> boundPoints　　　　バウンダリ 
//         double             gridR      ,       gridのedgeサイズ
bool TTriangleMesh::initAs2DGridPatch(const vector<TVector3> boundPoints2D, double gridR)
{
	fprintf( stderr, "start pach creation\n" );
	//bounding box計算
	double bound[4]; t_stroke2D_BoundBox( boundPoints2D, bound);
	bound[0] -= 4 * gridR; bound[1] -= 4 * gridR; bound[2] += 4 * gridR;   bound[3] += 4 * gridR;

	double    offset = 0.25 * gridR;
	const int n_X    = (int) ((bound[2] - bound[0]) / gridR );
	const int n_Y    = (int) ((bound[3] - bound[1]) / gridR );

	//grid作成
	vector<vector< TVector3>> g_Verts  ( n_Y, vector<TVector3>( n_X     ) );//grid頂点 位置
	vector<vector< int       >> g_bInside( n_Y, vector<    int   >( n_X,  0 ) );//grid頂点 boundの内部に入るかどうか 0:初期値 1:in -1:out
	vector<vector< int       >> g_VertIdx( n_Y, vector<    int   >( n_X, -1 ) );//grid頂点 内部頂点群でのindex
	//位置
	for( int iy = 0; iy < n_Y; ++iy, offset *= -1 )
	for( int ix = 0; ix < n_X; ++ix               ) g_Verts[iy][ix].Set( bound[0] + ix * gridR + offset, bound[1] + iy * gridR         , 0 );

	//内部の点をOn 1 にして外部の点をoff - 1 にする　//近い点を探す//
	for( int gy = 0; gy < n_Y; ++gy ) 
	for( int gx = 0; gx < n_X; ++gx ){
		for( int stri = 0; stri < (int) boundPoints2D.size(); ++stri)
		{
			const TVector3 &lineP0 = boundPoints2D[stri];
			const TVector3 &lineP1 = boundPoints2D[ (stri == 0 ) ? (int) boundPoints2D.size() - 1: stri - 1 ];
			if( t_distPointToLineSegment_2D( g_Verts[gy][gx], lineP0, lineP1) < gridR + 0.0001){
				g_bInside[gy][gx] = true;
				break;
			}
		}
	}

	bool isCounterCrockWise = false;

	for( int iy = 0; iy < n_Y; ++iy)
	for( int ix = 0; ix < n_X; ++ix) if( g_bInside[iy][ix] == 0 )
	{
		double angleSum = 0;
		bool isInside = t_stroke2D_pointInsideStroke( g_Verts[iy][ix], boundPoints2D, angleSum);
		
		g_bInside[iy][ix] = (isInside) ? 1 : -1;
		if( isInside ) isCounterCrockWise = (angleSum > 0);

		for( int iix = ix + 1; iix < n_X; ++iix){
			if( g_bInside[iy][iix] != 0 ) break;
			g_bInside[iy][iix] = (isInside) ? 1 : -1;
		}
	}

	//バウンダリ上で２近傍のみがOnのvertexはとる
	while( true )
	{
		bool isremoved = false;
		for( int iy = 1; iy < n_Y-1; ++iy) 
		for( int ix = 1; ix < n_X-1; ++ix) if( g_bInside[iy][ix] == 1 )
		{
			int sum = 0;
			if( iy%2 == 0) sum = ((g_bInside[ iy ][ix+1]==1)?1:0) + ((g_bInside[ iy ][ix-1]==1)?1:0)  +  //みぎに出っ張っている
								 ((g_bInside[iy+1][ ix ]==1)?1:0) + ((g_bInside[iy+1][ix+1]==1)?1:0)  +  
								 ((g_bInside[iy-1][ ix ]==1)?1:0) + ((g_bInside[iy-1][ix+1]==1)?1:0)  ;
			else           sum = ((g_bInside[ iy ][ix+1]==1)?1:0) + ((g_bInside[ iy ][ix-1]==1)?1:0)  +  
								 ((g_bInside[iy+1][ ix ]==1)?1:0) + ((g_bInside[iy+1][ix-1]==1)?1:0)  +  
								 ((g_bInside[iy-1][ ix ]==1)?1:0) + ((g_bInside[iy-1][ix-1]==1)?1:0)  ;
			if( sum == 2 ){
				g_bInside[iy][ix] = -1;
				isremoved = true;
			}
		}
		if(! isremoved ) break;
	}

	////////////////////////mesh生成////////////////////////////////////////////////////////////////////////////////////
	vector<TVector3> Vs;
	vector<TTriangle > Ps;
	for( int iy = 0; iy < n_Y; ++iy)
	for( int ix = 0; ix < n_X; ++ix) if( g_bInside[iy][ix]  == 1 ){
		g_VertIdx[iy][ix] = (int)Vs.size();
		Vs.push_back( g_Verts[iy][ix] );
	}
	//polygonを入れていく
	for( int iy = 0; iy < n_Y -1 ; ++iy){
		for( int ix = 0; ix < n_X-1; ++ix)
		{
			if( g_bInside[iy  ][ ix ] == 1 && g_bInside[iy+1][ ix ] == 1 &&  g_bInside[iy+1][ix+1] == 1) Ps.push_back( TTriangle(g_VertIdx[iy  ][ ix ], g_VertIdx[iy+1][ix+1], g_VertIdx[iy+1][ ix ]) );
			if( g_bInside[iy  ][ ix ] == 1 && g_bInside[iy+1][ix+1] == 1 &&  g_bInside[ iy ][ix+1] == 1) Ps.push_back( TTriangle(g_VertIdx[iy  ][ ix ], g_VertIdx[ iy ][ix+1], g_VertIdx[iy+1][ix+1]) );
		}
		++iy; if( iy >= n_Y -1 ) break;
		for( int ix = 1; ix < n_X; ++ix)
		{
			if( g_bInside[iy  ][ ix ] == 1 && g_bInside[ iy ][ix-1] == 1 &&  g_bInside[iy+1][ix-1] == 1 ) Ps.push_back( TTriangle( g_VertIdx[iy  ][ ix ], g_VertIdx[iy+1][ix-1], g_VertIdx[iy  ][ix-1]) );
			if( g_bInside[iy  ][ ix ] == 1 && g_bInside[iy+1][ix-1] == 1 &&  g_bInside[iy+1][ ix ] == 1 ) Ps.push_back( TTriangle( g_VertIdx[iy  ][ ix ], g_VertIdx[iy+1][ ix ], g_VertIdx[iy+1][ix-1]) );
		}
	}
	if( (int) Ps.size() == 0 ) {Clear(); return false; }

	if( !initFromVertsPolys( Vs, Ps ) ) return false;
	if( !isCounterCrockWise ) Flip();
	return true;
}
///////////initAsPatch_boundVCorrespond////////////////////////////////////////////////////////////////////
//2Dの輪郭 boundPointsをぴったり埋める grid状のpatchを生成する
//生成されるpatchの輪郭の頂点数がちょうど　boundPoints　の頂点数と同じになり、対応がmapPointsToVertsに入る
//     0   1   2   3   4   5   6   7   8   
//   0   1   2   3   4   5   6   7   8   9    
//     0   1   2   3   4   5   6   7   8
//   0   1   2   3   4   5   6   7   8   9    
//     0   1   2   3   4   5   6   7   8
//   0   1   2   3   4   5   6   7   8   9    
//
//   holeに対し boundaryは左回りになっている, 
//  生成されたmeshの normal はzの正方向を向くように
//
//   const vector<TVector3> boundPoints　　　　バウンダリ 
//         double             gridR      ,       gridのedgeサイズ
//		 vector<int>        &mapPointsToVerts　バウンダリの各頂点に対応するvertexのindex

bool TTriangleMesh::initAs2DGridPatch_boundVCorrespond( const vector<TVector3> boundPoints2D, double gridR, vector<int> &mapPointsToVerts)
{
	if( !initAs2DGridPatch( boundPoints2D, gridR ) ) return false;//construct slightly small grid patch

	//////////boundaryの頂点数を合わせる Split or Delete edges////////////////////////////////////////////////////////////////////////////
	const int boundarySize = (int) boundPoints2D.size();
	int pachBoundSize = 0;
	for( int i=0;i<m_eSize; ++i) if( m_edges[i].isBound() ) ++pachBoundSize;

	if( pachBoundSize < boundarySize ){ //pachBoundSizeが足りないので boundary edgeを分割する
		for( int i=0; i < boundarySize - pachBoundSize; ++i){
			int edgeIdx = -1;
			while(true){
				edgeIdx = (int) ((rand() / (double)RAND_MAX) * (m_eSize - 1));
				if( m_edges[edgeIdx].isBound()) break;
			}
			int newVertexId;
			splitEdge( edgeIdx, 0.5*(m_verts[ m_edges[edgeIdx].v[0]] + m_verts[m_edges[edgeIdx].v[1]]), newVertexId);
		}
	}else if( pachBoundSize > boundarySize ){ //pachBoundSizeが多すぎるのでboundary edgeを消去する
		for( int i=0; i < pachBoundSize - boundarySize; ++i){
			int edgeIdx = -1;
			while(true){
				edgeIdx = (int) ((rand() / (double)RAND_MAX) * (m_eSize - 1));
				if( m_edges[edgeIdx].isBound() ) break;
			}
			RemoveEdge_onBound( edgeIdx );
		}
	}

	//頂点間の対応付け////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int seedVIdx = -1;
	double dist  = DBL_MAX;
	for( int i=0; i<m_vSize; ++i)if( bBoundaryVertex(i) ){
		double d = t_distance( boundPoints2D.front(), m_verts[i] );
		if( d<dist ){dist=d; seedVIdx = i;}
	}

	vector<int> seqVs, seqEs;
	GetSequentialBoundVsEs( seedVIdx, seqVs, seqEs );//boundary(hole)を左手に見ながら進行するので、patchから見ると、法線を手前方向にとったら，右回りになる

	mapPointsToVerts.resize( seqVs.size() );
	mapPointsToVerts[0] = seqVs[0];
	for( int i=1; i<(int) seqVs.size(); ++i) mapPointsToVerts[i] = seqVs[seqVs.size()-i];
	return true;
}


void TTriangleMesh::SaveVsPsAsText( FILE *fp)
{
	fprintf( fp, "polgonModel %d %d\n", m_vSize, m_pSize );
	for( int i=0; i<m_vSize; ++i) fprintf( fp, "%f %f %f\n", m_verts[i].data[0], m_verts[i].data[1], m_verts[i].data[2] );
	for( int i=0; i<m_pSize; ++i) fprintf( fp, "%d %d %d\n", m_polys[i].idx [0], m_polys[i].idx [1], m_polys[i].idx [2] );
}

//古いバージョンのfileを読むときはエラーが出てfalseを返す
bool TTriangleMesh::LoadVsPsAsText( FILE *fp )
{
	int vSize, pSize;
	char buf[256];
	fscanf( fp, "%s%d%d", buf, &vSize, &pSize ); fprintf( stderr, "loading file %s vSize %d pSize%d\n", buf, vSize, pSize );
	vector< TVector3> Vs(vSize);
	vector< TTriangle > Ps(pSize);


	for( int i=0; i<vSize; ++i) if( fscanf( fp, "%lf%lf%lf", &Vs[i].data[0], &Vs[i].data[1], &Vs[i].data[2]) == EOF ) return false;
	for( int i=0; i<pSize; ++i) if( fscanf( fp, "%d%d%d"   , &Ps[i].idx [0], &Ps[i].idx [1], &Ps[i].idx [2]) == EOF ) return false;
	initFromVertsPolys( Vs, Ps );
	return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//global functions//////////////////////////////////////////////////////////////////////////////////////////////////////
static int TTRIMESH_SPHERE_RES = 7;

//Render sphere objects (center: p, radius:r) without surface normal
void t_drawSphere(double r, const TVector3 &p)
{
	static bool first = true;
	static TTriangleMesh sph;
	if( first ) { first = false; sph.initAsSphere( 1 , TTRIMESH_SPHERE_RES, TTRIMESH_SPHERE_RES ); }

	glPushMatrix();	
		glTranslated(  p.data[0],  p.data[1],  p.data[2] );
		glScaled( r,r,r);
		glBegin( GL_TRIANGLES );
		for( int i= 0; i < sph.getPnum(); ++i){
			glVertex3dv( sph.m_verts[ sph.m_polys[i].idx[0]].data );
			glVertex3dv( sph.m_verts[ sph.m_polys[i].idx[1]].data );
			glVertex3dv( sph.m_verts[ sph.m_polys[i].idx[2]].data );
		}	
		glEnd();
	glPopMatrix();
}
//Render edges of a sphere objects (center: p, radius:r) 
void t_drawSphereEdges(double r, const TVector3 &p)
{
	static bool first = true;
	static TTriangleMesh sph;
	if( first ) { first = false; sph.initAsSphere( 1 , TTRIMESH_SPHERE_RES, TTRIMESH_SPHERE_RES ); }	
	glPushMatrix();
		glTranslated(  p.data[0],  p.data[1],  p.data[2] );
		glScaled( r,r,r);
		sph.drawEdges();
	glPopMatrix();
}
//Render sphere objects (center: p, radius:r) with surface normal
void t_drawSphere_norm (double r, const TVector3 &p)
{
	static bool first = true;
	static TTriangleMesh sph;
	if( first ) { first = false; sph.initAsSphere( 1 , TTRIMESH_SPHERE_RES, TTRIMESH_SPHERE_RES ); }
	
	glPushMatrix();
		glTranslated(  p.data[0],  p.data[1],  p.data[2] );
		glScaled( r,r,r);
		glBegin( GL_TRIANGLES );
		for( int i= 0; i < sph.getPnum(); ++i)
		{
			const TVector3 &v0 = sph.m_verts[ sph.m_polys[i].idx[0]], &n0=sph.m_v_norm[ sph.m_polys[i].idx[0]];
			const TVector3 &v1 = sph.m_verts[ sph.m_polys[i].idx[1]], &n1=sph.m_v_norm[ sph.m_polys[i].idx[1]];
			const TVector3 &v2 = sph.m_verts[ sph.m_polys[i].idx[2]], &n2=sph.m_v_norm[ sph.m_polys[i].idx[2]];
			glNormal3d( n0.data[0] * r, n0.data[1] * r, n0.data[2] * r ); glVertex3dv( v0.data );
			glNormal3d( n1.data[0] * r, n1.data[1] * r, n1.data[2] * r ); glVertex3dv( v1.data );
			glNormal3d( n2.data[0] * r, n2.data[1] * r, n2.data[2] * r ); glVertex3dv( v2.data );
		}	
		glEnd();
	glPopMatrix();
}
//Render sphere objects (center: p, radius:r) with surface normal with specifed color info (ambi/diff/spec/shin)
void t_drawSphere_norm (double r, const TVector3 &p, const float ambi[4], const float diff[4], const float spec[4], const float shin[1])
{
	glMaterialfv( GL_FRONT, GL_SHININESS, shin );
	glMaterialfv( GL_FRONT, GL_SPECULAR,  spec );
	glMaterialfv( GL_FRONT, GL_DIFFUSE,   diff );
	glMaterialfv( GL_FRONT, GL_AMBIENT,   ambi );
	t_drawSphere_norm( r, p );
}
//Render multiple sphere objects (center: Ps, radius:r) with surface normal with specifed color info (ambi/diff/spec/shin)
void t_drawSpheres_norm(double r, const vector<TVector3> &Ps, const float ambi[4], const float diff[4], const float spec[4], const float shin[1])
{
	glMaterialfv( GL_FRONT, GL_SHININESS, shin );
	glMaterialfv( GL_FRONT, GL_SPECULAR,  spec );
	glMaterialfv( GL_FRONT, GL_DIFFUSE,   diff );
	glMaterialfv( GL_FRONT, GL_AMBIENT,   ambi );
	const int size = (int)Ps.size();
	for( int i=0; i<size; ++i) t_drawSphere_norm( r, Ps[i]);
}
//draw/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Marching cubes//////////////////////////////////////////////////////////////////////////////////////


/*------------------------------------------------------------------------------------
 Martching cubes
      
	 4_________5
     /        / |
	/        /  |
 7 /_______6/   |  vertex index 
   |        |   |
   |  0     |  / 1
   |        | /
  3|________|/ 2
    

	  ____4_____
     /        / |
	/7       /5 |9
   /___6___ /   |  edge index 
   |        |   |
 11|     0  |10/ 
   | 3      | /1
   |________|/ 
        2
-------------------------------------------------------------------------------------*/
int TTriangleMesh::MCsEdgeTable[256]={
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

//triangleの貼り方を指定
int TTriangleMesh::MCsTriTable[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
 {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
 {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
 {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
 {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
 {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
 {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
 {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
 {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
 {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
 {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
 {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
 {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
 {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
 {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
 {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
 {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
 {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
 {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
 {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
 {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
 {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
 {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
 {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
 {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
 {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
 {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
 {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
 {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
 {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
 {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
 {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
 {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
 {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
 {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
 {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
 {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
 {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
 {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
 {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
 {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
 {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
 {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
 {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
 {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
 {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
 {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
 {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
 {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
 {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
 {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
 {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
 {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
 {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
 {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
 {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
 { 9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
 {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
 {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
 {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
 { 1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 { 1,  3,  8,  9,  1,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 { 0,  9,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 { 0,  3,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};


/*  sampling cube edge
    y
    |
    |
    |e1                edge id は e0 = 0 残りの9本は隣のvoxelにある
    |                             e1 = 8
    |____e0_____x                 e2 = 3
    \
      \
        \e2
          \
            z
*/
struct SamplVolEdgeIds
{
	int e[3];
	SamplVolEdgeIds()
	{
		e[0]=-1;//x dir
		e[1]=-1;//y
		e[2]=-1;//z
	}
};


//W×H×D のbinary volume (vol)からmarching cubesによりsurfaceを生成する
//volumeを一回り大きいgrid voxelでsamplingする(各voxelの中心が、sampling gridの頂点になる)
void t_MartchingCubes_BinVolue(int W, int H, int D, const byte *vol, const TVector3 &volSize, TTriangleMesh &mesh, const int *minIdx, const int *maxIdx)
{
	int   *edgeTable     =  TTriangleMesh::MCsEdgeTable;
	int  (*triTable)[16] =  TTriangleMesh::MCsTriTable ;

	int preCount = 0;
	for( int i=0; i < W*H*D; ++i) if( vol[i] ) preCount++; 
	const int widhei = W* H;

	const int sv_W = W + 1;//sampling volume resolution
	const int sv_H = H + 1;
	const int sv_D = D + 1;
	//voxelの中心がsampling points
	const double xOff = volSize.data[0] / W ;
	const double yOff = volSize.data[1] / H;
	const double zOff = volSize.data[2] / D ;

	int* svEdgeX = new int[sv_W * sv_H * sv_D];
	int* svEdgeY = new int[sv_W * sv_H * sv_D];
	int* svEdgeZ = new int[sv_W * sv_H * sv_D];
	for( int i=0;i<sv_W * sv_H * sv_D;++i) { svEdgeX[i]= svEdgeY[i]= svEdgeZ[i]= -1;}

	vector< TVector3   > Vs; Vs.reserve( preCount );
	vector< TTriangle  > Ps; Ps.reserve( preCount );


	int xS, xE, yS, yE, zS, zE;
	if( minIdx != 0 && maxIdx != 0){
		xS = max(0, minIdx[0]);  xE = min( maxIdx[0]+2, sv_W); 
		yS = max(0, minIdx[1]);  yE = min( maxIdx[1]+2, sv_H);
		zS = max(0, minIdx[2]);  zE = min( maxIdx[2]+2, sv_D);
	}else{
		xS = 0; yS = 0;  zS= 0; xE = sv_W;  yE = sv_H; zE = sv_D;
	}


	short caseID = 0;
	byte p[8];
	int  v[12];
	for( int z = zS; z < zE; ++z)
	for( int y = yS; y < yE; ++y)
	for( int x = xS; x < xE; ++x) 
	{
		int svIdx  = x + y*sv_W + z*sv_W*sv_H;

		//8点サンプリング
		int pivVox = (x-1) + (y-1) * W + (z-1) * widhei;
		p[0] = (x==0     ||y==0	    ||z==0     ) ? 0 : vol[pivVox				  ];
		p[1] = (x==sv_W-1||y==0	    ||z==0     ) ? 0 : vol[pivVox + 1		  	  ];
		p[2] = (x==sv_W-1||y==0	    ||z==sv_D-1) ? 0 : vol[pivVox + 1     + widhei];
		p[3] = (x==0     ||y==0	    ||z==sv_D-1) ? 0 : vol[pivVox         + widhei];
		p[4] = (x==0     ||y==sv_H-1||z==0     ) ? 0 : vol[pivVox     + W         ];
		p[5] = (x==sv_W-1||y==sv_H-1||z==0     ) ? 0 : vol[pivVox + 1 + W         ];
		p[6] = (x==sv_W-1||y==sv_H-1||z==sv_D-1) ? 0 : vol[pivVox + 1 + W + widhei];
		p[7] = (x==0     ||y==sv_H-1||z==sv_D-1) ? 0 : vol[pivVox +   + W + widhei];

		caseID = 0;
		if( p[0] ) caseID |=   1;
		if( p[1] ) caseID |=   2;
		if( p[2] ) caseID |=   4;
		if( p[3] ) caseID |=   8;
		if( p[4] ) caseID |=  16;
		if( p[5] ) caseID |=  32;
		if( p[6] ) caseID |=  64;
		if( p[7] ) caseID |= 128;

		if (edgeTable[caseID] == 0) continue;

		//以下詳細はノート参照
		int svNexX = svIdx+1        ;
		int svNexY = svIdx+sv_W     ;
		int svNexZ = svIdx+sv_W*sv_H;
		
		int svNexXY = svIdx+1+sv_W          ;
		int svNexYZ = svIdx  +sv_W+sv_W*sv_H;
		int svNexZX = svIdx+1     +sv_W*sv_H;

		double xPos = x*xOff - xOff*0.5;
		double yPos = y*yOff - yOff*0.5;
		double zPos = z*zOff - zOff*0.5;
		if (edgeTable[caseID] &    1/* 0*/){ if( svEdgeX[svIdx  ] ==-1){svEdgeX[svIdx  ] =(int)Vs.size(); Vs.push_back( TVector3(xPos+xOff*0.5, yPos         , zPos          ) );} v[ 0] = svEdgeX[svIdx  ];}
		if (edgeTable[caseID] &    4/* 2*/){ if( svEdgeX[svNexZ ] ==-1){svEdgeX[svNexZ ] =(int)Vs.size(); Vs.push_back( TVector3(xPos+xOff*0.5, yPos         , zPos + zOff   ) );} v[ 2] = svEdgeX[svNexZ ];}
		if (edgeTable[caseID] &   16/* 4*/){ if( svEdgeX[svNexY ] ==-1){svEdgeX[svNexY ] =(int)Vs.size(); Vs.push_back( TVector3(xPos+xOff*0.5, yPos + yOff  , zPos          ) );} v[ 4] = svEdgeX[svNexY ];}
		if (edgeTable[caseID] &   64/* 6*/){ if( svEdgeX[svNexYZ] ==-1){svEdgeX[svNexYZ] =(int)Vs.size(); Vs.push_back( TVector3(xPos+xOff*0.5, yPos + yOff  , zPos + zOff   ) );} v[ 6] = svEdgeX[svNexYZ];}
		if (edgeTable[caseID] &    2/* 1*/){ if( svEdgeZ[svNexX ] ==-1){svEdgeZ[svNexX ] =(int)Vs.size(); Vs.push_back( TVector3(xPos+xOff    , yPos         , zPos +zOff*0.5) );} v[ 1] = svEdgeZ[svNexX ];}
		if (edgeTable[caseID] &    8/* 3*/){ if( svEdgeZ[svIdx  ] ==-1){svEdgeZ[svIdx  ] =(int)Vs.size(); Vs.push_back( TVector3(xPos         , yPos         , zPos +zOff*0.5) );} v[ 3] = svEdgeZ[svIdx  ];}
		if (edgeTable[caseID] &   32/* 5*/){ if( svEdgeZ[svNexXY] ==-1){svEdgeZ[svNexXY] =(int)Vs.size(); Vs.push_back( TVector3(xPos+xOff    , yPos+yOff    , zPos +zOff*0.5) );} v[ 5] = svEdgeZ[svNexXY];}
		if (edgeTable[caseID] &  128/* 7*/){ if( svEdgeZ[svNexY ] ==-1){svEdgeZ[svNexY ] =(int)Vs.size(); Vs.push_back( TVector3(xPos         , yPos+yOff    , zPos +zOff*0.5) );} v[ 7] = svEdgeZ[svNexY ];}
		if (edgeTable[caseID] &  256/* 8*/){ if( svEdgeY[svIdx  ] ==-1){svEdgeY[svIdx  ] =(int)Vs.size(); Vs.push_back( TVector3(xPos         , yPos+yOff*0.5, zPos          ) );} v[ 8] = svEdgeY[svIdx  ];}
		if (edgeTable[caseID] &  512/* 9*/){ if( svEdgeY[svNexX ] ==-1){svEdgeY[svNexX ] =(int)Vs.size(); Vs.push_back( TVector3(xPos+xOff    , yPos+yOff*0.5, zPos          ) );} v[ 9] = svEdgeY[svNexX ];}
		if (edgeTable[caseID] & 1024/*10*/){ if( svEdgeY[svNexZX] ==-1){svEdgeY[svNexZX] =(int)Vs.size(); Vs.push_back( TVector3(xPos+xOff    , yPos+yOff*0.5, zPos+zOff     ) );} v[10] = svEdgeY[svNexZX];}
		if (edgeTable[caseID] & 2048/*11*/){ if( svEdgeY[svNexZ ] ==-1){svEdgeY[svNexZ ] =(int)Vs.size(); Vs.push_back( TVector3(xPos         , yPos+yOff*0.5, zPos+zOff     ) );} v[11] = svEdgeY[svNexZ ];}
		//polygon生成
		for (int i=0; triTable[caseID][i]!=-1;i+=3) 
		{
			Ps.push_back( TTriangle( v[triTable[caseID][ i ]], 
				                     v[triTable[caseID][i+1]], 
								     v[triTable[caseID][i+2]] ) );
		}
	}

	mesh.initFromVertsPolys( Vs, Ps );
	fprintf( stderr, "v= %d poly = %d\n", Vs.size(), Ps.size() );
	delete[] svEdgeX; delete[] svEdgeY; delete[] svEdgeZ;
}

//volumeを一回り大きいgrid voxelでsamplingする(各voxelの中心が、sampling gridの頂点になる)
void t_MartchingCubes_floatVolume( int W, int H, int D,  const float *vol, 
															   float thresh,
													 const TVector3 &cubeSize, 
													    TTriangleMesh &mesh, const int *minIdx, const int *maxIdx)
{
	int   *edgeTable     =  TTriangleMesh::MCsEdgeTable;
	int  (*triTable)[16] =  TTriangleMesh::MCsTriTable ;
	                                                   
	int preCount = 0;
	for( int i=0; i < W*H*D; ++i) if( vol[i] > thresh ) preCount++; 
	const int widhei = W* H;
	
	vector< TVector3 > Vs; Vs.reserve( preCount );
	vector< TTriangle> Ps; Ps.reserve( preCount );
	//sampling volume resolution
	const int sv_W = W +1;
	const int sv_H = H +1;
	const int sv_D = D +1;
	//voxelの中心がsampling points
	const double xOff = cubeSize.data[0] / W ;
	const double yOff = cubeSize.data[1] / H;
	const double zOff = cubeSize.data[2] / D ;
	//実際に利用はしないが，pivotが端に来た時に踏みこえがあり得る
	int* svEdgeX = new int[(sv_W+1) * (sv_H+1) * (sv_D+1)];
	int* svEdgeY = new int[(sv_W+1) * (sv_H+1) * (sv_D+1)];
	int* svEdgeZ = new int[(sv_W+1) * (sv_H+1) * (sv_D+1)];
	for( int i=0; i < (sv_W+1) * (sv_H+1) * (sv_D+1); ++i) { svEdgeX[i]= svEdgeY[i]= svEdgeZ[i]= -1;}
	
	int xS, xE, yS, yE, zS, zE;
	if( minIdx != 0 && maxIdx != 0){
		xS = max(0, minIdx[0]);  xE = min( maxIdx[0]+2, sv_W); 
		yS = max(0, minIdx[1]);  yE = min( maxIdx[1]+2, sv_H);
		zS = max(0, minIdx[2]);  zE = min( maxIdx[2]+2, sv_D);
	}else{
		xS = 0; yS = 0;  zS= 0; xE = sv_W;  yE = sv_H; zE = sv_D;
	}

	short caseID = 0;
	double p[  8 ];
	int    v[ 12 ];                            
                                                                              	                                                                                     
	for( int z = zS; z < zE; ++z)
	for( int y = yS; y < yE; ++y)
	for( int x = xS; x < xE; ++x) 
	{
		int svIdx  = x + y*sv_W + z*sv_W*sv_H;

		//8点サンプリング
		int pivVox = (x-1) + (y-1) * W + (z-1) * widhei;
		p[0] = (x==0     ||y==0	    ||z==0     ) ? -1 : vol[pivVox				  ];
		p[1] = (x==sv_W-1||y==0	    ||z==0     ) ? -1 : vol[pivVox + 1		      ];
		p[2] = (x==sv_W-1||y==0	    ||z==sv_D-1) ? -1 : vol[pivVox + 1     + widhei];
		p[3] = (x==0     ||y==0	    ||z==sv_D-1) ? -1 : vol[pivVox         + widhei];
		p[4] = (x==0     ||y==sv_H-1||z==0     ) ? -1 : vol[pivVox     + W         ];
		p[5] = (x==sv_W-1||y==sv_H-1||z==0     ) ? -1 : vol[pivVox + 1 + W         ];
		p[6] = (x==sv_W-1||y==sv_H-1||z==sv_D-1) ? -1 : vol[pivVox + 1 + W + widhei];
		p[7] = (x==0     ||y==sv_H-1||z==sv_D-1) ? -1 : vol[pivVox +   + W + widhei];

		caseID = 0;
		if( p[0] > thresh ) caseID |=   1;
		if( p[1] > thresh ) caseID |=   2;
		if( p[2] > thresh ) caseID |=   4;
		if( p[3] > thresh ) caseID |=   8;
		if( p[4] > thresh ) caseID |=  16;
		if( p[5] > thresh ) caseID |=  32;
		if( p[6] > thresh ) caseID |=  64;
		if( p[7] > thresh ) caseID |= 128;

		if (edgeTable[caseID] == 0) continue;

		//以下詳細はノート参照
		int svNexX = svIdx+1        ;
		int svNexY = svIdx+sv_W     ;
		int svNexZ = svIdx+sv_W*sv_H;
		
		int svNexXY = svIdx+1+sv_W          ;
		int svNexYZ = svIdx  +sv_W+sv_W*sv_H;
		int svNexZX = svIdx+1     +sv_W*sv_H;

		double xPos = x*xOff - xOff*0.5;  double xNext = xPos+xOff;
		double yPos = y*yOff - yOff*0.5;  double yNext = yPos+yOff;
		double zPos = z*zOff - zOff*0.5;  double zNext = zPos+zOff;

		/*1 :0   p0 - p1*/
		/*4 :2   p3 - p2*/
		/*16:4   p4 - p5*/
		/*64:6   p7 - p6*/
		if(edgeTable[caseID]&   1 && svEdgeX[svIdx  ]==-1){svEdgeX[svIdx  ]=(int)Vs.size(); Vs.push_back( TVector3(xPos + xOff*(thresh-p[0])/(p[1]-p[0]), yPos , zPos ) );}
		if(edgeTable[caseID]&   4 && svEdgeX[svNexZ ]==-1){svEdgeX[svNexZ ]=(int)Vs.size(); Vs.push_back( TVector3(xPos + xOff*(thresh-p[3])/(p[2]-p[3]), yPos , zNext) );}
		if(edgeTable[caseID]&  16 && svEdgeX[svNexY ]==-1){svEdgeX[svNexY ]=(int)Vs.size(); Vs.push_back( TVector3(xPos + xOff*(thresh-p[4])/(p[5]-p[4]), yNext, zPos ) );}		
		if(edgeTable[caseID]&  64 && svEdgeX[svNexYZ]==-1){svEdgeX[svNexYZ]=(int)Vs.size(); Vs.push_back( TVector3(xPos + xOff*(thresh-p[7])/(p[6]-p[7]), yNext, zNext) );}		
		/*  2:1  p1-p2*/
		/*  8:3  p0-p3*/
		/* 32:5  p5-p6*/
		/*128:7  p4-p7*/
		if(edgeTable[caseID]&   2 && svEdgeZ[svNexX ]==-1){svEdgeZ[svNexX ]=(int)Vs.size(); Vs.push_back( TVector3(xNext, yPos , zPos +zOff*(thresh-p[1])/(p[2]-p[1])) );}		
		if(edgeTable[caseID]&   8 && svEdgeZ[svIdx  ]==-1){svEdgeZ[svIdx  ]=(int)Vs.size(); Vs.push_back( TVector3(xPos , yPos , zPos +zOff*(thresh-p[0])/(p[3]-p[0])) );}		
		if(edgeTable[caseID]&  32 && svEdgeZ[svNexXY]==-1){svEdgeZ[svNexXY]=(int)Vs.size(); Vs.push_back( TVector3(xNext, yNext, zPos +zOff*(thresh-p[5])/(p[6]-p[5])) );}			
		if(edgeTable[caseID]& 128 && svEdgeZ[svNexY ]==-1){svEdgeZ[svNexY ]=(int)Vs.size(); Vs.push_back( TVector3(xPos , yNext, zPos +zOff*(thresh-p[4])/(p[7]-p[4])) );}		
		/* 256:8   p0-p4*/
		/* 512:9   p1-p5*/
		/*1024:10  p2-p6*/
		/*2048:11  p3-p7*/
		if(edgeTable[caseID]& 256 && svEdgeY[svIdx  ]==-1){svEdgeY[svIdx  ]=(int)Vs.size(); Vs.push_back( TVector3(xPos , yPos+yOff*(thresh-p[0])/(p[4]-p[0]), zPos ) );}		
		if(edgeTable[caseID]& 512 && svEdgeY[svNexX ]==-1){svEdgeY[svNexX ]=(int)Vs.size(); Vs.push_back( TVector3(xNext, yPos+yOff*(thresh-p[1])/(p[5]-p[1]), zPos ) );}		
		if(edgeTable[caseID]&1024 && svEdgeY[svNexZX]==-1){svEdgeY[svNexZX]=(int)Vs.size(); Vs.push_back( TVector3(xNext, yPos+yOff*(thresh-p[2])/(p[6]-p[2]), zNext) );}		
		if(edgeTable[caseID]&2048 && svEdgeY[svNexZ ]==-1){svEdgeY[svNexZ ]=(int)Vs.size(); Vs.push_back( TVector3(xPos , yPos+yOff*(thresh-p[3])/(p[7]-p[3]), zNext) );}		
	
		v[0] = svEdgeX[svIdx  ];
		v[2] = svEdgeX[svNexZ ];
		v[4] = svEdgeX[svNexY ];
		v[6] = svEdgeX[svNexYZ];

		v[1] = svEdgeZ[svNexX ];
		v[3] = svEdgeZ[svIdx  ];
		v[5] = svEdgeZ[svNexXY];
		v[7] = svEdgeZ[svNexY ];

		v[8] = svEdgeY[svIdx  ];
		v[9] = svEdgeY[svNexX ];	
		v[10]= svEdgeY[svNexZX];
		v[11]= svEdgeY[svNexZ ];

		//polygon生成
		for (int i=0; triTable[caseID][i]!=-1;i+=3) 
		{
			Ps.push_back( TTriangle( v[triTable[caseID][ i ]], 
				                     v[triTable[caseID][i+1]], 
									 v[triTable[caseID][i+2]] ) );
		}
	}
	mesh.initFromVertsPolys( Vs, Ps );

	fprintf( stderr, "v= %d poly = %d\n", Vs.size(), Ps.size() );
	delete[] svEdgeX; 
	delete[] svEdgeY; 
	delete[] svEdgeZ;
}



//binary volumeからvoxel生成
void t_GenMeshFromBinaryVolume(int W, int H, int D, const byte *vol, const TVector3 &volSize, TTriangleMesh &mesh)
{
	double xOff = volSize.data[0] / W ;
	double yOff = volSize.data[1] / H;
	double zOff = volSize.data[2] / D ;

	int *voxId2vtxId = new int[ (W+1) * (H+1) * (D+1) ];
	for(int i=0; i<(W+1) * (H+1) * (D+1) ;++i) voxId2vtxId[i] = -1;

	int preCount = 0;
	int widhei = W* H;
	for( int z = 0; z < D ; ++z)
	for( int y = 0; y < H ; ++y)
	for( int x = 0; x < W ; ++x) 
	{
		int idx = x + y * W + z * widhei; if( !vol[idx] )continue;
		if( x == 0   || !vol[idx-1     ] ) ++preCount;
		if( x == W-1 || !vol[idx+1     ] ) ++preCount;
		if( y == 0   || !vol[idx-W     ] ) ++preCount;
		if( y == H-1 || !vol[idx+W     ] ) ++preCount;
		if( z == 0   || !vol[idx-widhei] ) ++preCount;
		if( z == 0   || !vol[idx+widhei] ) ++preCount;
	}

	fprintf( stderr, "polygon generation preCount %d\n", preCount);

	vector< TVector3 > Vs; Vs.reserve( preCount   );
	vector< TTriangle > Ps; Ps.reserve( preCount*2 );

	TVector3 p;
	for( int z = 0; z < D ; ++z)
	for( int y = 0; y < H; ++y)
	for( int x = 0; x < W ; ++x) 
	{
		int idx = x + y * W + z * widhei;	
		if( !vol[idx] ) continue;

		int p0id = idx          , p4id = idx         + widhei;
		int p1id = idx + 1      , p5id = idx + 1     + widhei;	
		int p2id = idx     + W  , p6id = idx     + W + widhei;	
		int p3id = idx + 1 + W  , p7id = idx + 1 + W + widhei;
		

		if( x == 0 || !vol[idx-1] ) //-xに面を貼る
		{
			if( voxId2vtxId[p0id]==-1){ voxId2vtxId[p0id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*  y  , zOff*  z  ); Vs.push_back( p); }
			if( voxId2vtxId[p2id]==-1){ voxId2vtxId[p2id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*(y+1), zOff*  z  ); Vs.push_back( p); }
			if( voxId2vtxId[p4id]==-1){ voxId2vtxId[p4id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*  y  , zOff*(z+1)); Vs.push_back( p); }
			if( voxId2vtxId[p6id]==-1){ voxId2vtxId[p6id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*(y+1), zOff*(z+1)); Vs.push_back( p); }
			Ps.push_back( TTriangle( voxId2vtxId[p0id], voxId2vtxId[p4id], voxId2vtxId[p2id]) );
			Ps.push_back( TTriangle( voxId2vtxId[p2id], voxId2vtxId[p4id], voxId2vtxId[p6id]) );
		}
		
		if( x == W-1 || !vol[idx+1]) //+xに面を貼る
		{
			if( voxId2vtxId[p1id]==-1){ voxId2vtxId[p1id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*  y  , zOff*  z  ); Vs.push_back( p); }
			if( voxId2vtxId[p3id]==-1){ voxId2vtxId[p3id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*(y+1), zOff*  z  ); Vs.push_back( p); }
			if( voxId2vtxId[p5id]==-1){ voxId2vtxId[p5id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*  y  , zOff*(z+1)); Vs.push_back( p); }
			if( voxId2vtxId[p7id]==-1){ voxId2vtxId[p7id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*(y+1), zOff*(z+1)); Vs.push_back( p); }
			Ps.push_back( TTriangle( voxId2vtxId[p1id], voxId2vtxId[p7id], voxId2vtxId[p5id]) );
			Ps.push_back( TTriangle( voxId2vtxId[p1id], voxId2vtxId[p3id], voxId2vtxId[p7id]) );
		}

		if( y == 0       || !vol[idx-W]) //-yに面を貼る
		{
			if( voxId2vtxId[p0id]==-1){ voxId2vtxId[p0id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*  y  , zOff*  z  ); Vs.push_back( p); }
			if( voxId2vtxId[p1id]==-1){ voxId2vtxId[p1id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*  y  , zOff*  z  ); Vs.push_back( p); }
			if( voxId2vtxId[p4id]==-1){ voxId2vtxId[p4id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*  y  , zOff*(z+1)); Vs.push_back( p); }
			if( voxId2vtxId[p5id]==-1){ voxId2vtxId[p5id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*  y  , zOff*(z+1)); Vs.push_back( p); }
			Ps.push_back( TTriangle( voxId2vtxId[p0id], voxId2vtxId[p1id], voxId2vtxId[p5id]) );
			Ps.push_back( TTriangle( voxId2vtxId[p0id], voxId2vtxId[p5id], voxId2vtxId[p4id]) );
		}
		
		if( y == H-1|| !vol[idx+W] ) //+yに面を貼る
		{
			if( voxId2vtxId[p2id]==-1){ voxId2vtxId[p2id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*(y+1), zOff*  z  ); Vs.push_back( p); }
			if( voxId2vtxId[p3id]==-1){ voxId2vtxId[p3id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*(y+1), zOff*  z  ); Vs.push_back( p); }
			if( voxId2vtxId[p6id]==-1){ voxId2vtxId[p6id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*(y+1), zOff*(z+1)); Vs.push_back( p); }
			if( voxId2vtxId[p7id]==-1){ voxId2vtxId[p7id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*(y+1), zOff*(z+1)); Vs.push_back( p); }
			Ps.push_back( TTriangle( voxId2vtxId[p2id], voxId2vtxId[p6id], voxId2vtxId[p3id]) );
			Ps.push_back( TTriangle( voxId2vtxId[p3id], voxId2vtxId[p6id], voxId2vtxId[p7id]) );
		}

		if( z == 0       || !vol[idx-widhei] ) //-zに面を貼る
		{			
			if( voxId2vtxId[p0id]==-1){ voxId2vtxId[p0id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*  y  , zOff*z); Vs.push_back( p); }
			if( voxId2vtxId[p1id]==-1){ voxId2vtxId[p1id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*  y  , zOff*z); Vs.push_back( p); }
			if( voxId2vtxId[p2id]==-1){ voxId2vtxId[p2id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*(y+1), zOff*z); Vs.push_back( p); }
			if( voxId2vtxId[p3id]==-1){ voxId2vtxId[p3id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*(y+1), zOff*z); Vs.push_back( p); }
			Ps.push_back( TTriangle( voxId2vtxId[p0id], voxId2vtxId[p3id], voxId2vtxId[p1id]) );
			Ps.push_back( TTriangle( voxId2vtxId[p0id], voxId2vtxId[p2id], voxId2vtxId[p3id]) );
		}
		
		if( z == D-1 || !vol[idx+widhei] ) //+zに面を貼る
		{	
			if( voxId2vtxId[p4id]==-1){ voxId2vtxId[p4id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*  y  , zOff*(z+1)); Vs.push_back( p); }
			if( voxId2vtxId[p5id]==-1){ voxId2vtxId[p5id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*  y  , zOff*(z+1)); Vs.push_back( p); }
			if( voxId2vtxId[p6id]==-1){ voxId2vtxId[p6id] = (int)Vs.size(); p.Set(xOff*  x  , yOff*(y+1), zOff*(z+1)); Vs.push_back( p); }
			if( voxId2vtxId[p7id]==-1){ voxId2vtxId[p7id] = (int)Vs.size(); p.Set(xOff*(x+1), yOff*(y+1), zOff*(z+1)); Vs.push_back( p); }
			Ps.push_back( TTriangle( voxId2vtxId[p4id], voxId2vtxId[p5id], voxId2vtxId[p6id]) );
			Ps.push_back( TTriangle( voxId2vtxId[p5id], voxId2vtxId[p7id], voxId2vtxId[p6id]) );
		}
	}
	mesh.initFromVertsPolys( Vs,Ps );
	delete[] voxId2vtxId;
}


void t_GenMeshFromBinaryVolume_coarseResNew( int W, int H, int D, 
											const byte        *vol_mask,
											const int          trgtMaskId , 
											const TVector3  &volSize, 
											TTriangleMesh &poly, int resolution)
{
	const int    r  = resolution;
	const int    WH = W*H;
	const double r3 = r*r*r;

	//volume coarse resolution
	int rW = W  / r,  rW_over = W % r;
	int rH = H  / r,  rH_over = H % r;
	int rD = D  / r,  rD_over = D % r;
	if(rW_over!=0) rW ++;
	if(rH_over!=0) rH ++;
	if(rD_over!=0) rD ++;

	double count, countAll;
	
	byte *vol = new byte[rW*rH*rD];
	for( int z = 0; z < rD; ++z)
	for( int y = 0; y < rH; ++y)
	for( int x = 0; x < rW; ++x) 
	{
		int volId  =  x + y * rW + z * rW * rH   ;
		int idx    = (x + y *  W + z * WH   )*r;

		count = countAll = 0;
		for( int zzz = 0; zzz < r; ++zzz ) if( z*r + zzz < D )
		for( int yyy = 0; yyy < r; ++yyy ) if( y*r + yyy < H )
		for( int xxx = 0; xxx < r; ++xxx ) if( x*r + xxx < W )
		{
			countAll += 1;
			if( vol_mask[idx+xxx+yyy*W + zzz*WH] == trgtMaskId ) count+=1;
		}
		vol[volId] = ( count/countAll > 0.501)? 255: 0 ;
	}
	TVector3 coarseVSize;
	coarseVSize.data[0] = ( rW_over == 0 ) ? volSize.data[0] : volSize.data[0]/W * r*rW;//volSize.data[0]/widthというvoxelが r*rW個あるのでこのサイズ
	coarseVSize.data[1] = ( rH_over == 0 ) ? volSize.data[1] : volSize.data[1]/H * r*rH;
	coarseVSize.data[2] = ( rD_over == 0 ) ? volSize.data[2] : volSize.data[2]/D * r*rD;
	t_MartchingCubes_BinVolue( rW, rH, rD, vol, coarseVSize, poly);

	delete[] vol;
	return;
}

void t_GenMeshFromBinaryVolume_coarseResNew( int W, int H, int D, const byte *vol_mask, const vector<int> &maskOnOff, const TVector3 &volSize, 
										 TTriangleMesh &mesh, int resolution)
{
	const int    r   = resolution;
	const int    WH  = W * H;
	const double r3  = r*r*r;

	//volume coarse resolution
	int rW = W / r,  rW_over = W % r;
	int rH = H / r,  rH_over = H % r;
	int rD = D / r,  rD_over = D % r;
	if(rW_over!=0) rW ++;
	if(rH_over!=0) rH ++;
	if(rD_over!=0) rD ++;

	double count, countAll;
	
	byte *vol = new byte[rW*rH*rD];
	for( int z = 0; z < rD; ++z)
	for( int y = 0; y < rH; ++y)
	for( int x = 0; x < rW; ++x) 
	{
		int volId  =  x + y * rW + z * rW * rH    ;
		int idx    = (x + y * W  + z * WH      )*r;

		count = countAll = 0;
		for( int zzz = 0; zzz < r; ++zzz ) if( z*r + zzz < D )
		for( int yyy = 0; yyy < r; ++yyy ) if( y*r + yyy < H )
		for( int xxx = 0; xxx < r; ++xxx ) if( x*r + xxx < W )
		{
			countAll += 1;
			if( maskOnOff[ vol_mask[idx+xxx+yyy*W+zzz*WH] ] ) count+=1;
		}

		vol[volId] = ( count/countAll > 0.501)? 255: 0 ;
	}
	TVector3 coarseVSize;
	coarseVSize.data[0] = ( rW_over == 0 ) ? volSize.data[0] : volSize.data[0]/W * r*rW;//volSize.data[0]/widthというvoxelが r*rW個あるのでこのサイズ
	coarseVSize.data[1] = ( rH_over == 0 ) ? volSize.data[1] : volSize.data[1]/H * r*rH;
	coarseVSize.data[2] = ( rD_over == 0 ) ? volSize.data[2] : volSize.data[2]/D * r*rD;
	t_MartchingCubes_BinVolue( rW, rH, rD, vol, coarseVSize, mesh);

	delete[] vol;
	return;
}

void t_GenMeshFromBinaryVolume_coarseResNew( int W, int H, int D, byte *volInOut, const TVector3 &volSize, TTriangleMesh &mesh, int resolution )
{
	const int    r  = resolution;
	const int    WH = W*H;
	const double r3 = r*r*r;

	//volume coarse resolution
	int rW = W  / r,  rW_over = W % r;
	int rH = H  / r,  rH_over = H % r;
	int rD = D  / r,  rD_over = D % r;
	if(rW_over!=0) rW ++;
	if(rH_over!=0) rH ++;
	if(rD_over!=0) rD ++;


	double count, countAll;
	
	byte *vol = new byte[rW*rH*rD];
	for( int z = 0; z < rD; ++z)
	for( int y = 0; y < rH; ++y)
	for( int x = 0; x < rW; ++x) 
	{
		int volId  =  x + y * rW    + z * rW * rH   ;
		int idx    = (x + y * W + z*WH   )*r;

		count = countAll = 0;
		for( int zzz = 0; zzz < r; ++zzz ) if( z*r + zzz < D )
		for( int yyy = 0; yyy < r; ++yyy ) if( y*r + yyy < H)
		for( int xxx = 0; xxx < r; ++xxx ) if( x*r + xxx < W )
		{
			countAll += 1;
			if( volInOut[ idx+xxx+yyy*W+zzz*WH ] ) count+=1;
		}

		vol[volId] = (count/countAll > 0.501) ? 255 : 0 ;
	}
	TVector3 coarseVSize;
	coarseVSize.data[0] = ( rW_over == 0 ) ? volSize.data[0] : volSize.data[0]/W * r*rW;//volSize.data[0]/widthというvoxelが r*rW個あるのでこのサイズ
	coarseVSize.data[1] = ( rH_over == 0 ) ? volSize.data[1] : volSize.data[1]/H* r*rH;
	coarseVSize.data[2] = ( rD_over == 0 ) ? volSize.data[2] : volSize.data[2]/D * r*rD;
	t_MartchingCubes_BinVolue( rW, rH, rD, vol, coarseVSize, mesh);

	delete[] vol;
	return;
}



///////////////////////////////////////////////////////////////////////////////////////////
void t_computeInternalVoxel( const int width, const int height, const int depth, TTriangleMesh &trgtObj, 
							 const TVector3     &volumeSize,  byte           *bInside )
{
	const int widhei = width * height;
	double voxelWidthX = volumeSize.data[0] / width ;
	double voxelWidthY = volumeSize.data[1] / height;
	double voxelWidthZ = volumeSize.data[2] / depth ;
	double offsetX     = 0.5 * voxelWidthX;
	double offsetY     = 0.5 * voxelWidthY;
	double offsetZ     = 0.5 * voxelWidthZ;

	memset( bInside, 0, sizeof( byte ) * width * height * depth );

	trgtObj.updateNormal();

	TVector3 BBmin, BBmax; 
	trgtObj.ComputeBBox( BBmin, BBmax );

	const TVector3 *verts = trgtObj.m_verts;
	const TTriangle  *polys = trgtObj.m_polys;

	/////////////////////////bin作成///////////////////////////////////////
	int binSize = 10;
	vector< vector<int> > polyIdBins( binSize*binSize, vector<int>() );
	vector< double      > Zpiv(binSize+1);
	vector< double      > Ypiv(binSize+1);

	for( int i=0; i <  binSize*binSize;++i) polyIdBins[i].reserve((int) ( trgtObj.getPnum()*0.5) );
	for( int i=0; i <= binSize        ;++i){
		Ypiv[i] = volumeSize.data[1] / binSize * i;
		Zpiv[i] = volumeSize.data[2] / binSize * i;
	}

	for( int i=0; i<trgtObj.getPnum(); ++i)
	{
		const TVector3 &V0 = verts[ polys[i].idx[0] ];
		const TVector3 &V1 = verts[ polys[i].idx[1] ];
		const TVector3 &V2 = verts[ polys[i].idx[2] ];

		for( int zi = 0; zi < binSize; ++zi)
		for( int yi = 0; yi < binSize; ++yi)
		{
			if( ( (Zpiv[zi]<=V0.data[2] && V0.data[2]<=Zpiv[zi+1]) || (Zpiv[zi]<=V1.data[2] && V1.data[2]<=Zpiv[zi+1]) || (Zpiv[zi]<=V2.data[2] && V2.data[2]<=Zpiv[zi+1]) ) && 
			    ( (Ypiv[yi]<=V0.data[1] && V0.data[1]<=Ypiv[yi+1]) || (Ypiv[yi]<=V1.data[1] && V1.data[1]<=Ypiv[yi+1]) || (Ypiv[yi]<=V2.data[1] && V2.data[1]<=Ypiv[yi+1]) ) )
			{
				polyIdBins[ zi * binSize + yi ].push_back(i);
			}
		}
	}

	for (int z_idx = 0; z_idx < depth ; ++z_idx)
	{
		double z = z_idx * voxelWidthZ + offsetZ;
		if( z < BBmin.data[2] || BBmax.data[2] < z ) continue;
		
		for (int y_idx = 0; y_idx < height; ++y_idx) 
		{
			double y = y_idx * voxelWidthY + offsetY;
			if( y < BBmin.data[1] || BBmax.data[1] < y ) continue;

			multimap<double, double> blist;// xPos & in or out
			vector< int > *trgtPolyBin = 0;

			for( int zi = 0; zi < binSize; ++zi)
			for( int yi = 0; yi < binSize; ++yi) 
				if( (Zpiv[zi] <= z && z <= Zpiv[zi+1]) && (Ypiv[yi] <= y && y <= Ypiv[yi+1]) )
				{
					trgtPolyBin = &polyIdBins[ zi * binSize + yi ];
					break;
				}
				
			if( trgtPolyBin == 0 ){ fprintf( stderr, "ray is out of volume\n"); continue;}

			for(int pi = 0; pi < (int) trgtPolyBin->size(); ++pi )
			{
				const int         *pVtx = polys[ (*trgtPolyBin)[pi] ].idx;
				const TVector3 &pNorm = trgtObj.m_p_norm[(*trgtPolyBin)[pi]];
				const TVector3 &V0 = verts[ pVtx[0] ];
				const TVector3 &V1 = verts[ pVtx[1] ];
				const TVector3 &V2 = verts[ pVtx[2] ];

				double x0 = V0.data[0], y0 = V0.data[1], z0 = V0.data[2];
				double x1 = V1.data[0], y1 = V1.data[1], z1 = V1.data[2];
				double x2 = V2.data[0], y2 = V2.data[1], z2 = V2.data[2];

				if( (y < y0 && y < y1 && y < y2) || (y > y0 && y > y1 && y > y2) ) continue;//pre-check
				if( (z < z0 && z < z1 && z < z2) || (z > z0 && z > z1 && z > z2) ) continue;
				
				double s,t;
				if( !t_solve2by2LinearEquation(y1-y0, y2-y0,z1-z0, z2-z0, y-y0, z-z0, s, t) ) continue; 
				
				if (s < 0 || t < 0 || s+t > 1) continue;
				
				double x = (1-s-t)*x0 + s*x1 + t*x2;//x 座標
				if( pNorm.data[0] != 0 ) blist.insert( pair<double,double>( x, pNorm.data[0]) );
			}

			if( blist.size() == 0 ) continue;

			//edgeにrayが交差すると同じ所に2点(or 3点)交点ができてしまう
			while( blist.size() != 0 )
			{
				if( blist.size() == 1 ){ blist.clear(); break;}

				bool found = false;
				multimap<double,double>::iterator it0 = blist.begin();
				multimap<double,double>::iterator it1 = blist.begin(); it1++;
				
				for(; it1 != blist.end(); ++it0, ++it1)
					if( it0->second * it1->second > 0)
					{
						blist.erase( it1 );
						found = true;
						break;
					}
				if( !found ) break;
			}				
		

			bool flag = false;
			int x_idx = 0;
			int pivIdx = y_idx * width + z_idx * widhei;
			for (multimap<double,double>::iterator it = blist.begin(); it != blist.end(); ++it) 
			{
				//rayと面との交点のちょうど手前のvoxel idx
				int pivIdx_pre = (int)( (it->first - offsetX) / voxelWidthX);
				for( ; x_idx <= pivIdx_pre && x_idx < width; ++x_idx ) bInside[ x_idx + pivIdx ] = flag; 

				flag = !flag;
			}

			if( flag == true)
			{
				// - - ++ や + + - - などmeshがひっくり返るとおかしなことが起こる
				fprintf( stderr, "fail      strange mesh shape !!!!!!!");
				flag = false;
			}
			//残りを塗る
			for ( ; x_idx < width; ++x_idx) bInside[ x_idx + pivIdx ] = flag;
		}
	}
}



void t_createCylinderSolid( double radius, double sizeY, int vertSizeR, int vertSizeY, TTriangleMesh &mesh)
{
	vector< TTriangle  > Ps;
	vector< TVector3 > Vs;
	
	double stepY = sizeY    / ( (double)vertSizeY - 1.0 );
	double stepR = 2 * M_PI / ( (double)vertSizeR - 1.0);

	for( int y = 0; y < vertSizeY; ++y)
	for( int r = 0; r < vertSizeR; ++r)
	{
		Vs.push_back( TVector3(  radius * cos( stepR * r ), stepY * y, 
			                         -radius * sin( stepR * r )) );

		if( y != 0 )
		{
//			//piv = vertSizeR * y + r
//        piv-1 +-----+ piv                                    piv  + vertSize - 1  +-----+ piv
//			    |    /|                                                             |    /|  
//			    |   / |                                                             |   / |  
//			    |  /  |                                                             |  /  |  
//			    | /   |                                                             | /   |  
//              +-----+ piv - vertSizeR                                             +-----+ piv - vertSizeR
// piv - vertSizeR - 1                                                            piv  - 1
			int piv = vertSizeR * y + r;
			if( r == 0 ){ Ps.push_back( TTriangle(piv, piv + vertSizeR -1, piv - 1             ) );
				          Ps.push_back( TTriangle(piv, piv             -1, piv - vertSizeR     ) );	}
			else{         Ps.push_back( TTriangle(piv, piv            - 1, piv - vertSizeR - 1 ) );
				          Ps.push_back( TTriangle(piv, piv - vertSizeR -1, piv - vertSizeR     ) ); }
		}
	}

	//ふた下
	Vs.push_back( TVector3( 0,0,0 ));
	int piv = (int) Vs.size()-1;
	for( int i = 0; i < vertSizeR ; ++i){
		if( i == 0 ) Ps.push_back( TTriangle(piv, i, vertSizeR-1) );
		else         Ps.push_back( TTriangle(piv, i, i-1      ) );
	}

	//ふた上
	Vs.push_back( TVector3( 0,sizeY,0 ));
	piv = (int) Vs.size()-1;
	for( int i = 0; i < vertSizeR ; ++i)
	{
		int k = vertSizeR * ( vertSizeY - 1);
		if( i == 0 ) Ps.push_back( TTriangle(piv, k + vertSizeR - 1,  k + i) );
		else         Ps.push_back( TTriangle(piv, k + i - 1    ,  k + i) );
	}
	mesh.initFromVertsPolys(Vs,Ps);
}

//center = 0,0,0
//r      = gRadの円(x-y平面)上に 
//generalized cylinderを生成する
void t_createRotatedCylinderSolid(double gRad, double radius,int vertSizeR, int vertSizeY, TTriangleMesh &mesh)
{
	const double stepR = 2 * M_PI / ( (double)vertSizeR - 1.0);
	vector< TTriangle  > Ps;
	vector< TVector3 > Vs;

	TMatrix16 frame;
	TVector3 v;

	for( int y = 0; y < vertSizeY; ++y)
	{
		double rot = y * 2 * M_PI / ((double) vertSizeY + 5);
		frame.SetIdentity();
		frame.SetAsRotateZ( -rot );
		for( int r = 0; r < vertSizeR; ++r)
		{
			v.Set( radius * cos( stepR * r ), 0, -radius * sin( stepR * r ) );
			v.data[0] -= gRad;
			t_MatMultVec( v, frame );
			Vs.push_back( v);

			if( y != 0 )
			{
				int piv = vertSizeR * y + r;
				if( r == 0 ){ Ps.push_back( TTriangle(piv, piv + vertSizeR -1, piv - 1             ) );
					          Ps.push_back( TTriangle(piv, piv             -1, piv - vertSizeR     ) );	}
				else {        Ps.push_back( TTriangle(piv, piv            - 1, piv - vertSizeR - 1 ) );
					          Ps.push_back( TTriangle(piv, piv - vertSizeR -1, piv - vertSizeR     ) );}
			}
		}
	}
	//ふたした
	Vs.push_back( TVector3( -gRad, 0, 0 ));
	int piv = (int) Vs.size()-1;
	
	for( int i = 0; i < vertSizeR ; ++i)
	{
		if( i == 0 ) Ps.push_back( TTriangle(piv, i, vertSizeR-1) );
		else         Ps.push_back( TTriangle(piv, i, i-1      ) );
	}


	//ふたうえ
	Vs.push_back( frame * TVector3( -gRad, 0, 0 ));
	piv = (int) Vs.size()-1;
	for( int i = 0; i < vertSizeR ; ++i)
	{
		int k = vertSizeR * ( vertSizeY - 1);
		if( i == 0 ) Ps.push_back( TTriangle(piv, k + vertSizeR - 1,  k + i) );
		else         Ps.push_back( TTriangle(piv, k + i - 1        ,  k + i) );
	}
	mesh.initFromVertsPolys( Vs, Ps );

}


//Todo vtracerとマージする時に 書きなおす必要があると思う
bool t_computeEdgeStrip( const int eSize, const TWingEdge *edges,
						 const int pSize, const TTriangle *polys,
						                  const TTriangle *polyEdges  ,
						 const vector<int>       &e_bCrossing,
							vector<vector<int>>  &e_strips)
{
	vector<int> e_flag(e_bCrossing.size(), 0);
	for( int i=0;i<(int)e_bCrossing.size(); ++i )if( e_bCrossing[i] ) e_flag[i] = 1;

	for( int i=0; i<pSize;++i) if( e_flag[ polyEdges[i].idx[0] ] && 
		                           e_flag[ polyEdges[i].idx[1] ] && 
								   e_flag[ polyEdges[i].idx[2] ] ) return false;

	while(true)
	{	
		int pivEdgeId = -1;
		int prePoly   = -1;
		for( int i=0; i<pSize;++i)
		{
			bool e0 = e_flag[polyEdges[i].idx[0]] == 1;
			bool e1 = e_flag[polyEdges[i].idx[1]] == 1;
			bool e2 = e_flag[polyEdges[i].idx[2]] == 1;
			if(  e0 && !e1 && !e2 ) {pivEdgeId = polyEdges[i].idx[0]; prePoly = i; break;}
			if( !e0 &&  e1 && !e2 ) {pivEdgeId = polyEdges[i].idx[1]; prePoly = i; break;}
			if( !e0 && !e1 &&  e2 ) {pivEdgeId = polyEdges[i].idx[2]; prePoly = i; break;}
		}

		if( pivEdgeId == -1 )
		{
			for( int i=0; i < eSize;++i) if( e_flag[i] == 1) 
			{
				pivEdgeId = i; 
				prePoly = edges[pivEdgeId].p[0];
				break;
			}			
		}

		if( pivEdgeId == -1 ) break;

		e_strips.push_back( vector<int>() );
		e_strips.back().push_back( pivEdgeId );
		e_flag[pivEdgeId] = 2;

		while( true )
		{
			int pivPoly = (prePoly != edges[ pivEdgeId ].p[0]) ? edges[ pivEdgeId ].p[0]: 
				                                                 edges[ pivEdgeId ].p[1] ;
			const TTriangle &pEdges = polyEdges[ pivPoly ];
			if     ( e_flag[ pEdges.idx[0] ] == 1 ) pivEdgeId = pEdges.idx[0];
			else if( e_flag[ pEdges.idx[1] ] == 1 ) pivEdgeId = pEdges.idx[1];
			else if( e_flag[ pEdges.idx[2] ] == 1 ) pivEdgeId = pEdges.idx[2];
			else                                          break;//fiberが途中で切れる
			
			e_strips.back().push_back( pivEdgeId );
			e_flag[pivEdgeId] = 2;
			prePoly = pivPoly;
		}
	}
	return true;
}






















