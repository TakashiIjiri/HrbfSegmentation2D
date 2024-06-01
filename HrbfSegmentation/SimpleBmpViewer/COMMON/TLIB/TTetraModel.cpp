#include "StdAfx.h"
#include "TTetraModel.h"

#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <algorithm>

const int TTetraModel::FACE_VTX_ORDER[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

TTetraModel::TTetraModel(void)
{
}

TTetraModel::~TTetraModel(void)
{
}

void TTetraModel::print() {
	using namespace std;
	fprintf( stderr,  "vertices: \n");
	for(int i=0; i<(int)m_verts .size(); ++i) fprintf( stderr, "%f  %f  %f   \n", m_verts[i].data[0], m_verts[i].data[1], m_verts[i].data[2] );
	for(int i=0; i<(int)m_tetras.size(); ++i) fprintf( stderr, "%d  %d  %d %d\n", m_tetras[i].v  [0], m_tetras[i].v  [1], m_tetras[i].v  [2], m_tetras[i].v[3]   );
}


void TTetraModel::calcNormal() 
{
	static int order[][3] = {{3, 2, 1}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};
	
	for (int i = 0; i < (int)m_tetras.size(); ++i) 
	{
		TTetra& tet = m_tetras[i];
		bool inverse = false;
		TVector3 v01, v02, v03, n;
		for (int j = 0; j < 4; ++j) 
		{
			v01.SetSubtract( m_verts[ tet.v[order[j][1]] ], m_verts[ tet.v[order[j][0]] ]);
			v02.SetSubtract( m_verts[ tet.v[order[j][2]] ], m_verts[ tet.v[order[j][0]] ]);
		
			tet.norm[j].SetCrossProd( v01, v02 );
			tet.norm[j].Normalize_Self();
		}
	}
}



void TTetraModel::calcNeighbor() 
{
	static int order[][3] = {{3, 2, 1}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};

	vector<vector<int>> vertexNeighbor( m_verts.size() );	
	for ( int i = 0; i < (int)m_tetras.size(); ++i) 
	{
		TTetra& tet = m_tetras[i];
		for (int j = 0; j < 4; ++j)  vertexNeighbor[ tet.v[j] ].push_back(i);
	}

	for (int i = 0; i < (int)m_tetras.size(); ++i) 
	{
		TTetra& tet = m_tetras[i];
		
		for (int j = 0; j < 4; ++j) 
		{
			vector<int>& n0 = vertexNeighbor[ tet.v[order[j][0]] ];//intersectionを計算
			vector<int>& n1 = vertexNeighbor[ tet.v[order[j][1]] ];
			vector<int>& n2 = vertexNeighbor[ tet.v[order[j][2]] ];
			vector<int> n01 ; set_intersection(n0 .begin(), n0 .end(), n1.begin(), n1.end(), back_inserter(n01 ));
			vector<int> n012; set_intersection(n01.begin(), n01.end(), n2.begin(), n2.end(), back_inserter(n012));

			if      ((int)n012.size() == 2) tet.n[j] = (n012.front() == i) ? n012.back() : n012.front(); 
			else if ((int)n012.size() == 1) tet.n[j] = -1;
			else {
				fprintf( stderr, "error in calcNeighbor()   %d  !!!!!!!!!!\n", i);
				return;
			}
		}
	}
}

void TTetraModel::checkStrangeTetra()
{
	for ( int i = 0; i < (int)m_tetras.size(); ++i) 
	{
		if( m_tetras[i].v[0] == m_tetras[i].v[1] || 
			m_tetras[i].v[0] == m_tetras[i].v[2] || 
			m_tetras[i].v[0] == m_tetras[i].v[3] || 
			m_tetras[i].v[1] == m_tetras[i].v[2] || 
			m_tetras[i].v[1] == m_tetras[i].v[3] || 
			m_tetras[i].v[2] == m_tetras[i].v[3] ) fprintf( stderr, "strange tetra is exist! %d\n", i);
	}


}
void TTetraModel::checkVtxOrder()
{
	TVector3 v2, v3, dir, n;
	for ( int i = 0; i < (int)m_tetras.size(); ++i) 
	{
		TTetra& tet = m_tetras[i];
		dir.SetSubtract( m_verts[ tet.v[0] ], m_verts[ tet.v[1] ] );
		v2 .SetSubtract( m_verts[ tet.v[0] ], m_verts[ tet.v[1] ] );
		v3 .SetSubtract( m_verts[ tet.v[0] ], m_verts[ tet.v[1] ] );

		n.SetCrossProd( v3, v2);
		if( n*dir < 0)
		{
			fprintf( stderr, "strange arrangement!!!!!!!!!!!!!!!!!!!!!!!!!");
			//ひっくり返してもいいかも
		}
	}	
}


void TTetraModel::calcEdgesFromTetras()
{
	static const int edg_mat[6][2] = {{0,1},{0,2},{0,3}, 
	                                  {1,2},{2,3},{3,1} } ;

	m_edges.clear();

	vector< list<int> > emanatEs( m_verts.size() );

	for ( int i = 0; i < (int)m_tetras.size(); ++i)
	{
		TTetra &tet = m_tetras[i];
		for( int i = 0; i < 6; i++ )
		{
			TTetEdge e( tet.v[ edg_mat[i][0] ], tet.v[ edg_mat[i][1] ]);
			if( e.v[0] > e.v[1] ) std::swap( e.v[0], e.v[1] ); 
			
			list<int>::iterator eIt ;
			for( eIt = emanatEs[ e.v[0] ].begin(); eIt != emanatEs[e.v[0]].end() ; eIt++ )
			{
				if( m_edges[*eIt].v[1] == e.v[1] ) break ;
			}


			if( eIt == emanatEs[ e.v[0] ].end() )
			{
				emanatEs[ e.v[0] ].push_back( (int)m_edges.size() ) ;
				emanatEs[ e.v[1] ].push_back( (int)m_edges.size() ) ;
				tet.e[i]  =                   (int)m_edges.size();
				m_edges.push_back( e ) ;
			}
			else
				tet.e[i] = *eIt;
		}
	}
}







double TTetraModel::calcVolume() const{
	return calcVolume( m_verts, m_tetras );	
}

double TTetraModel::calcVolume(const vector<TVector3> &verts, const vector<TTetra> &tetras)
{
	TVector3 v1,v2,v3;

	double res = 0;
	for( vector<TTetra>::const_iterator it = tetras.begin(); it != tetras.end(); ++it)
	{
		v1.SetSubtract( verts[ it->v[1] ], verts[ it->v[0] ] );
		v2.SetSubtract( verts[ it->v[2] ], verts[ it->v[0] ] );
		v3.SetSubtract( verts[ it->v[3] ], verts[ it->v[0] ] );
		
		double vol = t_V1crosV2_multV3( v2, v3, v1 );

		if( vol < 0 ){
			vol = 0 ;
			//vol *= -1.0; 
			//fprintf( stderr, "vol is negative\n" );
		}
		res += vol;
	}
	return res/6.0;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Compute Surface////////////////////////////////////////////////////////////////////////////////////////////
//map_tetV2polyV : tetraが変形された場合に、その変形からmeshを求めるためのmap
//map_tetV2polyV[polyVertexID] == tetraVertexID

void TTetraModel::toPolygonModel( const vector< TTetra   > &tetras  ,
		                          const vector< TVector3 > &tetVerts,

									    vector< TVector3  > &Vs,
									    vector< TTriangle > &Ps,
									    vector< int       > &map_tetV2polyV )
{
	map_tetV2polyV.clear();
	Vs.clear();
	Ps.clear();

	for (int t = 0; t< (int)tetras.size(); ++t) 
	{
		const TTetra& tet = tetras[t];

		for(int j = 0; j < 4; ++j) if( tet.n[j] == -1)
		{
			Ps.push_back( TTriangle() );

			for (int k = 0; k < 3; ++k) 
			{
				int vtxID = tet.v[ TTetraModel::FACE_VTX_ORDER[j][k] ];

				vector<int>::iterator findPos = find( map_tetV2polyV.begin(), map_tetV2polyV.end(), vtxID );
				if (findPos == map_tetV2polyV.end()) 
				{
					map_tetV2polyV.push_back( vtxID            );
					Vs            .push_back( tetVerts[vtxID]  );
					Ps.back().idx[k] = (int)map_tetV2polyV.size() - 1;
				}
				else 
					Ps.back().idx[k] = (int) distance( map_tetV2polyV.begin(), findPos );
			}
		}
	}
}

bool TTetraModel::toPolygonModel( const vector< TTetra  > &tetras  ,
								  const vector<TVector3 > &tetVerts, 
								        TTriangleMesh     &surface,
										vector<   int   > &map_tetV2polyV )
{
	vector< TVector3   > Vs;
	vector< TTriangle  > Ps;
	toPolygonModel(  tetras, tetVerts, Vs, Ps, map_tetV2polyV );
	return surface.initFromVertsPolys( Vs, Ps );
}

void TTetraModel::toPolygonModel( TTriangleMesh &resultMesh) const{
	vector<int> tmp;
	toPolygonModel( m_tetras, m_verts, resultMesh, tmp); 
}

void TTetraModel::toPolygonModel( TTriangleMesh &resultMesh, vector< int  >  &map_tetV2polyV) const{
	toPolygonModel( m_tetras, m_verts, resultMesh, map_tetV2polyV );
}

//頂点列はnewVersを利用
void TTetraModel::toPolygonModel( TTriangleMesh &resultMesh, const vector<TVector3> &newVerts ) const {
	if( m_verts.size() != newVerts.size() ) {fprintf( stderr, "cant create surface model\n") ; return;} 
	vector<int> tmp;
	toPolygonModel( m_tetras, newVerts, resultMesh,  tmp); 
}
//頂点列はnewVersを利用
void TTetraModel::toPolygonModel( TTriangleMesh &resultMesh, const vector<TVector3> &newVerts, vector<int> &map_tetV2polyV) const {
	if( m_verts.size() != newVerts.size() ) {fprintf( stderr, "cant create surface model\n") ; return;} 
	toPolygonModel( m_tetras, newVerts, resultMesh, map_tetV2polyV); 
}







inline static void subdiv_1(const TTetra &t, const vector<int> &map_e2v,  const vector<byte> &e_canDevide, vector<TTetra> &Ts)
{
	const int *v = t.v, *e = t.e;
	if     ( e_canDevide[e[0]] ){ Ts.push_back(TTetra( map_e2v[e[0]], v[0], v[3], v[2])); Ts.push_back(TTetra( map_e2v[e[0]], v[1], v[2], v[3])); }
	else if( e_canDevide[e[1]] ){ Ts.push_back(TTetra( map_e2v[e[1]], v[0], v[1], v[3])); Ts.push_back(TTetra( map_e2v[e[1]], v[1], v[2], v[3])); }
	else if( e_canDevide[e[2]] ){ Ts.push_back(TTetra( map_e2v[e[2]], v[0], v[2], v[1])); Ts.push_back(TTetra( map_e2v[e[2]], v[1], v[2], v[3])); }
	else if( e_canDevide[e[3]] ){ Ts.push_back(TTetra( map_e2v[e[3]], v[0], v[1], v[3])); Ts.push_back(TTetra( map_e2v[e[3]], v[0], v[3], v[2])); }
	else if( e_canDevide[e[4]] ){ Ts.push_back(TTetra( map_e2v[e[4]], v[0], v[1], v[3])); Ts.push_back(TTetra( map_e2v[e[4]], v[0], v[2], v[1])); }
	else if( e_canDevide[e[5]] ){ Ts.push_back(TTetra( map_e2v[e[5]], v[0], v[3], v[2])); Ts.push_back(TTetra( map_e2v[e[5]], v[0], v[2], v[1])); }
}


inline static void subdiv_2(const TTetra &t, const vector<int > &map_e2v, const vector<byte> &e_canDevide, vector<TTetra> &Ts)
{
	TTetra rt(t);
	bool oppoEdge = false;
	if     (e_canDevide[t.e[0]] && e_canDevide[t.e[1]] ){ }
	else if(e_canDevide[t.e[0]] && e_canDevide[t.e[2]] ){ rt.rotR(0); }
	else if(e_canDevide[t.e[0]] && e_canDevide[t.e[3]] ){ rt.rotR(3); }
	else if(e_canDevide[t.e[0]] && e_canDevide[t.e[4]] ){ oppoEdge = true;}
	else if(e_canDevide[t.e[0]] && e_canDevide[t.e[5]] ){ rt.rotL(2); rt.rotR(0);}
	else if(e_canDevide[t.e[1]] && e_canDevide[t.e[2]] ){ rt.rotL(0);}
	else if(e_canDevide[t.e[1]] && e_canDevide[t.e[3]] ){ rt.rotL(3);}
	else if(e_canDevide[t.e[1]] && e_canDevide[t.e[4]] ){ rt.rotR(1); rt.rotL(0);}
	else if(e_canDevide[t.e[1]] && e_canDevide[t.e[5]] ){ rt.rotL(0);             oppoEdge = true; }
	else if(e_canDevide[t.e[2]] && e_canDevide[t.e[3]] ){ rt.rotR(0);             oppoEdge = true; }
	else if(e_canDevide[t.e[2]] && e_canDevide[t.e[4]] ){ rt.rotL(1); rt.rotL(0);}
	else if(e_canDevide[t.e[2]] && e_canDevide[t.e[5]] ){ rt.rotR(2); rt.rotR(0);}
	else if(e_canDevide[t.e[3]] && e_canDevide[t.e[4]] ){ rt.rotL(3); rt.rotL(0);}
	else if(e_canDevide[t.e[3]] && e_canDevide[t.e[5]] ){ rt.rotR(3); rt.rotR(0);}
	else if(e_canDevide[t.e[4]] && e_canDevide[t.e[5]] ){ rt.rotL(1); rt.rotR(0);}

	int *v = rt.v, *e = rt.e; 
	for( int i=0; i<6; ++i) e[i] = map_e2v[ e[i] ];
	
	if( oppoEdge ){
		Ts.push_back( TTetra( v[0], v[2], e[4], e[0]));
		Ts.push_back( TTetra( v[0], v[3], e[0], e[4]));
		Ts.push_back( TTetra( v[1], v[2], e[0], e[4]));
		Ts.push_back( TTetra( v[1], v[3], e[4], e[0]));
		return;
	}
	else{
		Ts.push_back( TTetra( v[3], v[0], e[1], e[0] ) );
		if( e[0] < e[1]){ Ts.push_back( TTetra( v[3], e[0], e[1], v[2] ) );
						  Ts.push_back( TTetra( v[3], e[0], v[2], v[1] ) );}
		else{             Ts.push_back( TTetra( v[3], e[0], e[1], v[1] ) );
						  Ts.push_back( TTetra( v[3], e[1], v[2], v[1] ) );}
		}
}

inline static void subdiv_3(const TTetra &t, const vector<int > &map_e2v, const vector<byte> &e_canDevide, vector<TTetra> &Ts)
{
	TTetra rt(t);
	int caseid = -1;
	byte be[6]; for( int i=0; i<6; ++i) be[i] = e_canDevide[t.e[i]];

	if     (be[0] && be[1] && be[2] ){ caseid = 0;}
	else if(be[0] && be[1] && be[3] ){ caseid = 1;}
	else if(be[0] && be[1] && be[4] ){ caseid = 2;}
	else if(be[0] && be[1] && be[5] ){ caseid = 3;}
	else if(be[0] && be[2] && be[3] ){ caseid = 2; rt.rotR(0); }
	else if(be[0] && be[2] && be[4] ){ caseid = 3; rt.rotR(0); }
	else if(be[0] && be[2] && be[5] ){ caseid = 1; rt.rotR(0); }
	else if(be[0] && be[3] && be[4] ){ caseid = 3; rt.rotR(3); }
	else if(be[0] && be[3] && be[5] ){ caseid = 0; rt.rotR(3); }
	else if(be[0] && be[4] && be[5] ){ caseid = 2; rt.rotL(2); rt.rotR(0); }

	else if(be[1] && be[2] && be[3] ){ caseid = 3; rt.rotL(0); }
	else if(be[1] && be[2] && be[4] ){ caseid = 1; rt.rotL(0); }
	else if(be[1] && be[2] && be[5] ){ caseid = 2; rt.rotL(0); }
	else if(be[1] && be[3] && be[4] ){ caseid = 0; rt.rotL(3);}
	else if(be[1] && be[3] && be[5] ){ caseid = 2; rt.rotL(3);}
	else if(be[1] && be[4] && be[5] ){ caseid = 3; rt.rotR(1); rt.rotL(0);}

	else if(be[2] && be[3] && be[4] ){ caseid = 2; rt.rotL(1); rt.rotL(0);}
	else if(be[2] && be[3] && be[5] ){ caseid = 3; rt.rotR(2); rt.rotR(0);}
	else if(be[2] && be[4] && be[5] ){ caseid = 0; rt.rotL(1); }
	else if(be[3] && be[4] && be[5] ){ caseid = 1; rt.rotL(2); }

	int *e = rt.e, *v = rt.v;
	for( int i=0; i<6; ++i) e[i] = map_e2v[ e[i] ];
	if( caseid == 0)
	{
		Ts.push_back(TTetra( v[0], e[0], e[1], e[2]) );
		if     (e[0]<e[1]&&e[1]<e[2]){Ts.push_back(TTetra(v[1],v[3],v[2],e[0])); Ts.push_back(TTetra(v[2],v[3],e[1],e[0])); Ts.push_back(TTetra(v[3],e[2],e[1],e[0]));}
		else if(e[0]<e[2]&&e[2]<e[1]){Ts.push_back(TTetra(v[1],v[3],v[2],e[0])); Ts.push_back(TTetra(v[2],e[2],e[1],e[0])); Ts.push_back(TTetra(v[3],e[2],v[2],e[0]));}
		else if(e[1]<e[0]&&e[0]<e[2]){Ts.push_back(TTetra(v[1],v[3],e[1],e[0])); Ts.push_back(TTetra(v[2],v[3],e[1],v[1])); Ts.push_back(TTetra(v[3],e[2],e[1],e[0]));}
		else if(e[1]<e[2]&&e[2]<e[0]){Ts.push_back(TTetra(v[1],e[2],e[1],e[0])); Ts.push_back(TTetra(v[2],v[3],e[1],v[1])); Ts.push_back(TTetra(v[3],e[2],e[1],v[1]));}
		else if(e[2]<e[0]&&e[0]<e[1]){Ts.push_back(TTetra(v[1],e[2],v[2],e[0])); Ts.push_back(TTetra(v[2],e[2],e[1],e[0])); Ts.push_back(TTetra(v[3],e[2],v[2],v[1]));}
		else if(e[2]<e[1]&&e[1]<e[0]){Ts.push_back(TTetra(v[1],e[2],e[1],e[0])); Ts.push_back(TTetra(v[2],e[2],e[1],v[1])); Ts.push_back(TTetra(v[3],e[2],v[2],v[1]));}
	}
	else if( caseid == 1){
		Ts.push_back(TTetra( v[3], v[0], e[1], e[0]) );
		Ts.push_back(TTetra( v[3], v[2], e[3], e[1]) );
		Ts.push_back(TTetra( v[3], v[1], e[0], e[3]) );
		Ts.push_back(TTetra( v[3], e[0], e[1], e[3]) );
	}
	else if( caseid == 2)
	{
		Ts.push_back(TTetra( e[0], v[3], v[1], e[4]) );
		if( e[0] < e[1] ){ Ts.push_back(TTetra( e[4], e[0], e[1], v[2]) ); Ts.push_back(TTetra( e[4], e[0], v[2], v[1]) ); }
		else {             Ts.push_back(TTetra( e[4], e[0], e[1], v[1]) ); Ts.push_back(TTetra( e[4], e[1], v[2], v[1]) ); }
		if( e[1] < e[4] ){ Ts.push_back(TTetra( e[0], e[4], e[1], v[3]) ); Ts.push_back(TTetra( e[0], e[1], v[0], v[3]) ); }
		else {             Ts.push_back(TTetra( e[0], e[4], e[1], v[0]) ); Ts.push_back(TTetra( e[0], e[4], v[0], v[3]) ); }
	}
	else if( caseid == 3)
	{
		Ts.push_back(TTetra( e[1], v[2], v[3], e[5]) );
		if( e[0] < e[5] ){ Ts.push_back(TTetra( e[1], e[0], v[3], v[0]) ); Ts.push_back(TTetra( e[1], e[0], e[5], v[3]) ); }
		else {             Ts.push_back(TTetra( e[1], e[5], v[3], v[0]) ); Ts.push_back(TTetra( e[1], e[0], e[5], v[0]) ); }
		if( e[0] < e[1] ){ Ts.push_back(TTetra( e[5], e[0], e[1], v[2]) ); Ts.push_back(TTetra( e[5], e[0], v[2], v[1]) ); }
	    else {             Ts.push_back(TTetra( e[5], e[0], e[1], v[1]) ); Ts.push_back(TTetra( e[5], e[1], v[2], v[1]) ); }

	}
}




inline static void subdiv_4(const TTetra &t,const vector<int > &map_e2v, const vector<byte> &e_canDevide, vector<TTetra> &Ts)
{
	TTetra rt(t);
	bool oppoEdge = false;
	if     (!e_canDevide[t.e[0]] && !e_canDevide[t.e[1]] ){ }
	else if(!e_canDevide[t.e[0]] && !e_canDevide[t.e[2]] ){ rt.rotR(0); }
	else if(!e_canDevide[t.e[0]] && !e_canDevide[t.e[3]] ){ rt.rotR(3); }
	else if(!e_canDevide[t.e[0]] && !e_canDevide[t.e[4]] ){ oppoEdge = true;}
	else if(!e_canDevide[t.e[0]] && !e_canDevide[t.e[5]] ){ rt.rotL(2); rt.rotR(0);}
	else if(!e_canDevide[t.e[1]] && !e_canDevide[t.e[2]] ){ rt.rotL(0);}
	else if(!e_canDevide[t.e[1]] && !e_canDevide[t.e[3]] ){ rt.rotL(3);}
	else if(!e_canDevide[t.e[1]] && !e_canDevide[t.e[4]] ){ rt.rotR(1); rt.rotL(0);}
	else if(!e_canDevide[t.e[1]] && !e_canDevide[t.e[5]] ){ rt.rotL(0);             oppoEdge = true; }
	else if(!e_canDevide[t.e[2]] && !e_canDevide[t.e[3]] ){ rt.rotR(0);             oppoEdge = true; }
	else if(!e_canDevide[t.e[2]] && !e_canDevide[t.e[4]] ){ rt.rotL(1); rt.rotL(0);}
	else if(!e_canDevide[t.e[2]] && !e_canDevide[t.e[5]] ){ rt.rotR(2); rt.rotR(0);}
	else if(!e_canDevide[t.e[3]] && !e_canDevide[t.e[4]] ){ rt.rotL(3); rt.rotL(0);}
	else if(!e_canDevide[t.e[3]] && !e_canDevide[t.e[5]] ){ rt.rotR(3); rt.rotR(0);}
	else if(!e_canDevide[t.e[4]] && !e_canDevide[t.e[5]] ){ rt.rotL(1); rt.rotR(0);}

	int *e = rt.e, *v = rt.v;
	for( int i=0; i<6; ++i) e[i] = map_e2v[ e[i] ];
	if( oppoEdge )
	{
		TTetra T1, T2, T3, T4, T5, T6;
		//e1e2e3e5の全ての順序を考慮
		if( e[1]<e[2] && e[1]<e[3] && e[1]<e[5] )
		{
			T1.Set(e[1],e[3],e[5],v[1]); 			
			if( e[2] < e[5] ){ T2.Set(e[1],e[5],e[2],v[1]); T3.Set(e[1],e[2],v[0],v[1]); }
			else             { T2.Set(e[1],e[5],e[2],v[0]); T3.Set(e[1],e[5],v[0],v[1]); }
			T4.Set(e[1],e[2],e[5],v[3]); 
			if( e[3] < e[5] ){ T5.Set(e[1],e[5],e[3],v[3]); T6.Set(e[1],e[3],v[2],v[3]); }
			else             { T5.Set(e[1],e[5],e[3],v[2]); T6.Set(e[1],e[5],v[2],v[3]); }
		}
		else if( e[2]<e[3] && e[2]<e[5] )
		{
			T1.Set(e[2],e[3],e[5],v[1]); 			
			if( e[1] < e[3] ){ T2.Set(e[2],e[1],e[3],v[1]); T3.Set(e[2],e[1],v[1],v[0]); }
			else             { T2.Set(e[2],e[1],e[3],v[0]); T3.Set(e[2],e[3],v[1],v[0]); }
			T4.Set(e[2],e[3],e[1],v[2]); 
			if( e[3] < e[5] ){ T5.Set(e[2],e[5],e[3],v[3]); T6.Set(e[2],e[3],v[2],v[3]); }
			else             { T5.Set(e[2],e[5],e[3],v[2]); T6.Set(e[2],e[5],v[2],v[3]); }
		}
		else if( e[3]<e[5] )
		{
			T1.Set(e[3],e[2],e[1],v[0]); 			
			if( e[2] < e[5] ){ T2.Set(e[3],e[5],e[2],v[1]); T3.Set(e[3],e[2],v[0],v[1]); }
			else             { T2.Set(e[3],e[5],e[2],v[0]); T3.Set(e[3],e[5],v[0],v[1]); }
			T4.Set(e[3],e[2],e[5],v[3]); 
			if( e[1] < e[2] ){ T5.Set(e[3],e[1],e[2],v[3]); T6.Set(e[3],e[1],v[3],v[2]); }
			else             { T5.Set(e[3],e[1],e[2],v[2]); T6.Set(e[3],e[2],v[3],v[2]); }
		}
		else
		{
			T1.Set(e[5],e[2],e[1],v[0]); 			
			if( e[1] < e[3] ){ T2.Set(e[5],e[1],e[3],v[1]); T3.Set(e[5],e[1],v[1],v[0]); }
			else             { T2.Set(e[5],e[1],e[3],v[0]); T3.Set(e[5],e[3],v[1],v[0]); }
			T4.Set(e[5],e[3],e[1],v[2]); 
			if( e[1] < e[2] ){ T5.Set(e[5],e[1],e[2],v[3]); T6.Set(e[5],e[1],v[3],v[2]); }
			else             { T5.Set(e[5],e[1],e[2],v[2]); T6.Set(e[5],e[2],v[3],v[2]); }
		}
		Ts.push_back( T1 );
		Ts.push_back( T2 );
		Ts.push_back( T3 );
		Ts.push_back( T4 );
		Ts.push_back( T5 );
		Ts.push_back( T6 );


		return;
	}
	else{
		Ts.push_back(TTetra( v[3], e[2], e[4], e[5]) );
		Ts.push_back(TTetra( e[3], e[2], e[5], e[4]) );
		if( e[2] < e[4]){ Ts.push_back(TTetra( e[3], v[2], v[0], e[2])); Ts.push_back(TTetra( e[3], v[2], e[2], e[4] )); }
		else{             Ts.push_back(TTetra( e[3], v[2], v[0], e[4])); Ts.push_back(TTetra( e[3], v[0], e[2], e[4] )); }
		if( e[2] < e[5]){ Ts.push_back(TTetra( e[3], v[0], v[1], e[2])); Ts.push_back(TTetra( e[3], v[1], e[5], e[2] )); }
		else{             Ts.push_back(TTetra( e[3], v[0], v[1], e[5])); Ts.push_back(TTetra( e[3], v[0], e[5], e[2] )); }
	}
}



inline static void subdiv_5(const TTetra &t, const vector<int > &_map_e2v, const vector<byte> &e_canDevide, vector<TTetra> &Ts)
{
	TTetra rt(t);
	if(      !e_canDevide[ t.e[0] ] ){           }
	else if( !e_canDevide[ t.e[1] ] ){rt.rotL(0);}
	else if( !e_canDevide[ t.e[2] ] ){rt.rotR(0);}
	else if( !e_canDevide[ t.e[3] ] ){rt.rotR(3);}
	else if( !e_canDevide[ t.e[4] ] ){rt.rotR(2); rt.rotL(0);}
	else if( !e_canDevide[ t.e[5] ] ){rt.rotL(1);}

	int *v = rt.v;
	int  e[6]; for( int i=0; i<6; ++i) e[i] = _map_e2v[ rt.e[i] ];
	//const int *e = rt.e;
	Ts.push_back( TTetra( v[2], e[1], e[3], e[4] ) );
	Ts.push_back( TTetra( v[3], e[2], e[4], e[5] ) );

	if( e[1]<e[2] && e[1]<e[3] && e[1]<e[5] )
	{
		Ts.push_back( TTetra( e[4], e[1], e[5], e[2]) );
		Ts.push_back( TTetra( e[4], e[1], e[3], e[5]) );

		Ts.push_back( TTetra( e[1], e[3], e[5], v[1]) );
		if( e[2] < e[5] ){Ts.push_back( TTetra( e[1], e[5], e[2], v[1]) );  Ts.push_back( TTetra( e[1], e[2], v[0], v[1]) );}
		else             {Ts.push_back( TTetra( e[1], e[5], e[2], v[0]) );  Ts.push_back( TTetra( e[1], e[5], v[0], v[1]) );}
	}
	else if( e[2]<e[3] && e[2]<e[5] )
	{
		Ts.push_back( TTetra( e[4], e[2], e[1], e[3]) );
		Ts.push_back( TTetra( e[4], e[2], e[3], e[5]) );

		Ts.push_back( TTetra( e[2], e[3], e[5], v[1]) );
		if( e[1] < e[3] ){Ts.push_back( TTetra( e[2], e[1], e[3], v[1]) );  Ts.push_back( TTetra( e[2], e[1], v[1], v[0]) );}
		else             {Ts.push_back( TTetra( e[2], e[1], e[3], v[0]) );  Ts.push_back( TTetra( e[2], e[3], v[1], v[0]) );}
	}
	else if( e[3]<e[5] )
	{
		Ts.push_back( TTetra( e[4], e[2], e[1], e[3]) );
		Ts.push_back( TTetra( e[4], e[2], e[3], e[5]) );

		Ts.push_back( TTetra( e[3], e[2], e[1], v[0]) );
		if( e[2] < e[5] ){Ts.push_back( TTetra( e[3], e[5], e[2], v[1]) );  Ts.push_back( TTetra( e[3], e[2], v[0], v[1]) );}
		else             {Ts.push_back( TTetra( e[3], e[5], e[2], v[0]) );  Ts.push_back( TTetra( e[3], e[5], v[0], v[1]) );}
	}
	else
	{
		Ts.push_back( TTetra( e[4], e[1], e[5], e[2]) );
		Ts.push_back( TTetra( e[4], e[1], e[3], e[5]) );

		Ts.push_back( TTetra( e[5], e[2], e[1], v[0]) );
		if( e[1] < e[3] ){Ts.push_back( TTetra( e[5], e[1], e[3], v[1]) );  Ts.push_back( TTetra( e[5], e[1], v[1], v[0]) );}
		else             {Ts.push_back( TTetra( e[5], e[1], e[3], v[0]) );  Ts.push_back( TTetra( e[5], e[3], v[1], v[0]) );}
	}
}

inline static void subdiv_6(const TTetra       &t       , 
							const vector<int > &_map_e2v, 
							const vector<byte> &e_canDevide,
							vector<TTetra    > &Ts,
							const vector<TVector3> &Vs )
{
	int e[6]; for( int i=0; i<6; ++i) e[i] = _map_e2v[ t.e[i] ];

	Ts.push_back( TTetra( t.v[0], e[0], e[1], e[2] ));
	Ts.push_back( TTetra( t.v[1], e[0], e[5], e[3] ));
	Ts.push_back( TTetra( t.v[2], e[1], e[3], e[4] ));
	Ts.push_back( TTetra( t.v[3], e[2], e[4], e[5] ));
	
	double e0e4 = t_distance_sq( Vs[e[0]], Vs[e[4]]);
	double e1e5 = t_distance_sq( Vs[e[1]], Vs[e[5]]);
	double e2e3 = t_distance_sq( Vs[e[2]], Vs[e[3]]);
	
	if( e0e4 <= e1e5 && e0e4 <= e2e3 ) 
	{
		Ts.push_back( TTetra( e[0], e[1], e[2], e[4] ));
		Ts.push_back( TTetra( e[0], e[2], e[5], e[4] ));
		Ts.push_back( TTetra( e[0], e[5], e[3], e[4] ));
		Ts.push_back( TTetra( e[0], e[3], e[1], e[4] ));
	}
	else if( e1e5 <= e2e3 ) 
	{
		Ts.push_back( TTetra( e[1], e[2], e[0], e[5] ));
		Ts.push_back( TTetra( e[1], e[0], e[3], e[5] ));
		Ts.push_back( TTetra( e[1], e[3], e[4], e[5] ));
		Ts.push_back( TTetra( e[1], e[4], e[2], e[5] ));
	}
	else
	{
		Ts.push_back( TTetra( e[2], e[0], e[1], e[3] ));
		Ts.push_back( TTetra( e[2], e[1], e[4], e[3] ));
		Ts.push_back( TTetra( e[2], e[4], e[5], e[3] ));
		Ts.push_back( TTetra( e[2], e[5], e[0], e[3] ));
	}

}



void TTetraModel::subdivision(vector<byte> &trgtTet, bool doEdgeThresholding, double edgeLength)
{
	vector<TTetra  >  oldTetras = m_tetras;
	vector<TTetEdge>  oldEdges  = m_edges ;
	vector<TVector3>  oldVerts  = m_verts ;

	//分割できるedgeを計算
	vector<byte> e_canDevide( oldEdges.size(), true);

	if( doEdgeThresholding ){
		for( int i=0; i<(int)oldEdges.size(); ++i)
			if( t_distance( oldVerts[ oldEdges[i].v[0]] , oldVerts[ oldEdges[i].v[1]] ) < edgeLength )  e_canDevide[i] = false;
			else                                                                                        e_canDevide[i] = true ;
	}

	for( int i=0; i< (int)oldTetras.size(); ++i) if( !trgtTet[i] ){
		e_canDevide[ oldTetras[i].e[0] ] = e_canDevide[ oldTetras[i].e[1] ] = e_canDevide[ oldTetras[i].e[2] ] = false;
		e_canDevide[ oldTetras[i].e[3] ] = e_canDevide[ oldTetras[i].e[4] ] = e_canDevide[ oldTetras[i].e[5] ] = false;
	}

	int DivideEdgeNum = 0;
	for( int i=0; i<(int)oldEdges.size(); ++i) if( e_canDevide[i] ) DivideEdgeNum++;


	//edgeを分割
	vector<int > map_e2v( oldEdges.size(), -1);//edge 分割後に生成されるvertexId
	for( int i=0; i<(int)oldEdges.size(); ++i)if( e_canDevide[i] )
	{
		map_e2v[i] = (int)m_verts.size();
		m_verts.push_back( 0.5 * ( m_verts[oldEdges[i].v[0]] + m_verts[oldEdges[i].v[1]] ) );
	}

	//tetraを挿入
	m_tetras.clear();
	for( int i=0; i< (int)oldTetras.size(); ++i)
	{
		TTetra &t = oldTetras[i];
		int num = 0;
		for( int kk=0; kk<6; ++kk) if( e_canDevide[ t.e[kk] ] ) num++;
		if     ( num == 0 ) m_tetras.push_back( TTetra( t.v[0], t.v[1], t.v[2], t.v[3] ) );
		else if( num == 1 ) subdiv_1(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 2 ) subdiv_2(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 3 ) subdiv_3(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 4 ) subdiv_4(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 5 ) subdiv_5(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 6 ) subdiv_6(t, map_e2v, e_canDevide, m_tetras, m_verts );
	}

	fprintf( stderr, "devided Edges = %d    new v %d  new t %d\n", DivideEdgeNum, m_verts.size(), m_tetras.size());
}


/*
edge_toDivide[i] = trueのedgeをsubdivideする
分割点位置は　edge_divPos[i] 
 
新しい頂点  は m_verts の後ろに積まれる
新しいテトラは m_tetrasの後ろに積まれる
*/
void TTetraModel::subdivision_divideEdge( const vector<byte    > &edge_toDivide, //edge[i]を割るかどうかのフラグ 
	                                      const vector<TVector3> &edge_divPos   )//edge[i]をdivideしたときの頂点一
{
	calcNeighbor();
	calcEdgesFromTetras();
	vector<TTetra  >  oldTetras = m_tetras;
	vector<TTetEdge>  oldEdges  = m_edges ;
	vector<TVector3>  oldVerts  = m_verts ;

	//分割できるedgeを計算

	int DivideEdgeNum = 0;
	for( int i=0; i<(int)oldEdges.size(); ++i) if( edge_toDivide[i] ) DivideEdgeNum++;


	//edgeを分割
	vector<int > map_e2v( oldEdges.size(), -1);//edge 分割後に生成されるvertexId
	for( int i=0; i<(int)oldEdges.size(); ++i)if( edge_toDivide[i] )
	{
		map_e2v[i] = (int)m_verts.size();
		m_verts.push_back( edge_divPos[i] );
	}

	//tetraを挿入
	m_tetras.clear();
	for( int i=0; i< (int)oldTetras.size(); ++i)
	{
		TTetra &t = oldTetras[i];
		int num = 0;
		for( int kk=0; kk<6; ++kk) if( edge_toDivide[ t.e[kk] ] ) num++;
		if     ( num == 0 ) m_tetras.push_back( TTetra( t.v[0], t.v[1], t.v[2], t.v[3] ) );
		else if( num == 1 ) subdiv_1(t, map_e2v, edge_toDivide, m_tetras );
		else if( num == 2 ) subdiv_2(t, map_e2v, edge_toDivide, m_tetras );
		else if( num == 3 ) subdiv_3(t, map_e2v, edge_toDivide, m_tetras );
		else if( num == 4 ) subdiv_4(t, map_e2v, edge_toDivide, m_tetras );
		else if( num == 5 ) subdiv_5(t, map_e2v, edge_toDivide, m_tetras );
		else if( num == 6 ) subdiv_6(t, map_e2v, edge_toDivide, m_tetras, m_verts );
	}

	calcNeighbor();
	calcEdgesFromTetras();
	fprintf( stderr, "devided Edges = %d    new v %d  new t %d\n", DivideEdgeNum, m_verts.size(), m_tetras.size());
}




void TTetraModel::subdivision_notDivBoundary(vector<byte> &trgtTet)
{
	vector<TTetra    > oldTetras = m_tetras;
	vector<TTetEdge  > oldEdges  = m_edges ;
	vector<TVector3> oldVerts  = m_verts ;

	//分割できるedgeを計算
	vector<byte> e_canDevide( oldEdges.size(), true);//edge 分割後に生成されるvertexId
	for( int i=0; i< (int)oldTetras.size(); ++i) {
		TTetra &t = m_tetras[i];
		if( !trgtTet[i] ){
			e_canDevide[ t.e[0] ] = e_canDevide[ t.e[1] ] = e_canDevide[ t.e[2] ] = false;
			e_canDevide[ t.e[3] ] = e_canDevide[ t.e[4] ] = e_canDevide[ t.e[5] ] = false;
		}else{
			if( t.n[0] ==-1 ){ e_canDevide[ t.e[3] ] = e_canDevide[ t.e[4] ] = e_canDevide[ t.e[5] ] = false; }
			if( t.n[1] ==-1 ){ e_canDevide[ t.e[1] ] = e_canDevide[ t.e[2] ] = e_canDevide[ t.e[4] ] = false; }
			if( t.n[2] ==-1 ){ e_canDevide[ t.e[0] ] = e_canDevide[ t.e[2] ] = e_canDevide[ t.e[5] ] = false; }
			if( t.n[3] ==-1 ){ e_canDevide[ t.e[0] ] = e_canDevide[ t.e[1] ] = e_canDevide[ t.e[3] ] = false; }
		}
	}

	//edgeを分割
	vector<int > map_e2v( oldEdges.size(), -1);//edge 分割後に生成されるvertexId
	for( int i=0; i<(int)oldEdges.size(); ++i)if( e_canDevide[i] )
	{
		map_e2v[i] = (int)m_verts.size();
		m_verts.push_back( 0.5 * ( m_verts[oldEdges[i].v[0]] + m_verts[oldEdges[i].v[1]] ) );
	}

	m_tetras.clear();
	for( int i=0; i< (int)oldTetras.size(); ++i)
	{
		TTetra &t = oldTetras[i];
		int num = 0;
		for( int kk=0; kk<6; ++kk) if( e_canDevide[ t.e[kk] ] ) num++;
		if     ( num == 0 ) m_tetras.push_back( TTetra( t.v[0], t.v[1], t.v[2], t.v[3] ) );
		else if( num == 1 ) subdiv_1(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 2 ) subdiv_2(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 3 ) subdiv_3(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 4 ) subdiv_4(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 5 ) subdiv_5(t, map_e2v, e_canDevide, m_tetras );
		else if( num == 6 ) subdiv_6(t, map_e2v, e_canDevide, m_tetras, m_verts );
	}
}





static int numOfIntersection( set<int> &s1, set<int> &s2)
{
	vector<int> is; 
	set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(), back_inserter(is));
	return (int)is.size();

}


void TTetraModel::sortVerts_By1RingRegionGrow()
{
	{
		double ave_1ring  = 0, ave_intersec = 0;
		vector<set<int>> ring1(m_verts.size() );
		for( int i=0; i<(int)m_tetras.size(); ++i) 
			for( int j=0; j<4; ++j) for( int k=0; k<4; ++k) ring1[ m_tetras[i].v[j] ].insert( m_tetras[i].v[k] );
		for( int i=0; i<(int)m_verts.size(); ++i)
		{
			ave_1ring    += ring1[i].size();
			ave_intersec += (i!=0) ? numOfIntersection( ring1[i], ring1[i-1] ) : 0;
		}
		fprintf( stderr, "1 ring頂点数 　　　 = %f\n", ave_1ring   / m_verts.size());
		fprintf( stderr, "前の1ringとの重複数 = %f\n", ave_intersec/ m_verts.size());
	}

	vector<set<int>> ring1(m_verts.size() );
	for( int i=0; i<(int)m_tetras.size(); ++i) 
		for( int j=0; j<4; ++j) for( int k=0; k<4; ++k) ring1[ m_tetras[i].v[j] ].insert( m_tetras[i].v[k] );

	//region grow //適当にpivotを決めてregion growingし，その順に入れていく
	const int vSize = (int) m_verts.size();
	vector< int > sortedVid     ( vSize, -1);
	vector< int > map_oldV2sortV( vSize, -1);

	int insertAt = 0;
	while( true )
	{
		//pivot
		int pivVertex = -1;
		for( int i=0; i<vSize; ++i) if( map_oldV2sortV[i] == -1 ) { pivVertex = i; break; }
		if( pivVertex == -1 ) break;

		//growth
		set<int> growthFront; 
		growthFront.insert( pivVertex );
		while( !growthFront.empty() )
		{	
			//popする物を選択・・前の領域とintersectionの一番大きいものを選択
			set<int>::iterator it;
			if     ( insertAt           == 0 ) it = growthFront.begin();
			else if( growthFront.size() == 1 ) it = growthFront.begin();
			else{
				int maxNum = -1;
				for( set<int>::iterator tmpIt = growthFront.begin(); tmpIt != growthFront.end(); ++tmpIt ){
					int n = numOfIntersection( ring1[*tmpIt], ring1[ sortedVid[insertAt-1]] );
					if( n > maxNum ){maxNum = n; it = tmpIt;}
				}
			}

			//*itを挿入
			int vid = *it;
			map_oldV2sortV[vid     ] = insertAt;
			sortedVid     [insertAt] = vid;
			growthFront.erase( it );
			insertAt++;

			for( set<int>::iterator tmpIt = ring1[vid].begin(); tmpIt != ring1[vid].end(); ++tmpIt) if( map_oldV2sortV[*tmpIt] == -1 ) growthFront.insert( *tmpIt ); 
		}
	}

	vector< TVector3 > newVs( m_verts.size() );
	for( int i=0; i<vSize; ++i) newVs[i] = m_verts[sortedVid[i]];
	swap( newVs, m_verts );

	const int tSize = (int) m_tetras.size();
	for( int i=0; i<tSize; ++i){
		m_tetras[i].v[0] = map_oldV2sortV[ m_tetras[i].v[0] ];
		m_tetras[i].v[1] = map_oldV2sortV[ m_tetras[i].v[1] ];
		m_tetras[i].v[2] = map_oldV2sortV[ m_tetras[i].v[2] ];
		m_tetras[i].v[3] = map_oldV2sortV[ m_tetras[i].v[3] ];
	}
	checkVtxOrder();
	calcNeighbor();
	calcNormal();
	
	{
		double ave_1ring  = 0, ave_intersec = 0;
		vector<set<int>> ring1(m_verts.size() );
		for( int i=0; i<(int)m_tetras.size(); ++i) 
			for( int j=0; j<4; ++j) for( int k=0; k<4; ++k) ring1[ m_tetras[i].v[j] ].insert( m_tetras[i].v[k] );
		for( int i=0; i<(int)m_verts.size(); ++i)
		{
			ave_1ring    += ring1[i].size();
			ave_intersec += (i!=0) ? numOfIntersection( ring1[i], ring1[i-1] ) : 0;
		}
		fprintf( stderr, "1 ring頂点数 　　　 = %f\n", ave_1ring   / m_verts.size());
		fprintf( stderr, "前の1ringとの重複数 = %f\n", ave_intersec/ m_verts.size());
	}
}







struct comp_sortVertsPolys_byXvalue{
	bool operator()(const pair<double,int>& a, const pair<double,int>&b) const { return b.first > a.first; }
};

//memory効率を上げるためにvertexをソート(x座標の値を利用)
void TTetraModel::sortVerts_ByXvalue()
{
	{
		vector<set<int>> ring1(m_verts.size() );
		for( int i=0; i<(int)m_tetras.size(); ++i){
			int *v = m_tetras[i].v;
			for( int j=0; j<4; ++j) for( int k=0; k<4; ++k) ring1[ v[j] ].insert( v[k] );
		}
		double ave_1ring  = 0, ave_intersec = 0;
		for( int i=0; i<(int)m_verts.size(); ++i){
			ave_1ring += ring1[i].size();
			if( i!=0){
				vector<int> is; 
				set_intersection( ring1[i].begin(), ring1[i].end(), ring1[i-1].begin(), ring1[i-1].end(), back_inserter(is));
				ave_intersec += is.size();
			}
		}
		fprintf( stderr, "1 ring頂点数 　　　 = %f\n", ave_1ring   / m_verts.size());
		fprintf( stderr, "前の1ringとの重複数 = %f\n", ave_intersec/ m_verts.size());
	}

	//sort//
	const int vSize = (int) m_verts.size();
	vector< pair<double,int> > sortedVid( vSize );
	for( int i=0; i<vSize; ++i) {
		sortedVid[i].first  = m_verts[i].data[0];
		sortedVid[i].second = i;
	}
	sort( sortedVid.begin(), sortedVid.end(), comp_sortVertsPolys_byXvalue()); 

	//map //
	vector< int > map_oldV2sortV( vSize, -1);
	for( int i=0; i<vSize; ++i) map_oldV2sortV[ sortedVid[i].second ] = i;

	//swap
	vector< TVector3 > newVs( m_verts.size() );
	for( int i=0; i<vSize; ++i) newVs[i] = m_verts[ sortedVid[i].second ];
	swap( newVs, m_verts );

	const int tSize = (int) m_tetras.size();
	for( int i=0; i<tSize; ++i){
		m_tetras[i].v[0] = map_oldV2sortV[ m_tetras[i].v[0] ];
		m_tetras[i].v[1] = map_oldV2sortV[ m_tetras[i].v[1] ];
		m_tetras[i].v[2] = map_oldV2sortV[ m_tetras[i].v[2] ];
		m_tetras[i].v[3] = map_oldV2sortV[ m_tetras[i].v[3] ];
	}

	checkVtxOrder();
	calcNeighbor();
	calcNormal();

	{
		vector<set<int>> ring1(m_verts.size() );
		for( int i=0; i<(int)m_tetras.size(); ++i){
			int *v = m_tetras[i].v;
			for( int j=0; j<4; ++j) for( int k=0; k<4; ++k) ring1[ v[j] ].insert( v[k] );
		}
		double ave_1ring  = 0, ave_intersec = 0;
		for( int i=0; i<(int)m_verts.size(); ++i){
			ave_1ring += ring1[i].size();
			if( i!=0){
				vector<int> is; 
				set_intersection( ring1[i].begin(), ring1[i].end(), ring1[i-1].begin(), ring1[i-1].end(), back_inserter(is));
				ave_intersec += is.size();
			}
		}
		fprintf( stderr, "1 ring頂点数 　　　 = %f\n", ave_1ring   / m_verts.size());
		fprintf( stderr, "前の1ringとの重複数 = %f\n", ave_intersec/ m_verts.size());
	}

}













//Project fllagelum用 initialization//
/*
by kenshi
v1-v2-v3は外から見て左回り
     v0--v2
    / \ .|
   /  .\ |
  / .   \|
v1.______v3

    p3 ----------------p5
     | \             /
	 |   \          / |
     |     \       /  |    tet1 = p0-p1 p2 p4 
	 |       \ p4 /   |    tet2 = p0-p2 p5 p4
	 |	       \ /    |    tet3 = p0-p3 p4 p5
	 |          |     |      
	 |          |     |
   p0 ----------|----- p2
       \        |    /
	     \      |   /
		   \    |  /
		     \  | /
			   \ /
                p1


	cylinderの断面は vNumCが6の時
     
	y+4
	 *-----* y+3
    /       \
y+5/   y     \
  *     *     * y+2 (r = 1) 
   \         /
    \       /
	 *-----* y+1 (r = 0)
   y+6

   一断面に vNum + 1個の頂点がある.
   k番目の時は, y = k * (vNum + 1)
 */
void TTetraModel::initializeAsCylinder(int vNumC, int vNumL, double rad, double length   )
{
	m_verts .clear();
	m_tetras.clear();

	double stepZ = length    / ( (double)vNumL - 1.0 );
	double stepR = 2 * M_PI  / ( (double)vNumC       );

	int p0 = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0;

	for( int k = 0; k < vNumL; ++k)
	{
		double zPos = stepZ * k;

		m_verts.push_back( TVector3( 0, 0, zPos )) ;
		for( int r = 0; r < vNumC; ++r)
		{
			m_verts.push_back( TVector3( rad * cos( stepR * r ),  rad * sin( stepR * r ), zPos) );
			
			//tetraを追加
			if( k != 0 )
			{
				p0 = (k-1) * (vNumC + 1);
				p3 =   k   * (vNumC + 1);
				if( r == 0 )
				{
					p1 = p0 + vNumC; p4 = p3 + vNumC;
					p2 = p0 + 1; 	 p5 = p3 + 1;
				}else{
					p1 = p0 + r;     p4 = p3 + r;
					p2 = p0 + r+1;   p5 = p3 + r+1;
				}
				m_tetras.push_back( TTetra(p0, p1, p2, p4) ); 
				m_tetras.push_back( TTetra(p0, p2, p5, p4) );
				m_tetras.push_back( TTetra(p0, p3, p4, p5) );
			}
		}
	}
}

void TTetraModel::initializeAsRotatedCylinder(int vNumC, int vNumL, double rad, double circleRad)
{
	TMatrix16 frame;
	TVector3 v;

	m_verts .clear();
	m_tetras.clear();

	const double stepR     = 2 * M_PI / ( (double)vNumC);
	const double stepAngle = 2 * M_PI / ( (double)vNumL + 5);

	int p0 = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0;
	
	for( int k = 0; k < vNumL; ++k)
	{
		frame.SetAsRotateY( -k * stepAngle );
		
		//軸上の点//
		v.Set( circleRad, 0, 0 );
		m_verts.push_back( frame * v);

		for( int r = 0; r < vNumC; ++r)
		{
			v.Set( rad * cos( stepR * r ) , rad * sin( stepR * r ), 0 );
			v.data[0] += + circleRad;
			m_verts.push_back( frame * v);
			
			//tetraを追加
			if( k != 0 )
			{
				p0 = (k-1) * (vNumC + 1);
				p3 =   k   * (vNumC + 1);
				if( r == 0 )
				{
					p1 = p0 + vNumC; p4 = p3 + vNumC;
					p2 = p0 + 1; 	 p5 = p3 + 1;
				}else{
					p1 = p0 + r;     p4 = p3 + r;
					p2 = p0 + r+1;   p5 = p3 + r+1;
				}
				m_tetras.push_back( TTetra(p0, p1, p2, p4) ); 
				m_tetras.push_back( TTetra(p0, p2, p5, p4) );
				m_tetras.push_back( TTetra(p0, p3, p4, p5) );
			}
		}
	}
}




void TTetraModel::initializeAsRotatedCylinderWithFlagellumBase( int vNumC, int vNumL, double rad, double circleRad, double rootRad)
{
	const double stepR      = 2 * M_PI / ( (double)vNumC);
	const double stepAngle  = 2 * M_PI / ( (double)vNumL);

	const int    rootThresh = (int)(vNumL * 0.08 );

	m_verts .clear();
	m_tetras.clear();
	
	vector<TVector3> axis;
	TMatrix16 frame;
	TVector3  v    ;

	//まず、axis pointsを生成
	for( int k = 0; k < vNumL; ++k)
	{
		frame.SetAsRotateY( -k * stepAngle );	
		v.Set( circleRad, 0, 0 );
		axis.push_back( frame * v );
		
		if( k < rootThresh ) 
		{
			//半径を狭める
			double l    = axis.back().Length();
			//double a    = (rootThresh - k) / (double) rootThresh;//だんだん０になる
			double a    = (rootThresh - k) * (rootThresh - k) / ((double) rootThresh * (double) rootThresh);//だんだん０になる
			double newL = a * rootRad + (1.0-a) * circleRad;
			axis.back() *= newL / l;
		}
	}

	int p0 = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0;
	TVector3 dir, Z(0,0,1);
	
	for( int k = 0; k < vNumL; ++k)
	{
		dir.Set(0,0,0);
		if( k != 0       ) dir.AddSubtract( axis[ k ], axis[k-1] );
		if( k != vNumL-1 ) dir.AddSubtract( axis[k+1], axis[ k ] );

		m_verts.push_back( axis[k] );
		
		
		t_getRotMat_V1ToV2( Z, dir, frame );
		

		//周囲の点と,
		for( int r = 0; r < vNumC; ++r)
		{
			v.Set( rad * cos( stepR * r ) , rad * sin( stepR * r ), 0 );
			v = frame * v;
			v += axis[k];
			m_verts.push_back( v );
			
			//tetraを追加
			if( k != 0 )
			{
				p0 = (k-1) * (vNumC + 1);
				p3 =   k   * (vNumC + 1);
				if( r == 0 )
				{
					p1 = p0 + vNumC; p4 = p3 + vNumC;
					p2 = p0 + 1; 	 p5 = p3 + 1;
				}else{
					p1 = p0 + r;     p4 = p3 + r;
					p2 = p0 + r+1;   p5 = p3 + r+1;
				}
				m_tetras.push_back( TTetra(p0, p1, p2, p4) ); 
				m_tetras.push_back( TTetra(p0, p2, p5, p4) );
				m_tetras.push_back( TTetra(p0, p3, p4, p5) );
			}
		}
	}
}


void TTetraModel::initializeSimpleTetraMesh()
{
	clear();
	m_verts.push_back( TVector3(0,0,0  ) );
	m_verts.push_back( TVector3( 2,3,-2) );
	m_verts.push_back( TVector3(-2,3,-2) );
	m_verts.push_back( TVector3( 0,3, 2) );
	m_verts.push_back( TVector3( 0,6,0 ) );
	m_tetras.push_back( TTetra( 0,1,2,3) );
	m_tetras.push_back( TTetra( 4,1,3,2) );

}






bool TTetraModel::loadfile_EleNode( const char *fname_noExt )
{
	clear();
	fprintf( stderr, "load tetMesh filepath = %s \n", fname_noExt);
	vector<TVector3> tmpVs;
	vector<TTetra    > tmpTs;

	string fnameNode( fname_noExt ); fnameNode += ".node";
	string fnameEle ( fname_noExt ); fnameEle  += ".ele";

	string line;
	ifstream ifs;
	
	ifs.open(fnameNode.c_str());
	if (!ifs.is_open()) { 
		fprintf(stderr, "Could not open file: %s\n", fnameNode.c_str()); 
		return false; 
	}
	else
	{
		int numVtx;  ifs >> numVtx;
		int tmp   ;  ifs >> tmp >> tmp >> tmp;
		tmpVs.resize( numVtx );
		for (int i = 0; i < numVtx; ++i) ifs >> tmp >> tmpVs[i].data[0] >> tmpVs[i].data[1] >> tmpVs[i].data[2];
		ifs.close();
	}

	ifs.open(fnameEle.c_str());
	if (!ifs.is_open())
	{
		fprintf(stderr, "Could not open file: %s\n", fnameEle.c_str()); 
		return false;
	}
	else
	{
		int numTet; ifs >> numTet    ;
		int tmp   ; ifs >> tmp >> tmp;
		tmpTs.resize(numTet);
		for (int i = 0; i < numTet; ++i) 
			ifs >> tmp >> tmpTs[i].v[0] >> tmpTs[i].v[1] >> tmpTs[i].v[2] >> tmpTs[i].v[3];
		ifs.close();
	}	
	m_verts = tmpVs;
	m_tetras= tmpTs;
	return true;
}




bool TTetraModel::savefile_EleNode( const char *fname_noExt )
{
	string fnameNode( fname_noExt ); fnameNode += ".node";
	string fnameEle ( fname_noExt ); fnameEle  += ".ele";

	string line;
	ofstream ofs;
	ofs.open(fnameNode.c_str(), ios::out);
	if (!ofs.is_open()) { 
		fprintf(stderr, "Could not open file: %s\n", fnameNode.c_str()); 
		return false; 
	}
	else
	{
		ofs << (int) m_verts.size() << " " << 3 << " " << 0 << " " << 0 << endl;
		for (int i = 0; i <(int)m_verts.size(); ++i) ofs << i << " " << m_verts[i].data[0] << " " 
			                                                         << m_verts[i].data[1] << " " 
																     << m_verts[i].data[2] << endl;
		ofs.close();
	}

	ofs.open(fnameEle.c_str());
	if (!ofs.is_open()){
		fprintf(stderr, "Could not open file: %s\n", fnameEle.c_str()); 
		return false;
	}else{
		ofs << (int) m_tetras.size()<< " " << 4 << " " << 0  << endl;
		for (int i = 0; i <(int)m_tetras.size(); ++i) ofs << i << " " << m_tetras[i].v[0] << " " << m_tetras[i].v[1] 
		                                                       << " " << m_tetras[i].v[2] << " " << m_tetras[i].v[3] << endl;
		ofs.close();
	}	
	return true;
}



/*
.ttet file (bit)
vertexNum tetNum
vertex
*/
bool TTetraModel::loadfile_ttet( const char *fname )
{
	FILE *fp;
	fp = fopen(fname,"rb");
	if( fp == 0 ) return false;

	clear();
	int tSize, vSize;
	fread( &vSize, sizeof(int), 1, fp);   
	fread( &tSize, sizeof(int), 1, fp);

	m_verts.resize( vSize ); 
	for( int i=0; i<vSize; ++i) fread( m_verts [i].data, sizeof(double), 3, fp);
	m_tetras.resize( tSize ); 
	for( int i=0; i<tSize; ++i) fread( m_tetras[i].v   , sizeof(int   ), 4, fp);

	fclose(fp);
	return true;
}

bool TTetraModel::savefile_ttet( const char *fname )
{
	FILE *fp;
	fp = fopen(fname,"wb");
	if( fp == 0 ) return false;

	int tSize = (int)m_tetras.size(), vSize = (int)m_verts.size();
	fwrite( &vSize, sizeof(int), 1, fp);   
	fwrite( &tSize, sizeof(int), 1, fp);
	for( int i=0; i<vSize; ++i) fwrite( m_verts [i].data, sizeof(double), 3, fp);
	for( int i=0; i<tSize; ++i) fwrite( m_tetras[i].v   , sizeof(int   ), 4, fp);
	fclose(fp);
	return true;
}




//x0x1x2は左まわりが正
inline static bool getCrossingPointRay2Triangle(const TVector3 &rayP, const TVector3 &rayD,  
	                                            const TVector3 &x0  , const TVector3 &x1, const TVector3 &x2, TVector3 &p)
{
	TMatrix9 M;
	TVector3 norm, b, stk;
	//norm.Set_V1subtV2_outmult_V3subtV4( x1, x0, x2, x0 );
	//norm.Normalize_Self();

	//if( norm * rayD >= 0) return false;

	M.Set( x1.data[0]-x0.data[0], x2.data[0]-x0.data[0], - rayD.data[0], 
		   x1.data[1]-x0.data[1], x2.data[1]-x0.data[1], - rayD.data[1], 
		   x1.data[2]-x0.data[2], x2.data[2]-x0.data[2], - rayD.data[2] );

	if( !M.getInvertSelf() ) return false;
	
	b.SetSubtract( rayP, x0);
	t_MatMultVec( M, b, stk);
	double s = stk.data[0];
	double t = stk.data[1];
	double k = stk.data[2];

	if (s < 0 || 1 < s || 
		t < 0 || 1 < t || 1 < s + t) return false;
	p.Set_V1_Add_CoefMultV2( rayP, k, rayD);

	return true;

}


bool TTetraModel::pickOneTetra( const TVector3 &rayP, const TVector3 &rayD, TVector3 &pos, int tetId)
{
	if( tetId                <  0     ) return false;
	if( (int)m_tetras.size() <= tetId ) return false;

	double minDist = 1000000000;

	TVector3 p;
	const TVector3 &x0 = m_verts[ m_tetras[tetId].v[0] ];
	const TVector3 &x1 = m_verts[ m_tetras[tetId].v[1] ];
	const TVector3 &x2 = m_verts[ m_tetras[tetId].v[2] ];
	const TVector3 &x3 = m_verts[ m_tetras[tetId].v[3] ];
	bool found = false;
	if( getCrossingPointRay2Triangle( rayP, rayD, x0, x2, x1, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); found = true;}
	if( getCrossingPointRay2Triangle( rayP, rayD, x0, x3, x2, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); found = true;}
	if( getCrossingPointRay2Triangle( rayP, rayD, x0, x1, x3, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); found = true;}
	if( getCrossingPointRay2Triangle( rayP, rayD, x1, x2, x3, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); found = true;}
	return found;
}


bool TTetraModel::pickAllTetra( const TVector3 &rayP, const TVector3 &rayD, TVector3 &pos, int &tetId) const
{
	tetId = -1;
	double minDist = DBL_MAX;

	TVector3 p;
	for( int i=0; i< (int) m_tetras.size(); ++i)
	{
		const TVector3 &x0 = m_verts[ m_tetras[i].v[0] ];
		const TVector3 &x1 = m_verts[ m_tetras[i].v[1] ];
		const TVector3 &x2 = m_verts[ m_tetras[i].v[2] ];
		const TVector3 &x3 = m_verts[ m_tetras[i].v[3] ];

		if( getCrossingPointRay2Triangle( rayP, rayD, x0, x2, x1, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); tetId = i; fprintf( stderr, "1  ");}
		if( getCrossingPointRay2Triangle( rayP, rayD, x0, x3, x2, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); tetId = i; fprintf( stderr, "2  ");}
		if( getCrossingPointRay2Triangle( rayP, rayD, x0, x1, x3, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); tetId = i; fprintf( stderr, "3  ");}
		if( getCrossingPointRay2Triangle( rayP, rayD, x1, x2, x3, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); tetId = i; fprintf( stderr, "4  ");}
	}
	return tetId != -1;
}


bool TTetraModel::pickAllTetra( const TVector3 &rayP, const TVector3 &rayD, const vector<byte> &ignorTetList, TVector3 &pos, int &tetId, double &minDist) const
{
	tetId   = -1;
	minDist = DBL_MAX;

	TVector3 p;
	for( int i=0; i< (int) m_tetras.size(); ++i) if( !ignorTetList[i] ) 
	{
		const TVector3 &x0 = m_verts[ m_tetras[i].v[0] ];
		const TVector3 &x1 = m_verts[ m_tetras[i].v[1] ];
		const TVector3 &x2 = m_verts[ m_tetras[i].v[2] ];
		const TVector3 &x3 = m_verts[ m_tetras[i].v[3] ];

		if( getCrossingPointRay2Triangle( rayP, rayD, x0, x2, x1, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); tetId = i; fprintf( stderr, "1  ");}
		if( getCrossingPointRay2Triangle( rayP, rayD, x0, x3, x2, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); tetId = i; fprintf( stderr, "2  ");}
		if( getCrossingPointRay2Triangle( rayP, rayD, x0, x1, x3, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); tetId = i; fprintf( stderr, "3  ");}
		if( getCrossingPointRay2Triangle( rayP, rayD, x1, x2, x3, p) && t_distance_sq( p, rayP) < minDist) { pos = p; minDist = t_distance_sq( p, rayP); tetId = i; fprintf( stderr, "4  ");}
	}
	return tetId != -1;
}




bool TTetraModel::do23Flip(int tId)
{
	return false;
}


bool TTetraModel::do32Flip(int tId, const vector<vector<int>> &e_neiborTets)
{
	static int oppoT1[6] = {2,3,1,3,1,2};
	static int oppoT2[6] = {3,1,2,0,0,0};


	TTetra &tet = m_tetras[tId];

	set<int> tetV; 
	tetV.insert( tet.v[0] );
	tetV.insert( tet.v[1] );
	tetV.insert( tet.v[2] );
	tetV.insert( tet.v[3] );

	for( int k=0; k<6; ++k) if( e_neiborTets[ tet.e[k] ].size() == 3 )
	{

		if( tet.n[oppoT1[k]] == -1 || tet.n[oppoT2[k]] == -1 ) continue;

		TTetra &ot1 = m_tetras[ tet.n[oppoT1[k]] ];  
		TTetra &ot2 = m_tetras[ tet.n[oppoT2[k]] ];  
		int oppoTet1V, oppoTet2V ;
		for( int v=0; v<4; ++v) if( tetV.find( ot1.v[v] ) == tetV.end() ) oppoTet1V = ot1.v[v];
		for( int v=0; v<4; ++v) if( tetV.find( ot2.v[v] ) == tetV.end() ) oppoTet2V = ot2.v[v];

	/*
		if( oppoTet1V != oppoTet2V ) {
			fprintf( stderr, "\na2 %d %d %d %d\n", ot1.n[0], ot1.n[1], ot1.n[2], ot1.n[3] );
			fprintf( stderr, "a2 %d %d %d %d\n"  , ot2.n[0], ot2.n[1], ot2.n[2], ot2.n[3] );
		}
		*/
		if( oppoTet1V != oppoTet2V ) continue;
		fprintf( stderr, "a1 ");

		TVector3 &p = m_verts[oppoTet1V];
		if( !t_sphereTest_bIn( m_verts[tet.v[0]], m_verts[tet.v[1]], m_verts[tet.v[2]], m_verts[tet.v[3]], p ) ) continue;


		//do 3-2 flip!
		fprintf( stderr, "found flip able tetra and edge!%d\n", tId);

	}

	return false;

}

void TTetraModel::faceFlipping(int maxNum)
{

	calcEdgesFromTetras();
	vector<vector<int> > e_neiTets( m_edges.size(), vector<int>());
	map<double,int> t_sortTets ;

	TVector3 center, v1,v2,v3;
	double     rad;
	for( int i=0; i<(int)m_tetras.size(); ++i)
	{
		TTetra &t = m_tetras[i];
		e_neiTets[ t.e[0] ].push_back( i ); e_neiTets[ t.e[1] ].push_back( i ); e_neiTets[ t.e[2] ].push_back( i );
		e_neiTets[ t.e[3] ].push_back( i ); e_neiTets[ t.e[4] ].push_back( i ); e_neiTets[ t.e[5] ].push_back( i );

		t_computeCircumSphere( m_verts[t.v[0]], m_verts[t.v[1]], m_verts[t.v[2]], m_verts[t.v[3]], center, rad);
		v1.SetSubtract( m_verts[ t.v[1] ], m_verts[ t.v[0] ] );
		v2.SetSubtract( m_verts[ t.v[2] ], m_verts[ t.v[0] ] );
		v3.SetSubtract( m_verts[ t.v[3] ], m_verts[ t.v[0] ] );
		double vol = t_V1crosV2_multV3( v2, v3, v1 );
		t_sortTets.insert( make_pair( vol/(rad*rad*rad), i) );
	}

	for( map<double, int>::iterator it=t_sortTets.begin(); it!=t_sortTets.end(); ++it)
	{
		bool canFlip = false;
		TTetra &t = m_tetras[ it->second ];

		//3-2 flip出来る？
		if( do32Flip( it->second, e_neiTets) ) 
		{
			canFlip = true;
			break;
		}


		//2-3 flipできる？
		


		if( canFlip ) break;
	}



}



void TTetraModel::laplacianSmoothing( int iterationN )
{
	vector< byte > v_canMove( m_verts.size(), true );
	vector<set<int>> v_1Ring( m_verts.size() );


	for( int i=0; i<(int) m_tetras.size(); ++i)
	{
		TTetra &t = m_tetras[i];
		if( t.n[0] == -1 ) v_canMove[t.v[1]] = v_canMove[t.v[2]] = v_canMove[t.v[3]] = false;
		if( t.n[1] == -1 ) v_canMove[t.v[0]] = v_canMove[t.v[2]] = v_canMove[t.v[3]] = false;
		if( t.n[2] == -1 ) v_canMove[t.v[0]] = v_canMove[t.v[1]] = v_canMove[t.v[3]] = false;
		if( t.n[3] == -1 ) v_canMove[t.v[0]] = v_canMove[t.v[1]] = v_canMove[t.v[2]] = false;
		
		for( int k=0; k<4; ++k)
		for( int l=0; l<4; ++l) if( k != l) v_1Ring[t.v[k]].insert( t.v[l]); 
	}

	for( int kkk=0; kkk < iterationN; ++kkk)
	{
		vector< TVector3 > newVs( m_verts.size() );
		TVector3 gc;
		for( int v = 0; v<(int) m_verts.size(); ++v)
		{
			if( ! v_canMove[v] ){ newVs[v] = m_verts[v]; continue;}

			gc.Set( 0,0,0 );
			for( set<int>::iterator it = v_1Ring[v].begin(); it != v_1Ring[v].end(); ++it) gc += m_verts[*it];
			gc /= (double) v_1Ring[v].size();

			newVs[v] = gc;
		}
		m_verts = newVs;

		checkStrangeTetra();
	}

}


//注意 too slow implimentation
int  TTetraModel::getEdgeIdx( const int &v0, const int &v1){
	for( int i=0,s=(int)m_edges.size(); i<s; ++i)
	{
		if( m_edges[i].v[0] == v0 && m_edges[i].v[1] == v1 ) return i;
		if( m_edges[i].v[0] == v1 && m_edges[i].v[1] == v0 ) return i;
	}
	return -1;
}
