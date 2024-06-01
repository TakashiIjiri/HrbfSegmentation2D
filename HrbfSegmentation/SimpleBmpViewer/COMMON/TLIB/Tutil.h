#ifndef __TUTIL_H_INCLUDED__
#define __TUTIL_H_INCLUDED__


/*-----------------------------------------------
TUtil.h is a high level function library.
Written by Takashi Ijiri @ riken / Oct. 2011
This library depends on TMath.h TOGL TTriangleMesh, and TOGLImages
------------------------------------------------*/

#include "tmath.h"
#include "TOGL.h"
#include "TTriangleMesh.h"
#include "TTetraModel.h"
#include "TOGLImage.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//CUT surface by 2D stroke/////////////////////////////////////////////////////////////////////////////////////

//t_computeCrossingEdgesToStroke/////////////////////////////////////////////////////////////
//compute cross sections between a Triangle Mesh and a 2D stroke (drawn on screen)

//const vector<TVector3>   &stroke2D, (input) : drawn 2D stroke  stroke2D[i].data[0] = 0 
//const TTriangleMesh      &mesh    , (input) : target 3D triangle mesh
//const TOGL               &ogl     , (input) : current open gl instance
//vector<int>              &e_isCrossing     , (output) : size is equal to mesh.m_eSize, true if the edge has intersection //そのエッジが交点を持つかどうか
//vector<TVector3>         &e_crossingPos    , (output) : size is equal to mesh.m_eSize, intersect position                //交点を持つ場合の交点の座標,
//vector<double>           &e_crsPosOnStroke , (output) : size is equal to mesh.m_eSize, relative position on the 2D stroke//交点ノ対応する2Dストローク上の位置

//edgeと2Dstrokeの交差点をすべて発見する
//計算はscreenに射影した2D上で行う
inline void t_computeCrossingEdgesToStroke( const vector<TVector3> &stroke2D, 
	                                 const TTriangleMesh    &mesh    ,
									 const TOGL             &togl     ,
								           vector<int>      &e_isCrossing ,   //そのエッジが交点を持つかどうか
								           vector<TVector3> &e_crossingPos,   //交点を持つ場合の交点の座標,
								           vector<double>   &e_crsPosOnStroke)//交点ノ対応する2Dストローク上の位置
{
	const TVector3 eyeP = togl.getEyePoint();
	const TVector3 eyeF = togl.getEyeFocus();
	const TVector3 eDir = (eyeF - eyeP).Normalize(); 

	//mesh vertex をスクリーンへ射影
	vector<TVector3> verts2D;
	t_projectVerticesOn2D( mesh.getVnum(), mesh.m_verts, verts2D, togl );

	vector<double> strokePtLen( stroke2D.size() );
	double len = 0;
	for( int i=0; i<(int) stroke2D.size()-1; ++i){ strokePtLen[i] = len;  len += t_distance( stroke2D[i], stroke2D[i+1] ); }
	strokePtLen.back() = len;

	const int vSize = mesh.getVnum(); const TVector3 *verts = mesh.m_verts;
	const int eSize = mesh.getEnum(); const TWingEdge  *edges = mesh.m_edges;

	e_isCrossing    .resize( eSize, 0);
	e_crossingPos   .resize( eSize   );
	e_crsPosOnStroke.resize( eSize, 0);
	
	//stroke v.s. edgeの交差判定
	TVector3 vD, strD, dir;
	for( int i=0; i<eSize; ++i)
	{
		//視点より後ろのedgeは飛ばす
		if( t_V1subtV2_multV3( verts[ edges[i].v[0] ], eyeP, eDir) < 0 ) continue;
		if( t_V1subtV2_multV3( verts[ edges[i].v[1] ], eyeP, eDir) < 0 ) continue; //(vP1 - eye) * eDir < 0

		const TVector3 &vP0 = verts2D[ edges[i].v[0] ];
		const TVector3 &vP1 = verts2D[ edges[i].v[1] ];
		vD.SetSubtract( vP1, vP0 );

		for( int kk=1; kk<(int) stroke2D.size(); ++kk)
		{
			const TVector3 &strP0 = stroke2D[kk-1];
			const TVector3 &strP1 = stroke2D[ kk ];
			strD.SetSubtract( strP1, strP0 );
			if( !t_intersect2D_Lines_preCheck(vP0, vP1, strP0, strP1) ) continue;
			
			//h = strP0 + s*strD = v0 + t*vD
			//(strD -vD)(s,t) = -StrP0 + vP0
			double s, t;
			if( t_solve2by2LinearEquation(   strD.data[0], -vD.data[0], strD.data[1], -vD.data[1],
										   -strP0.data[0] + vP0.data[0], -strP0.data[1] + vP0.data[1], s,t) &&  
										    0<=s && s<= 1 && 0< t && t<  1)
			{
				e_isCrossing[i] = 1;
				dir.SetSubtract( verts[edges[i].v[1]] , verts[edges[i].v[0]]);
				e_crossingPos[i].Set_V1_Add_CoefMultV2( verts[edges[i].v[0]], t, dir);

				e_crsPosOnStroke[i] = strokePtLen[kk-1] + s * strD.Length();//h = strP0 + s * strD
				break;
			}
		}
	}
}


//t_computeStripEdgeSeq///////////////////////////////////////////////////////////////
//compute sequential order of target edges 
//const TTriangleMesh     &mesh     ,(input): terget mesh
//const vector<int>       &e_bTerget,(input): true : target edge
//vector<vector<int>>     &e_strips ,(output) : multiple edge strips 
inline bool t_computeEdgeStrip( const TTriangleMesh     &mesh,
						 const vector<int>       &e_bTerget,
							vector<vector<int>>  &e_strips)
{
	const int eSize = mesh.getEnum(); const TWingEdge* edges   = mesh.m_edges  ;
	const int pSize = mesh.getPnum(); const TTriangle* p_edges = mesh.m_p_edges;
	                                  const TTriangle* polys   = mesh.m_polys  ;

	vector<int> e_flag( eSize, 0);
	for( int i=0;i< eSize; ++i )if( e_bTerget[i] ) e_flag[i] = 1;

	for( int i=0; i<pSize;++i) if( e_flag[ p_edges[i].idx[0] ] && 
		                           e_flag[ p_edges[i].idx[1] ] && 
								   e_flag[ p_edges[i].idx[2] ] ) return false;

	while(true)
	{	
		int pivEdgeId = -1;
		int prePoly   = -1;
		for( int i=0; i<pSize;++i)
		{
			bool e0 = e_flag[p_edges[i].idx[0]] == 1;
			bool e1 = e_flag[p_edges[i].idx[1]] == 1;
			bool e2 = e_flag[p_edges[i].idx[2]] == 1;
			if(  e0 && !e1 && !e2 ) {pivEdgeId = p_edges[i].idx[0]; prePoly = i; break;}
			if( !e0 &&  e1 && !e2 ) {pivEdgeId = p_edges[i].idx[1]; prePoly = i; break;}
			if( !e0 && !e1 &&  e2 ) {pivEdgeId = p_edges[i].idx[2]; prePoly = i; break;}
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
			const TTriangle &pEdges = p_edges[ pivPoly ];
			if     ( e_flag[ pEdges.idx[0] ] == 1 ) pivEdgeId = pEdges.idx[0];
			else if( e_flag[ pEdges.idx[1] ] == 1 ) pivEdgeId = pEdges.idx[1];
			else if( e_flag[ pEdges.idx[2] ] == 1 ) pivEdgeId = pEdges.idx[2];
			else break;//fiberが途中で切れる
			
			e_strips.back().push_back( pivEdgeId );
			e_flag[pivEdgeId] = 2;
			prePoly = pivPoly;
		}
	}
	return true;
}



//Dikstra method for 2D grid image///////////////////////
// gridVerts: 3D height field vertices aliged in W*H grid
// ringSize : neighborhood ring size
// distField:
inline void t_computeDistField_Dikstra(const TVector3 *gridVerts,
								const int W    , const int H,
								const int seedX, const int seedY, 
								const int ringSize, 
								float* distField/*allocated already size is width*height*/)
{
	const int seedVtx = seedX + W * seedY;

	set< int > fronteer                             ;//set<int>で現在の境界のset
	for(int i=0; i< W*H; ++i) distField[i] = FLT_MAX; //距離自体は，distFieldにそのまま入れる	


	fronteer.insert( seedVtx );
	distField[ seedVtx ] = 0;

	while( !fronteer.empty() )
	{
		//最小頂点の検索
		double minDist = DBL_MAX;
		int    trgtIdx = - 1    ;

		for( set<int>::iterator it = fronteer.begin(); it != fronteer.end(); ++it)
			if( distField[ *it ] < minDist ) { minDist = distField[ *it ]; trgtIdx = *it;}

		//値の確定，fronteerから削除しone ringをfronteerに挿入 //ringSize * ringSize windowを見る
		const int pivY = trgtIdx / W;
		const int pivX = trgtIdx - pivY*W;
		for( int y = -ringSize; y <= ringSize; ++y) if( 0 <= y+pivY && y+pivY < H) 
		for( int x = -ringSize; x <= ringSize; ++x) if( 0 <= x+pivX && x+pivX < W)
		{
			int idx = pivX+x + (pivY+y)*W;
			double d = minDist + t_distance( gridVerts[ trgtIdx], gridVerts[ idx ] );
			if( d < distField[ idx ] ){
				fronteer.insert( idx );
				distField[ idx ] = (float)d;
			}
		}
		fronteer.erase( fronteer.find( trgtIdx ) );
	}
}


inline void t_renderBoundaryTetras( const vector<TTetra> tetras, const vector<TVector3> verts )
{
	TVector3 n0,n1, n2,n3;
	glBegin( GL_TRIANGLES );
	for( int i=0, s = (int)tetras.size(); i< s; ++i)
	{
		const TTetra &t = tetras[i];
		const TVector3 &v0 = verts[ t.v[0] ];
		const TVector3 &v1 = verts[ t.v[1] ];
		const TVector3 &v2 = verts[ t.v[2] ];
		const TVector3 &v3 = verts[ t.v[3] ];
		n0.Set_V1subtV2_crsprd_V3subtV4( v2,v1,v3,v1); n0.Normalize_Self();
		n1.Set_V1subtV2_crsprd_V3subtV4( v3,v0,v2,v0); n1.Normalize_Self();
		n2.Set_V1subtV2_crsprd_V3subtV4( v0,v3,v1,v3); n2.Normalize_Self();
		n3.Set_V1subtV2_crsprd_V3subtV4( v1,v2,v0,v1); n3.Normalize_Self();
		if( t.n[0] == -1 ){ glNormal3dv( n0.data ); glVertex3dv( v1.data );  glVertex3dv( v2.data );  glVertex3dv( v3.data ); }
		if( t.n[1] == -1 ){ glNormal3dv( n1.data ); glVertex3dv( v0.data );  glVertex3dv( v3.data );  glVertex3dv( v2.data ); }
		if( t.n[2] == -1 ){ glNormal3dv( n2.data ); glVertex3dv( v3.data );  glVertex3dv( v0.data );  glVertex3dv( v1.data ); }
		if( t.n[3] == -1 ){ glNormal3dv( n3.data ); glVertex3dv( v2.data );  glVertex3dv( v1.data );  glVertex3dv( v0.data ); }
	}
	glEnd();	
}

#endif	//__TUTIL_H_INCLUDED__