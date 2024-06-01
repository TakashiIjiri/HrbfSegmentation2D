#ifndef __T_TRIANGLEMESH_H_INCLUDED__
#define __T_TRIANGLEMESH_H_INCLUDED__

#include "TOGL.h"
#include "tmath.h"
#include <set>

/*TWingEdge-------------------------------------------------------
			      \    /
			  e[2] \  / e[0]
				    \/
				    . v[0] ( < v[1] )
				    |
			 p[1]   |    p[0]
				    |
				    . v[1]
				    /\
			  e[3] /  \ e[1]
			      /    \  
----------------------------------------------------------------*/	
class TWingEdge 
{     
public:
	int e[4] ;// winged edge id
	int p[2] ;// polygon IDs
	int v[2] ;// vertex ids v[0] < v[1]
	TWingEdge(){ e[0]=e[1]=e[2]=e[3]=p[0]=p[1]=v[0]=v[1]=-1 ; }
	TWingEdge(const TWingEdge& we){ Set(we) ;}
	void  Set(const TWingEdge& we){ e[0] = we.e[0];  e[1] = we.e[1] ; e[2] = we.e[2] ; e[3] = we.e[3] ; 
		                            p[0] = we.p[0];  p[1] = we.p[1] ; v[0] = we.v[0] ; v[1] = we.v[1] ;}
	TWingEdge operator=(const TWingEdge& we){ Set(we); return *this ;	}
	bool ReplacePolyID( int src,int tgt ){
		if     (p[0]==src) p[0]=tgt;
		else if(p[1]==src) p[1]=tgt;
		else return false;
		return true ;
	}	
	bool ReplaceVtxID( int src,int tgt ){
		if     (v[0]==src) v[0]=tgt;
		else if(v[1]==src) v[1]=tgt;
		else return false;
		return true ;
	}
	bool ReplaceEdgeID( int src,int tgt ){
		if     (e[0]==src) e[0]=tgt;
		else if(e[1]==src) e[1]=tgt;
		else if(e[2]==src) e[2]=tgt;
		else if(e[3]==src) e[3]=tgt;
		else return false;
		return true ;
	}

	inline void invert(){ swap(v[0], v[1]); swap(p[0], p[1]);
						  swap(e[0], e[3]); swap(e[1], e[2]);}
	inline void FlipLR(){ swap(p[0], p[1]); swap(e[0], e[2]); swap(e[1], e[3]);}
	inline bool isBound()const{ return (p[0] == -1 || p[1] == -1); }

	inline void trace(){
		fprintf( stderr, "v[0/1] %d %d\n ", v[0], v[1]);
		fprintf( stderr, "p[0/1] %d %d\n ", p[0], p[1]);
		fprintf( stderr, "e[0/1/2/3] %d %d %d %d\n ", e[0], e[1], e[2], e[3]);
	}

	inline bool hasV0V1(const int &v0, const int v1)const{
		return (v[0]==v0 && v[1]==v1) || 
			   (v[0]==v1 && v[1]==v0);
	}
} ; 

class TTriangle 
{
public:
	int idx[3];//vertices or edges IDs  stored in the counter-clockwise
	
	TTriangle( int v0=0, int v1=0, int v2=0  ){ Set(v0,v1,v2); }
	TTriangle( const TTriangle& p ){ Set(p) ; }
	TTriangle& operator=( const TTriangle& p ){	Set(p) ; return *this ;	}

	inline void Set(int v0,int v1,int v2){idx[0] =      v0 ; idx[1] =      v1 ; idx[2] =      v2 ; }
	inline void Set(const int v[3]      ){idx[0] =     v[0]; idx[1] =     v[1]; idx[2] =     v[2]; }
	inline void Set(const TTriangle &p  ){idx[0] = p.idx[0]; idx[1] = p.idx[1]; idx[2] = p.idx[2]; }

	inline void Flip(int i0, int i1) { swap( idx[i0], idx[i1] ); }
	bool ReplaceID(int src,int tgt){
		if     (src == idx[0]) idx[0] = tgt;
		else if(src == idx[1]) idx[1] = tgt;
		else if(src == idx[2]) idx[2] = tgt;
		else return false;
		return true;
	}
	void trace( ){fprintf( stderr, "%d %d %d\n", idx[0], idx[1], idx[2] );}
} ;

/*----Useage-------------------------------------------------------------
1 + initial contruction by a obj file  or  vertices and polygons 
	+ verts polys and normal will be constructed
2 + construct edges   by updateWingEdge function *optional
3 + construct oneRing by updateOneRing  function *optional
4 + update normal     by updateNormal   function when the model is deformed. 
The user is supporsed NOT to modify vertices/polygons/edges directry from outside of this class//
------------------------------------------------------------------------*/
class TTriangleMesh  
{
protected:
	int m_vSize, m_eSize, m_pSize;
public:
	TVector3  *m_verts  ;// constructed by init* functions	
	TVector3  *m_v_norm ;// constructed by init* functions  &  updated by updateNormal function
	TVector3  *m_v_uv   ;

	TTriangle *m_polys  ;// constructed by init* functions 
	TVector3  *m_p_norm ;// constructed by init* functions  &  updated by updateNormal function

	TWingEdge  *m_edges  ;//constructed by CalcWingEdgeFromPolygon() method
	int        *m_v_edge ;//constructed by CalcWingEdgeFromPolygon() method	(reference for one wing edge from this verted)
	TTriangle  *m_p_edges;//constructed by CalcWingEdgeFromPolygon() method
	vector<vector<int>> m_v_1ring;//constructed by updateOneRing() method //各頂点のone ring neighbor vertices(sizeの保持が面倒なのでvector利用)

	static int MCsEdgeTable[256]   ;//Marching cubes用table
	static int MCsTriTable[256][16];//Marching cubes用table

	bool       m_WingEdgeIsAvailable;

public:   
	TTriangleMesh() ;
	TTriangleMesh(const TTriangleMesh &src);
	virtual ~TTriangleMesh();
	virtual void Clear();
	inline int getVnum()const{return m_vSize;}
	inline int getEnum()const{return m_eSize;}
	inline int getPnum()const{return m_pSize;}
	inline TTriangleMesh& operator=(const TTriangleMesh& src){ init( src ); return *this; }
	//initialization  and  construction of normal / edge / onering 
	void init              (const TTriangleMesh &pm);
	bool initFromVertsPolys(const vector<TVector3> &verts, const vector<TTriangle > &polys);
	bool initFromObjFile   (const char* fname);
	void initAsIcosahedron (double r );
	void initAsSphere      (double r, int resH, int resV);
	void initAsCube        (double sizeX, double sizeY, double sizeZ);
	bool initAs2DGridPatch(const vector<TVector3> bPoints_2D, double gridR);
	bool initAs2DGridPatch_boundVCorrespond(const vector<TVector3> bPoints_2D, double gridR, vector<int> &mapPointsToVerts);

	void updateNormal  ();// this construct "m_p_norm / m_v_norm" from "m_verts / m_polys"  

	//WingedEdge operation////////////////
	bool checkWingedEdge();
	bool checkVtxMoreThan4BoundEdge();//more than boundary edges / no polygons
	void GetVsPsEsAroundVertex(const int vID, vector<int> &Verts, vector<int> &Polys, vector<int> &Edges) const;
	void GetPsAroundVertex(const int vID, vector<int> &Polys) const;
	void GetVsAroundVertex(const int vID, vector<int> &Verts) const;
	void GetEsAroundVertex(const int vID, vector<int> &Edges) const;
	int  GetDegreeOfVertex(const int vID);		
	bool bBoundaryVertex  (const int vID)const;		
	//subdivision////////////////////////////////////////////////////////////////////////////////////////
	void subdivisionRoot3(const vector<byte> &trgtPs);
	void subdivisionRoot3(const        byte  *trgtPs);
	void subdivisionRoot3();
	bool subdivision_divideEdges( const vector<byte> e_bSubdiv, const vector<TVector3> &e_vPos );

	//Smoothing/////////////////////////////////////////////////////////////////////////////////////////
	void SmoothingByLaplacianNormal(int time);
	void SmoothingByLaplaceMagDiff (int time);
	void SmoothingByLaplacianZero  (int time);
	void SmoothingByLaplacianZero_fixBound     (int time);
	void SmoothingByLaplacianZero_OnlyShapeEdge(int time = 1, double dotThresh = 0.01);
	void SmoothingAtOneVertex      (int vIdx);
	//Topology modifications/////////////////////////////////////////////////////////////////////////////
	void splitEdge(int trgtEdgeId, const TVector3 &newVertexPos, int &newVertexId);
	void RemoveUnUsedVertices();
	bool RemoveEdge        (int edgeid, bool topologyOnly = false );
	bool RemoveEdge_onBound(int edgeId );
	void RemoveShortEdges( const vector<short> &v_removable, vector< int >  &v_mapPreToNew, double edgeThreshRate);
	//region growing and labeling///////////////////////////////////////////////////////////////////////////////////////////
	int  calcLinkedLabel     ( vector<int> &v_labels );
	int  calcLinkedLabel_poly( vector<int> &p_labels );
	void calcRoiMesh         ( const vector<int> &seedVerts, int roiSize, vector<TVector3> &roiVs, vector<TTriangle> &roiPs, 
		                                                                  vector<int> &vmap_ROItoFUL, vector<int> &vmap_FULtoROI);
	void getRingNeighborhood( int vIdx, int ringSize, set<int> &ringVerts);
	bool calcNearestPolygon_insideCheck(const TVector3 &p, TVector3 &nearPos, int &polyIdx,  double &dist) const;
	//Hole filling///////////////////////////////////////////////////////////////////////////////////////
	bool GetSequentialBoundVsEs( int seedVtx, vector<int> &seqVs, vector<int> &seqEs);
	bool fillHole ( int seedVtx );
	bool fillHoles(  );
	bool fillHole_noExtraVertex(int seedIdx);
	bool fillHoles_noExtraVertex();
	//MISCS///////////////////////////////
	void   drawPolygons   () const;
	void   drawPolygons_flatShading() const;
	void   drawEdges      ();
	void   drawBoundEdges ();
	void   drawBoundPoints();
	void   trace          ();
	void   SaveVsPsAsText(FILE *fp);
	bool   LoadVsPsAsText(FILE *fp);
	bool   pickSurface                 (const TVector3 &rayP, const TVector3 &rayD, TVector3 &pos, int &polyID, int bothFrontBack = 0/*0both 1front 2back*/) const;
	bool   getCrossingPolygonToRay     (const TVector3 &rayP, const TVector3 &rayD, int &polyID, TVector3 &pos, int bothFrontBack = 0/*0both 1front 2back*/) const;
	bool   getCrossingPolygonToLineSegm(const TVector3 &p0  , const TVector3 &p1  , int &polyID, TVector3 &pos, double &distInP0P1) const;
	double   CalcAverageEdgeLength();
	TVector3 CalcAverageNormal    ();

	void   sortVerts_byXvalue();//ランダムアクセスが多少改善される
	void   sortPolys_byXvalue();//ランダムアクセスが多少改善され
	inline int searchEdgeId( int vId1, int vId2) const{
		if( vId1 > vId2 ) swap( vId1, vId2);
		vector<int> Es; GetEsAroundVertex( vId1, Es );
		for( int i=0, s = (int)Es.size(); i<s; ++i)  if( m_edges[ Es[i] ].v[0] == vId1 && m_edges[ Es[i] ].v[1] == vId2 ) return Es[i];
		return -1;
	}
	inline void ComputeBBox(TVector3& minB,TVector3& maxB) const { t_verts_calcBoundary( m_vSize, m_verts, minB, maxB); }
	inline void TranslateAllVertices(const TVector3 &shift ){ for( int i=0; i<m_vSize; ++i) m_verts[i] += shift; }
	inline void ScaleAllVertices    (const double     &scale ){ for( int i=0; i<m_vSize; ++i) m_verts[i] *= scale; }
	inline void TransformAllVertices(const TMatrix16 &M     ){ for( int i=0; i<m_vSize; ++i) t_MatMultVec( m_verts[i], M);}

	inline void Flip          (){ for( int i=0; i<m_pSize; ++i){ m_polys  [i].Flip(1,2);
	                                                             m_p_edges[i].Flip(0,2);
																 m_p_norm [i] *= -1; }
		                          for( int i=0; i<m_eSize; ++i)  m_edges  [i].FlipLR(); 
		                          for( int i=0; i<m_vSize; ++i)  m_v_norm [i] *= -1; }
	inline void FlipNormalOnly(){ for( int i=0; i<m_pSize; ++i)  m_p_norm[i] *= -1;
		                          for( int i=0; i<m_vSize; ++i)  m_v_norm[i] *= -1;}
	inline void CalcNormalOfOnePolygon(int pol_id, TVector3 &norm) const {
		int* pVtx = m_polys[pol_id].idx; 
		norm.Set_V1subtV2_crsprd_V3subtV4( m_verts[ pVtx[1] ], m_verts[ pVtx[0] ],m_verts[ pVtx[2] ], m_verts[ pVtx[0] ]);
		norm.Normalize_Self();
	}
	inline double GetPolygonArea(int pol_id){
		int* pVtx = m_polys[pol_id].idx; TVector3 n;
		n.Set_V1subtV2_crsprd_V3subtV4( m_verts[ pVtx[1] ], m_verts[ pVtx[0] ],m_verts[ pVtx[2] ], m_verts[ pVtx[0] ]);
		return 0.5 * n.Length();
	}
	inline double GetSurfaceArea(){ double ret=0;   for( int i=0; i<m_pSize; ++i) ret += GetPolygonArea(i);   return ret ; }
	inline double getNearestVertexIdx( const TVector3 &p, int &vIdx )const{ double dist; t_verts_getNearestPoint( m_vSize, m_verts, p, vIdx, dist); return dist;}	
	inline void   GetVerts(vector<TVector3> &trgt)const{trgt.resize( m_vSize ); for( int i=0; i<m_vSize; ++i) trgt[i]=m_verts[i]; }
	inline void   GetPolys(vector<TTriangle > &trgt)const{trgt.resize( m_pSize ); for( int i=0; i<m_pSize; ++i) trgt[i]=m_polys[i]; }
	inline void   GetEdges(vector<TWingEdge > &trgt)const{trgt.resize( m_eSize ); for( int i=0; i<m_eSize; ++i) trgt[i]=m_edges[i]; }

protected:
	void pushbackNewEdge( int num = 1 );//m_edges       の末尾に初期化された edgeをnum個追加する
	void pushbackNewPoly( int num = 1 );//m_polysとm_p_*の末尾に初期化された edgeをnum個追加する
	void pushbackNewVert( int num = 1 );//m_vertsとm_v_*の末尾に初期化された edgeをnum個追加する
	bool RemoveEdgePossible(int edgeid )const ;
	void RemoveVsPsEs   (const set<int> &Vs, const set<int> &Ps, const set<int> &Es ) ;
	void updateWingEdge();// this construct "m_edges / m_p_edges / m_v_edge / m_v_1ring  from "m_verts / m_polys"  
	void updateOneRing ();
};

bool t_intersect_RayToPolygons     (const TVector3 &rayP, const TVector3 &rayD, const TTriangleMesh &mesh, int &polyID, TVector3 &pos , int bothFrontBack);
bool t_intersect_LineSegmToPolygons(const TVector3 &p0  , const TVector3 &p1  , const TTriangleMesh &mesh, int &polyID, TVector3 &pos, double &distInP0P1);

////////////////////////////////////////////////////////////////////////////////
//global functions related to ILPOlygonModel
void t_drawSphere      (double r, const TVector3 &p);
void t_drawSphereEdges (double r);
void t_drawSphere_norm (double r, const TVector3 &p);
void t_drawSphere_norm (double r, const TVector3 &p         , const float ambi[4], const float diff[4], const float spec[4], const float shin[1]);
void t_drawSpheres_norm(double r, const vector<TVector3> &Ps, const float ambi[4], const float diff[4], const float spec[4], const float shin[1]);
///////////////////////////////////////////////////////////////////////////////////
////Marching cubes/////////////////////////////////////////////////////////////////
/*   vol[idx] = true or false --> inside or outside*/
void t_MartchingCubes_BinVolue             ( int W, int H, int D, const byte   *vol,                                    const TVector3 &cubeSize, TTriangleMesh &poly, const int *minIdx=0, const int *maxIdx=0);
void t_MartchingCubes_floatVolume          ( int W, int H, int D, const float  *vol, float thresh                     , const TVector3 &cubeSize, TTriangleMesh &poly, const int *minIdx=0, const int *maxIdx=0);
void t_MartchingCubes_floatVolume          ( int W, int H, int D, const float  *vol, float thresh, const float *vol_InOut , const TVector3 &cubeSize, TTriangleMesh &poly);
void t_GenMeshFromBinaryVolume             ( int W, int H, int D, const byte   *vol                                   , const TVector3 &cubeSize, TTriangleMesh &poly);
void t_GenMeshFromBinaryVolume_coarseResNew( int W, int H, int D, const byte   *vol_mask, int trgtMaskId              , const TVector3 &cubeSize, TTriangleMesh &poly, int resolution = 3);
void t_GenMeshFromBinaryVolume_coarseResNew( int W, int H, int D, const byte   *vol_mask, const vector<int> &maskOnOff, const TVector3 &cubeSize, TTriangleMesh &poly, int resolution = 3 );
void t_GenMeshFromBinaryVolume_coarseResNew( int W, int H, int D,       byte   *volInOut, const TVector3 &volSize                               , TTriangleMesh &poly, int resolution = 3 );
void t_computeInternalVoxel                ( const int W, const int H, const int D, TTriangleMesh &mesh, const TVector3     &volCubeSize,  byte *bInside );//bInsideはnewされている

///////////////////////////////////////////////////////////////////////////////////
////MISCS//////////////////////////////////////////////////////////////////////////
void t_createCylinderSolid       (              double radius, double sizeY, int vertSizeR, int vertSizeY, TTriangleMesh &mesh);
void t_createRotatedCylinderSolid( double gRad, double radius,               int vertSizeR, int vertSizeY, TTriangleMesh &mesh);








#endif //__T_TRIANGLEMESH_H_INCLUDED__