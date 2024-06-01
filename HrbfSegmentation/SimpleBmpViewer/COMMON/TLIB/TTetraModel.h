#pragma once

#include "tmath.h"
#include "TTriangleMesh.h"


/*
ILTetraModelの改良
+edgeデータ追加
+subdivision対応

v1-v2-v3は外から見て左回り

        v0  e1   v2
         *-------*
        / \    . |
       /   \ .   |
    e0/    .\    |
     /   .  e2   |
    / .       \  |
   /.          \ |
  *_____________\*
  v1             v3

        v0       v2
         *-------*
        / \    . |
       /   \ .   |
      /    .\    |e4
     /  e3   \   | 
    / .       \  |
   /.          \ |
  *_____________\*
  v1     e5      v3
*/



class TTetra
{	
public:
	int v[4];//vertices
	int n[4];//neighbors
	int e[6];//edges
	TVector3 norm[4];

	TTetra(int v0=-1, int v1=-1, int v2=-1, int v3=-1) {
		v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
		n[0] = n[1] = n[2] = n[3] = -1;
		e[0] = e[1] = e[2] = e[3] = e[4] = e[5] = -1; 
	}

	inline void Set( int v0, int v1, int v2, int v3){
		v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
	}
	inline void Set( const TTetra &t ){
		memcpy( v, t.v, sizeof(int)*4); 
		memcpy( n, t.n, sizeof(int)*4); 
		memcpy( e, t.e, sizeof(int)*6); 
		norm[0].Set(t.norm[0]); norm[1].Set(t.norm[1]);
		norm[2].Set(t.norm[2]); norm[3].Set(t.norm[3]);
	}
	TTetra                  (const TTetra &src){ Set( src ); }
	inline TTetra& operator=(const TTetra &src){ Set( src ); return *this; }

	inline void rotR(int pivV)
	{
		static const int vorder[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
		static const int eorder[4][6] = {{0,1,2, 3,4,5}, {0,5,3, 2,4,1}, {1,3,4, 0,5,2}, {2,4,5, 1,3,0}};
		const int *vd = vorder[pivV];
		const int *ed = eorder[pivV];

		swap( v[ vd[2] ], v[ vd[1] ]);   swap( v[ vd[1] ], v[ vd[0] ]); 
		swap( n[ vd[2] ], n[ vd[1] ]);   swap( n[ vd[1] ], n[ vd[0] ]); swap( norm[ vd[2] ], norm[ vd[1] ]);   swap( norm[ vd[1] ], norm[ vd[0] ]); 
		swap( e[ ed[2] ], e[ ed[1] ]);   swap( e[ ed[1] ], e[ ed[0] ]); 
		swap( e[ ed[5] ], e[ ed[4] ]);   swap( e[ ed[4] ], e[ ed[3] ]);
	}
	inline void rotL(int pivV)
	{
		static const int vorder[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
		static const int eorder[4][6] = {{0,1,2, 3,4,5}, {0,5,3, 2,4,1}, {1,3,4, 0,5,2}, {2,4,5, 1,3,0}};
		const int *vd = vorder[pivV];
		const int *ed = eorder[pivV];
		swap( v[ vd[0] ], v[ vd[1] ]);   swap( v[ vd[1] ], v[ vd[2] ]); 
		swap( n[ vd[0] ], n[ vd[1] ]);   swap( n[ vd[1] ], n[ vd[2] ]); swap( norm[ vd[0] ], norm[ vd[1] ]);   swap( norm[ vd[1] ], norm[ vd[2] ]); 
		swap( e[ ed[0] ], e[ ed[1] ]);   swap( e[ ed[1] ], e[ ed[2] ]); 
		swap( e[ ed[3] ], e[ ed[4] ]);   swap( e[ ed[4] ], e[ ed[5] ]);
	}

	inline double calcVol(vector<TVector3> &Vs)
	{
		TVector3 v1,v2,v3;
		v1.SetSubtract( Vs[ v[1] ], Vs[ v[0] ] );
		v2.SetSubtract( Vs[ v[2] ], Vs[ v[0] ] );
		v3.SetSubtract( Vs[ v[3] ], Vs[ v[0] ] );
		return t_V1crosV2_multV3( v2, v3, v1 ) / 6.0;
	}
	inline bool isOnBoundary(){
		return n[0]==-1 || n[1] == -1 || 
			   n[2]==-1 || n[3] == -1;
	}


};


class TTetEdge
{
public:
	int v[2];//vertices  v0 < v1

	TTetEdge(int v0=-1, int v1=-1) { v[0] = v0; v[1] = v1; }

	inline void Set           (const TTetEdge &t ){ memcpy( v, t.v, sizeof(int)*2); }
	TTetEdge                  (const TTetEdge &src){ Set( src ); }
	inline TTetEdge& operator=(const TTetEdge &src){ Set( src ); return *this; }
};



class TTetBaryCoord
{
public:
	int     vIdx[4];
	double  coef[4];
	TTetBaryCoord(void){
		vIdx[0] = vIdx[1] = vIdx[2] = vIdx[3] = 0;
		coef[0] = coef[1] = coef[2] = coef[3] = 0;
	}
	TTetBaryCoord( const TTetBaryCoord &src){ Set( src ); }

	inline void Set( const TTetBaryCoord &src ){
		memcpy( vIdx, src.vIdx, sizeof( int  ) * 4 ); 
		memcpy( coef, src.coef, sizeof(double) * 4 );
	}
	inline void setAsPos( const vector<TVector3> &verts, const TTetra &tet, const TVector3 &pos)
	{
		for(int i=0;i<4; ++i) vIdx[i] = tet.v[i];
		t_calcBaryCoordTetra( pos, verts[ tet.v[0] ], verts[ tet.v[1] ], 
								   verts[ tet.v[2] ], verts[ tet.v[3] ], coef);	
	}

	inline void calcInterpolation( const vector<TVector3> &verts, TVector3 &trgtP)const{ 
		trgtP.SetAdditionWithCoef( coef[0], verts[ vIdx[0] ], 
			                       coef[1], verts[ vIdx[1] ],
								   coef[2], verts[ vIdx[2] ], 
								   coef[3], verts[ vIdx[3] ]);
	}
	inline void calcInterpolation( const TVector3 *verts, TVector3 &trgtP)const{ 
		trgtP.SetAdditionWithCoef( coef[0], verts[ vIdx[0] ], 
			                       coef[1], verts[ vIdx[1] ],
								   coef[2], verts[ vIdx[2] ], 
								   coef[3], verts[ vIdx[3] ]);
	}
	inline void calcInterpolation( const vector< double   >  &vals,  double    &val  )const{ 
		val =  coef[0] * vals[ vIdx[0] ] + 
			   coef[1] * vals[ vIdx[1] ] + 
			   coef[2] * vals[ vIdx[2] ] + 
			   coef[3] * vals[ vIdx[3] ] ; 
	}
};




class TTetraModel
{
public:
	static const int FACE_VTX_ORDER[4][3];
	vector<TVector3> m_verts ;
	vector<TTetra    > m_tetras;
	vector<TTetEdge  > m_edges ;

	TTetraModel(void);
	~TTetraModel(void);
	inline void clear(){
		m_verts .clear();
		m_tetras.clear();
		m_edges .clear();
	}

	bool loadfile_EleNode( const char *fname_noExt );
	bool savefile_EleNode( const char *fname_noExt );
	bool loadfile_ttet   ( const char *fname       );
	bool savefile_ttet   ( const char *fname       );

	void print        ();
	void calcNormal   ();
	void calcNeighbor ();
	void calcEdgesFromTetras();
	void checkVtxOrder    ();
	void checkStrangeTetra();

	void initializeSimpleTetraMesh();
	void initializeAsCylinder       (int vNumC, int vNumL, double rad, double length   );
	void initializeAsRotatedCylinder(int vNumC, int vNumL, double rad, double circleRad);
	void initializeAsRotatedCylinderWithFlagellumBase
		                            (int vNumC, int vNumL, double rad, double circleRad, double rootRad);

	double calcVolume() const;
	static double calcVolume(const vector<TVector3> &verts, const vector<TTetra> &tetras);

	///////////////////////////////////////////////////////////////////////////////////////
	//表面surfaceを計算////////////////////////////////////////////////////////////////////
	static void toPolygonModel( const vector< TTetra   > &tetras  ,
		                        const vector< TVector3 > &tetVerts,

									  vector< TVector3 > &Vs,
									  vector< TTriangle> &Ps,
									  vector< int      > &map_tetV2polyV );

	static bool toPolygonModel( const vector< TTetra   > &tetras, const vector< TVector3 > &tetVerts,  TTriangleMesh &surface,  vector< int > &map_tetV2polyV );

	void toPolygonModel( TTriangleMesh &resultMesh) const;
	void toPolygonModel( TTriangleMesh &resultMesh, const vector< TVector3> &newVerts ) const;
	void toPolygonModel( TTriangleMesh &resultMesh,                                      vector< int  >  &map_tetV2polyV) const;

	void toPolygonModel( TTriangleMesh &resultMesh, const vector< TVector3> &newVerts, vector< int  >  &map_tetV2polyV) const;




	void TranslateAllVertices(const TVector3 &t){for(int i=0;i<(int)m_verts.size();++i) m_verts[i]+=t;}
	void ScaleAllVertices    (const double     &d){for(int i=0;i<(int)m_verts.size();++i) m_verts[i]*=d;}

	void subdivision               (vector<byte> &trgtTetFlg , bool doEdgeThresholding = false, double edgeLength= 1.0);
	void subdivision_divideEdge    (const vector<byte> &trgtEdgeFlg, const vector<TVector3> &edgeDivPos );
	void subdivision_notDivBoundary(vector<byte> &trgtTet);
	void sortVerts_By1RingRegionGrow(); //memory効率を上げるためにvertexをソート(1ringを利用してregion growする)
	void sortVerts_ByXvalue         (); //memory効率を上げるためにvertexをソート(1ringを利用してregion growする)



	bool pickAllTetra( const TVector3 &rayP, const TVector3 &rayD, const vector<byte> &ignorTetList, TVector3 &pos, int &tetId, double &minDist) const;
	bool pickAllTetra( const TVector3 &rayP, const TVector3 &rayD, TVector3 &pos, int &tetId) const;
	bool pickOneTetra( const TVector3 &rayP, const TVector3 &rayD, TVector3 &pos, int  tetId);

	int  getEdgeIdx( const int &v0, const int &v1);//注意 too slow implimentation

	void laplacianSmoothing( int iterationN= 10);
void faceFlipping(int maxNum);
private:
	bool do23Flip(int tId);
	bool do32Flip(int tId, const vector<vector<int>> &e_neiborTets);
};



inline void t_computeCircumSphere( const TVector3 &x0, 
	                               const TVector3 &x1, 
							       const TVector3 &x2, 
							       const TVector3 &x3, TVector3 &p, double &r)
{
	TVector3 a, b, c, ab, bc, ca;
	a.SetSubtract( x1,x0);
	b.SetSubtract( x2,x0);
	c.SetSubtract( x3,x0);
	ab.SetCrossProd( a, b );
	bc.SetCrossProd( b, c );
	ca.SetCrossProd( c, a );

	double aa = a*a;
	double bb = b*b;
	double cc = c*c;

	p.SetAdditionWithCoef( aa, bc, bb, ca, cc, ab);
	p /= (2*a*bc);

	p+=x0;


	r = (aa*bc + bb*ca + cc*ab).Length() /(2*fabs(a * bc));
	/*
	fprintf( stderr, "%f %f %f %f\n", (x0-p).Length() - r, 
		                              (x1-p).Length() - r, 
		                              (x2-p).Length() - r, 
		                              (x3-p).Length() - r);	
									  */
}



inline bool t_sphereTest_bIn( const TVector3 &x0,  const TVector3 &x1, 
							  const TVector3 &x2,  const TVector3 &x3, const TVector3 &p)
{
	TVector3 center;
	double r;
	t_computeCircumSphere( x0, x1,x2,x3, center, r);
	return t_distance( center, p) <= r;
}








inline void t_calcNearestTetraBarycenter(const vector<TTetra>     &tetras, 
	                                     const vector<TVector3> &verts,
										 const TVector3         &pos, 
										       TTetBaryCoord      &bary)
{
	//最近傍tetra検索
	int nearestTetIdx = -1;
	double dist       = DBL_MAX;

	TVector3 gCenter;
	int idx = 0;
	for( vector< TTetra >::const_iterator tetI = tetras.begin(); tetI != tetras.end(); ++tetI, ++idx)
	{
		gCenter.SetGravityCenter(verts[ tetI->v[0]], verts[ tetI->v[1]], verts[ tetI->v[2]], verts[ tetI->v[3]]);
		double d = t_distance_sq( gCenter, pos );
		if( d < dist ) {dist = d; nearestTetIdx = idx;}
	}

	//barycentric coordinateにする
	bary.setAsPos( verts, tetras[nearestTetIdx], pos );
}
