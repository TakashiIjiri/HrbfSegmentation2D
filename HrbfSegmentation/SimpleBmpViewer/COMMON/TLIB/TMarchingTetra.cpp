#include "StdAfx.h"
#include "TMarchingTetra.h"

#include <map>


static const int pairs[][2] = {{0, 3}   , {0, 2}   , {1, 2}   , {1, 3}   };
static const int mapping_pattern[][4] = 
{
	{-1, -1, -1, -1},	//0
	{0, 1, 2, 3},	//1 +
	{1, 2, 0, 3},	//2 +
	{0, 1, 2, 3},	//3
	{2, 0, 1, 3},	//4 +
	{0, 2, 3, 1},	//5
	{1, 2, 0, 3},	//6
	{3, 2, 1, 0},	//7
	{3, 2, 1, 0},	//8 +
	{0, 3, 1, 2},	//9
	{1, 3, 2, 0},	//10
	{2, 0, 1, 3},	//11
	{2, 3, 0, 1},	//12
	{1, 2, 0, 3},	//13
	{0, 1, 2, 3},	//14
	{0, 1, 2, 3}	//15
};


//calc index key
//key is a  unique key for each surface vertices
// if the surface vertex corresponds to tetra edge idx0-idx1 (idx0<idx1) : key is idx0 * tetra vertex size + idx1
// if the surface vertex corresponds to a tetra vertex idx0              : key is idx0 * tetra vertex size + idx0
inline static int t_calcVtxKeyOnEdge( const int idx0, const int idx1, const int vtxSize)
{
	//int index0 = min(idx0, idx1);
	//int index1 = max(idx0, idx1);
	//wKey[j] = vtxSize * index0 + index1;
	return vtxSize * min(idx0, idx1) + max(idx0, idx1);
}

inline static double t_calc_innerDivCoef( double v0, double v1, double thresh)
{
	if( v0<v1 ) return (thresh-v0)/(v1-v0);
	else        return (thresh-v1)/(v0-v1);
}


//t_addPoly    : This method adds vertices (if necessary ) and a polygon to new surface (Vs//Ps);
//Vs, Ps       : surface vertices and polygons
//existVtxIdMap: a data structure to check where a surface vertex is already pushed back or not 
//               a surface vertex that corresponds to tetra edge (v0-v1) or to tetra vertex (v0)
//                   was already pushed back --> existVtxIdMap[v0].find( v1 ) return true and existVtxIdMap[v0].find( v1 )->second is a vertex ID
inline static void t_addPoly(  vector< TVector3       > &Vs, 
	                           vector< TTriangle        > &Ps,  
							   vector< map<int, int>    > &existVtxIdMap, 
					           vector< TMapSurfV_TetV> &sV_info, 
							   const TVector3     &vPos0, const TVector3     &vPos1, const TVector3     &vPos2, 
							   const TMapSurfV_TetV &vKey0, const TMapSurfV_TetV &vKey1, const TMapSurfV_TetV &vKey2)
{
	const TVector3*     vPos_p[3] = {&vPos0, &vPos1, &vPos2};
	const TMapSurfV_TetV* vKey_p[3] = {&vKey0, &vKey1, &vKey2};
		
	int vID[3];//polygon vertex ids // search existing vertex (hit--> use the found id / not exist push back new vertex)
	for (int i = 0; i < 3; ++i) 
	{
		const int &idx0 = vKey_p[i]->idx0;
		const int &idx1 = vKey_p[i]->idx1;
		map<int, int>::iterator it = existVtxIdMap[ idx0 ].find( idx1 );

		if ( it != existVtxIdMap[idx0].end() ) {vID[i] = it->second; continue; }
		
		vID[i] = (int)Vs.size();
		Vs    .push_back( *vPos_p[i] );
		existVtxIdMap[idx0].insert( make_pair( idx1, vID[i]) );
		existVtxIdMap[idx1].insert( make_pair( idx0, vID[i]) );
				
		sV_info.push_back( *vKey_p[i] );
	}

	Ps.push_back( TTriangle(vID[0], vID[1], vID[2]) );
}


//////t_calcMarchingTetra///////////////////////////////////////////////////
//Marching tetraÇ…ÇÊÇËsurface meshÇê∂ê¨Ç∑ÇÈ
//tetras/tetVertsÇ…tetraMeshèÓïÒÇ™ì¸ÇÈ
void t_calcMarchingTetra(const vector<TTetra    > &tetras,
						 const vector<TVector3> &tetVerts,  
						 const vector<double    > &vtxValue, 
						 const double               threshold, 
						 TTriangleMesh            &trgt,
						 vector<TMapSurfV_TetV>   &surfVinfo)
{	
	TVector3     vPos[4], wPos[4];//vPos, vKey is vertex of original tetra mesh
	TMapSurfV_TetV vKey[4], wKey[4];//wPos, wKey is interpolated vertex on the edge of tetra mesh
	double         val [4]         ;//tetra vertex vallue

	surfVinfo.clear();
	vector< TVector3 > Vs;
	vector< TTriangle  > Ps;
	vector< map<int, int> >  existVtxIdMap( tetVerts.size() );
	
	const int vtxSize = (int)tetVerts.size();
	const int tetSize = (int)tetras  .size();

	for (int i = 0; i < tetSize; ++i) 
	{
		const TTetra& tet = tetras[i];
		int code = 0;//0x0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111ÇÃ16í ÇË
		for (int j = 0; j < 4; ++j) {
			val[j] = vtxValue[ tet.v[j] ];
			if (val[j] < threshold) code += (1 << j);
		}

		if (code == 0) continue;

		const int* mapping = mapping_pattern[code];
		
		for (int j = 0; j < 4; ++j) {
			vPos[j] = tetVerts[ tet.v[ mapping[j] ]];
			vKey[j].Set(        tet.v[ mapping[j] ]);
		}
		if (code == 15)// code = 0x1111 all vertex is inside
		{
			if(tet.n[0] == -1) t_addPoly(Vs,Ps,existVtxIdMap, surfVinfo, vPos[1], vPos[2], vPos[3], vKey[1], vKey[2], vKey[3]);
			if(tet.n[1] == -1) t_addPoly(Vs,Ps,existVtxIdMap, surfVinfo, vPos[3], vPos[2], vPos[0], vKey[3], vKey[2], vKey[0]);
			if(tet.n[2] == -1) t_addPoly(Vs,Ps,existVtxIdMap, surfVinfo, vPos[1], vPos[3], vPos[0], vKey[1], vKey[3], vKey[0]);
			if(tet.n[3] == -1) t_addPoly(Vs,Ps,existVtxIdMap, surfVinfo, vPos[2], vPos[1], vPos[0], vKey[2], vKey[1], vKey[0]);
		}
		else if( code == 1 || code == 2 || code == 4 || code == 8) //code = 0x0001 0010 0100 1000 one vertex is inside
		{
			for (int j = 0; j < 3; ++j) 
			{
				int idx0 = mapping[  0  ];
				int idx1 = mapping[j + 1];
				double u = (threshold - val[ idx0 ]) / (val[ idx1 ] - val[ idx0 ]);
				wPos[j].SetAdditionWithCoef( 1-u, vPos[ 0 ], u, vPos[j+1] );
				wKey[j].Set( tet.v[idx0], tet.v[idx1], u);
			}

			t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, wPos[0], wPos[1], wPos[2], wKey[0], wKey[1], wKey[2]);

			if (tet.n[mapping[1]] == -1) t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, wPos[1], vPos[0], wPos[2],  wKey[1], vKey[0], wKey[2]);
			if (tet.n[mapping[2]] == -1) t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, wPos[2], vPos[0], wPos[0],  wKey[2], vKey[0], wKey[0]);
			if (tet.n[mapping[3]] == -1) t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, wPos[0], vPos[0], wPos[1],  wKey[0], vKey[0], wKey[1]);
		}
		else if( code == 3 || code == 5 || code == 6 || code == 9 || code == 10 || code == 12)
		{
			for (int j = 0; j < 4; ++j) 
			{
				const int& i0 = pairs[j][0], idx0 = mapping[ i0 ];
				const int& i1 = pairs[j][1], idx1 = mapping[ i1 ];
				double u = (threshold - val[ idx0 ]) / (val[ idx1 ] - val[ idx0 ]);
				wPos[j].SetAdditionWithCoef( 1-u, vPos[ i0 ], u, vPos[i1] );
				wKey[j].Set( tet.v[idx0], tet.v[idx1], u);
			}
			t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo,wPos[0], wPos[2], wPos[1], wKey[0], wKey[2], wKey[1]);
			t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo,wPos[0], wPos[3], wPos[2], wKey[0], wKey[3], wKey[2]);
			
			if (tet.n[mapping[0]] == -1)  t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, vPos[1], wPos[2], wPos[3],  vKey[1], wKey[2], wKey[3]);
			if (tet.n[mapping[1]] == -1)  t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, vPos[0], wPos[0], wPos[1],  vKey[0], wKey[0], wKey[1]);
			if (tet.n[mapping[2]] == -1) {t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, vPos[1], wPos[3], wPos[0],  vKey[1], wKey[3], wKey[0]);
				                          t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, vPos[1], wPos[0], vPos[0],  vKey[1], wKey[0], vKey[0]);}
			if (tet.n[mapping[3]] == -1) {t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, vPos[0], wPos[1], wPos[2],  vKey[0], wKey[1], wKey[2]);
				                          t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo, vPos[0], wPos[2], vPos[1],  vKey[0], wKey[2], vKey[1]);}
		}
		else if( code == 7 || code == 11 || code == 13 || code == 14)
		{
			for (int j = 0; j < 3; ++j)
			{
				int idx0 = mapping[  0  ];
				int idx1 = mapping[j + 1];
				double u = (threshold - val[ idx0 ]) / (val[ idx1 ] - val[ idx0 ]);
				wPos[j].SetAdditionWithCoef( 1-u, vPos[ 0 ], u, vPos[j+1] );
				wKey[j].Set( tet.v[idx0], tet.v[idx1], u);
			}
			t_addPoly(Vs,Ps,existVtxIdMap,surfVinfo,wPos[2], wPos[1], wPos[0], wKey[2], wKey[1], wKey[0]);
		
			if (tet.n[mapping[0]] == -1) t_addPoly( Vs,Ps,existVtxIdMap, surfVinfo, vPos[1], vPos[2], vPos[3], vKey[1], vKey[2], vKey[3]);
			if( tet.n[mapping[1]] == -1){t_addPoly( Vs,Ps,existVtxIdMap, surfVinfo, vPos[2], wPos[1], wPos[2], vKey[2], wKey[1], wKey[2]);
										 t_addPoly( Vs,Ps,existVtxIdMap, surfVinfo, vPos[2], wPos[2], vPos[3], vKey[2], wKey[2], vKey[3]);}
			if( tet.n[mapping[2]] == -1){t_addPoly( Vs,Ps,existVtxIdMap, surfVinfo, vPos[3], wPos[2], wPos[0], vKey[3], wKey[2], wKey[0]);
										 t_addPoly( Vs,Ps,existVtxIdMap, surfVinfo, vPos[3], wPos[0], vPos[1], vKey[3], wKey[0], vKey[1]);}
			if( tet.n[mapping[3]] == -1){t_addPoly( Vs,Ps,existVtxIdMap, surfVinfo, vPos[1], wPos[0], wPos[1], vKey[1], wKey[0], wKey[1]);
		    	                         t_addPoly( Vs,Ps,existVtxIdMap, surfVinfo, vPos[1], wPos[1], vPos[2], vKey[1], wKey[1], vKey[2]);}
		}
	}
	trgt.initFromVertsPolys( Vs, Ps ) ; 
}











