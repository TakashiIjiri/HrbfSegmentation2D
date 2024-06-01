#pragma once

#include "TTriangleMesh.h"
#include "TTetraModel.h"
#include <map>

//TMapSurfV_TetV//////////////////////////////////////////
//mapping information from surface vertex to tetra vertex
//idx0 and idx1 are corresponding tetra vertex ids (if surface vertex exist on a tetra vertex idx0 == idx1 and isOnEdge=false)
//surface vertexからtetra vertexへのマッピング　
//tetra のエッジに対応する場合は 内分する係数も保存
class TMapSurfV_TetV
{
public:
	bool isOnEdge ;
	int idx0, idx1;
	double u;
	TMapSurfV_TetV( ) {
		isOnEdge = false;
		idx0 = 0;
		idx1 = 0;
		u    = 0;
	}
	TMapSurfV_TetV( bool _isOnEdge, int _idx0, int _idx1 = -1, double _u = 0){
		isOnEdge = _isOnEdge;
		idx0     = _idx0;
		idx1     = _idx1;
		u        = _u;
	}
	void Set( int _idx0){
		isOnEdge = false;
		idx0 = idx1 = _idx0;
		u=0;
	}
	void Set( int _idx0, int _idx1, double _u)
	{
		isOnEdge = true;
		idx0 = _idx0;
		idx1 = _idx1;
		u = _u;
	}
};

////////////////////////////////////////////////////////////////////////////////////
//tetras, tetVerts : tetrahedral mesh
//vtxValue         : double value at all vertices
//threshold        : threshold for marching tetra
//trgt             : resulting mesh
//surfVinfo        : surface vertex info (mapping from surface vertices to tetra vertices (2 vertex ids and coeficients)
void t_calcMarchingTetra(const vector<TTetra    > &tetras,
						 const vector<TVector3> &tetVerts,  
						 const vector<double    > &vtxValue, 
						 const double threshold, 
						 TTriangleMesh &trgt,
						 vector<TMapSurfV_TetV> &surfVinfo);