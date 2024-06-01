#pragma once

#include "tmath.h"
#include "TOGL.h"

class T_FlowFieldByRBF3D
{
	double *m_X_lambda_abc, *m_Y_lambda_abc, *m_Z_lambda_abc;
	int m_N;//num of constraints
	vector<TVector3> m_constPos;
	vector<TVector3> m_constDir;

public:
	T_FlowFieldByRBF3D(void);
	~T_FlowFieldByRBF3D(void);

	void addConstraints   ( const vector< pair<TVector3, TVector3>> &constPosAndDir );
	void addConstraints   ( const vector<TVector3> &constPos, const vector< TVector3 > &contDir);
	void modifyConstraints( const vector<TVector3> &constPos, const vector< TVector3 > &contDir);
	void deleteAllConstraints();

	void getValue( const TVector3 &v, TVector3& orientation );
	void renderVis( );
private:
	void updateRBFWeight();

	inline double radFunc( double distSquare)
	{
		return sqrt( distSquare );
		//r*r*r
		//return distSquare  * sqrt( distSquare );
		//r*r*logr
		//if( distSquare == 0) return 0;
		//return distSquare * log( sqrt(distSquare) );
	}
};
