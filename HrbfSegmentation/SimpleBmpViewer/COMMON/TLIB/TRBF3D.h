#pragma once

#include "tmath.h"
#include <vector>

class TRBF3D
{
	int m_N;//num of constraints
	double *m_lambda_abc;
	vector<TVector3> m_constPos;
	vector<double    > m_constVal;


public:
	TRBF3D(void);
	~TRBF3D(void);

	void   addConstraints   ( const vector<TVector3> &constPos, const vector< double > &contVal);
	void   deleteAllConstraints();
	void   renderVis();
	void   modifyConstraints( const vector<TVector3> &constPos, const vector< double > &constVal);

private:
	void updateRBFWeight();
	inline double radFunc( const double &distSquare)
	{
		return sqrt( distSquare );
		//r*r*r
		//return distSquare  * sqrt( distSquare );
		//r*r*logr
		//if( distSquare == 0) return 0;
		//return distSquare * log( sqrt(distSquare) );
	}
public:		
	inline double getValue( const TVector3 &pos)
	{
		if( m_N == 0) return 0;

		double val = 0;
		for( int i = 0; i < m_N; ++i)
			val += m_lambda_abc[i] * radFunc( t_distance_sq( pos, m_constPos[i] ) );
		val += m_lambda_abc[m_N];
		return val;
	}


};
