#include "StdAfx.h"

#include "T_FlowFieldByRBF3D.h"


T_FlowFieldByRBF3D::T_FlowFieldByRBF3D(void)
{
	m_X_lambda_abc = m_Y_lambda_abc = m_Z_lambda_abc = 0;
	m_N = 0;//numOfConst
}



T_FlowFieldByRBF3D::~T_FlowFieldByRBF3D(void)
{
	if( m_X_lambda_abc != 0) delete[] m_X_lambda_abc;
	if( m_Y_lambda_abc != 0) delete[] m_Y_lambda_abc;
	if( m_Z_lambda_abc != 0) delete[] m_Z_lambda_abc;
}



void T_FlowFieldByRBF3D::addConstraints   ( const vector< pair<TVector3, TVector3> > &constPosAndDir )
{
	m_N += (int) constPosAndDir.size();

	for( int i = 0; i < (int) constPosAndDir.size(); ++i)
	{
		m_constPos.push_back( constPosAndDir[i].first  ); 
		m_constDir.push_back( constPosAndDir[i].second ); 
	}

	//update field//////////////////////////////////////////////////////////
	if( m_X_lambda_abc != 0) delete[] m_X_lambda_abc;  	m_X_lambda_abc = new double[m_N + 1];
	if( m_Y_lambda_abc != 0) delete[] m_Y_lambda_abc;   m_Y_lambda_abc = new double[m_N + 1];
	if( m_Z_lambda_abc != 0) delete[] m_Z_lambda_abc;	m_Z_lambda_abc = new double[m_N + 1];
	memset( m_X_lambda_abc, 0 , sizeof( double ) * ( m_N + 1 ) );
	memset( m_Y_lambda_abc, 0 , sizeof( double ) * ( m_N + 1 ) );
	memset( m_Z_lambda_abc, 0 , sizeof( double ) * ( m_N + 1 ) );

	updateRBFWeight();
}


void T_FlowFieldByRBF3D::addConstraints( const vector<TVector3> &constPos, const vector< TVector3 > &constDir)
{
	if( constPos.size() != constDir.size()) {fprintf(stderr, "error1232213213"); return;}
	m_N += (int) constPos.size();

	for( int i = 0; i < (int) constPos.size(); ++i)
	{
		m_constPos.push_back( constPos[i] ); 
		m_constDir.push_back( constDir[i] ); 
	}

	//update field//////////////////////////////////////////////////////////
	if( m_X_lambda_abc != 0) delete[] m_X_lambda_abc;  	m_X_lambda_abc = new double[m_N + 1];
	if( m_Y_lambda_abc != 0) delete[] m_Y_lambda_abc;   m_Y_lambda_abc = new double[m_N + 1];
	if( m_Z_lambda_abc != 0) delete[] m_Z_lambda_abc;	m_Z_lambda_abc = new double[m_N + 1];
	memset( m_X_lambda_abc, 0 , sizeof( double ) * ( m_N + 1 ) );
	memset( m_Y_lambda_abc, 0 , sizeof( double ) * ( m_N + 1 ) );
	memset( m_Z_lambda_abc, 0 , sizeof( double ) * ( m_N + 1 ) );

	updateRBFWeight();
}




void T_FlowFieldByRBF3D::modifyConstraints( const vector<TVector3> &constPos, const vector< TVector3 > &constDir)
{
	if( m_constPos.size() != constPos.size() || m_constDir.size() != constDir.size() ) {fprintf(stderr, "error134214321"); return;}
	m_constPos = constPos;
	m_constDir = constDir;

	updateRBFWeight();
}







// rbf( x ) = sigrma ( ai * f(x-ci) ) + c //–{“–‚Í ax + by + cz + d 
void T_FlowFieldByRBF3D::updateRBFWeight()
{
	double **A;
	A = new double*[m_N + 1];
	for( int i = 0; i < m_N + 1; ++i) A[i] = new double[m_N +1];

	double *m_X_b = new double[m_N + 1];
	double *m_Y_b = new double[m_N + 1];
	double *m_Z_b = new double[m_N + 1];


	//set up A////////////////////////////////////////////////////////////
	for( int i = 0; i < m_N; ++i)
	for( int j = i; j < m_N; ++j)
		A[i][j] = A[j][i] = radFunc( t_distance_sq( m_constPos[i], m_constPos[j] ) );
	
	for( int i = 0; i < m_N; ++i)
		A[ i ][ m_N   ] = A[ m_N   ][ i ] = 1;

	A[m_N][m_N]=0;


	//set up constraints to b/////////////////////////////////////////////
	for( int i = 0; i < m_N; ++i){
		m_X_b[i] = m_constDir[i].data[0];
		m_Y_b[i] = m_constDir[i].data[1];
		m_Z_b[i] = m_constDir[i].data[2];
	}
	m_X_b[ m_N ] = 0;
	m_Y_b[ m_N ] = 0;
	m_Z_b[ m_N ] = 0;
	
	//solve for X Y Z Dir///////////////////////////////////////////////////
	t_solveLinearEquationByCG( m_N + 1, A, m_X_b, m_X_lambda_abc);
	t_solveLinearEquationByCG( m_N + 1, A, m_Y_b, m_Y_lambda_abc);
	t_solveLinearEquationByCG( m_N + 1, A, m_Z_b, m_Z_lambda_abc);

	delete[] m_X_b;	
	delete[] m_Y_b;
	delete[] m_Z_b;
	for( int i = 0; i < m_N + 1; ++i) delete[] A[i];
	delete[] A;
}



void T_FlowFieldByRBF3D::getValue( const TVector3 &v, TVector3& orientation )
{
	if( m_N == 0) return;
	orientation.Set( 0,0,0);

	for( int i = 0; i < m_N; ++i)
	{
		double f = radFunc( t_distance_sq( v, m_constPos[i] ) );
		orientation.data[0] += m_X_lambda_abc[i] * f;
		orientation.data[1] += m_Y_lambda_abc[i] * f;
		orientation.data[2] += m_Z_lambda_abc[i] * f;
	}
	orientation.data[0] += m_X_lambda_abc[m_N];
	orientation.data[1] += m_Y_lambda_abc[m_N];
	orientation.data[2] += m_Z_lambda_abc[m_N];
	orientation.Normalize_Self();
}


void T_FlowFieldByRBF3D::deleteAllConstraints()
{
	m_N = 0;
	m_constDir.clear();
	m_constPos.clear();
	if( m_X_lambda_abc != 0) delete[] m_X_lambda_abc; m_X_lambda_abc = 0;
	if( m_Y_lambda_abc != 0) delete[] m_Y_lambda_abc; m_Y_lambda_abc = 0;
	if( m_Z_lambda_abc != 0) delete[] m_Z_lambda_abc; m_Z_lambda_abc = 0;
}



void T_FlowFieldByRBF3D::renderVis()
{
	glPushMatrix();
	glTranslated( 0,0,0.1);

	glDisable( GL_LIGHTING );
	glLineWidth(  8   );
	glColor3d(  0,1,1 );
	glBegin( GL_LINES );
	for( int i = 0 ; i < m_N; ++i)
	{
		glVertex3dv( m_constPos[i].data );
		glVertex3dv( (m_constPos[i] + m_constDir[i]).data  );
	}
	glEnd();
	glPopMatrix();
}


