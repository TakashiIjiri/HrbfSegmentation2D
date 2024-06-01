#include "StdAfx.h"
#include "TRBF3D.h"
#include "TOGL.h"


TRBF3D::~TRBF3D(void){
	if( m_lambda_abc != 0) delete[] m_lambda_abc;
}

TRBF3D::TRBF3D(void)
{
	m_lambda_abc = 0;
	m_N          = 0;//numOfConst
}

void TRBF3D::addConstraints( const vector<TVector3> &constPos, const vector< double > &constVal)
{
	if( constPos.size() != constVal.size()) {fprintf(stderr, "error1232213213"); return;}
	
	m_N += (int) constPos.size();
	for( int i = 0; i < (int) constPos.size(); ++i){
		m_constPos.push_back( constPos[i] ); 
		m_constVal.push_back( constVal[i] ); 
	}

	//update field//////////////////////////////////////////////////////////
	if( m_lambda_abc != 0) delete[] m_lambda_abc;  	m_lambda_abc = new double[m_N + 1];
	memset( m_lambda_abc, 0 , sizeof( double ) * ( m_N + 1 ) );

	updateRBFWeight();
}



void TRBF3D::modifyConstraints( const vector<TVector3> &constPos, const vector< double > &constVal)
{
	if( m_constPos.size() != constPos.size() || m_constVal.size() != constVal.size() ) {fprintf(stderr, "error134214321"); return;}
	m_constPos = constPos;
	m_constVal = constVal;

	updateRBFWeight();
}



void TRBF3D::deleteAllConstraints()
{
	m_N = 0;
	m_constVal.clear();
	m_constPos.clear();
	if( m_lambda_abc != 0) delete[] m_lambda_abc; m_lambda_abc = 0;
}




// rbf( x ) = sigrma ( ai * f(x-ci) ) + c //–{“–‚Í ax + by + cz + d 
void TRBF3D::updateRBFWeight()
{

	const double coef = 1; //Ax = b‚ðˆÀ’è‚É‰ð‚­‚½‚ß, cA x  = c b ‚Æ‚·‚é
	
	double **A;
	A = new double*[m_N + 1];
	for( int i = 0; i < m_N + 1; ++i) A[i] = new double[m_N +1];

	double *m_b = new double[m_N + 1];

	//set up A////////////////////////////////////////////////////////////
	for( int i = 0; i < m_N; ++i)
	for( int j = i; j < m_N; ++j)
		A[i][j] = A[j][i] = coef*radFunc( t_distance_sq( m_constPos[i], m_constPos[j] ) );
	
	for( int i = 0; i < m_N; ++i) 
		A[i][m_N] = A[m_N][i] = coef*1;
	
	A[m_N][m_N]=0;

	//set up constraints to b/////////////////////////////////////////////
	for( int i = 0; i < m_N; ++i) m_b[i] = coef*m_constVal[i];
	m_b[ m_N ] = 0;

	//solve for X Y Z Dir///////////////////////////////////////////////////
	t_solveLinearEquationByCG( m_N + 1, A, m_b, m_lambda_abc);

	delete[] m_b;	
	for( int i = 0; i < m_N + 1; ++i) delete[] A[i];
	delete[] A;
}



void TRBF3D::renderVis()
{
	glPushMatrix();
	glTranslated( 0,0,0.1);

	glDisable( GL_LIGHTING );
	glPointSize( 5 );

	glBegin( GL_POINTS );
	for( int i = 0 ; i < m_N; ++i){
		if( m_constVal[i] < 0 ) glColor3f( 0, (float)m_constVal[i], 0);
		else                    glColor3f( (float)m_constVal[i], 0, 0);
		glVertex3dv( m_constPos[i].data );
	}
	glEnd();

	glPopMatrix();
}

