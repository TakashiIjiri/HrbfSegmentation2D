#include "stdafx.h"

#include "Tstroke3D.h"
#include "TMatrix.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


//////////////////////////////////////////////////////////////////////////////////////////
//Draw stroke/////////////////////////////////////////////////////////////////////////////
#ifdef USE_TOGL
void Tstroke3D::drawLineStrip(float width) const
{
	glLineWidth( width );
	glBegin( GL_LINE_STRIP );
	for( int i= 0; i<(int)m_points.size(); ++i) glVertex3dv( m_points[i].data );
	if( m_isClosed && !m_points.empty()) glVertex3dv( m_points.front().data );
	glEnd();
}
#endif 



void t_stroke_Smooth_BSpline3Interpolation(const vector<TVector3> &CPs, vector<TVector3> &target, bool isClosed, int NumOfPoint)
{
	//p(t) = Σ_{j=[-1, N-1+1} N3(t-j) Qj   0 <= t <= N-1 
	//p(i) = CPs[i]を満たすようにQiを決定
	const double c  = 1/6.0;
	const double c4 = 4/6.0;
	const int N = (int)CPs.size();

	if     ( N == 0 ) {target.resize( NumOfPoint        ); return;}
	else if( N == 1 ) {target.resize( NumOfPoint, CPs[0]); return;}
	else if( N == 2 ) {t_stroke_devideEquals( NumOfPoint-1, CPs, target); return;} 

	TSparseMatrix M( N, N);

	for( int i=1; i< N-1; ++i){ 
		M.push_back( i, i-1, c);
		M.push_back( i, i  , c4);
		M.push_back( i, i+1, c);
	}
	if( isClosed ){ M.push_back( 0  , 0  ,  c4); M.push_back( 0  , 1  ,   c); M.push_back(   0, N-1, c);
		            M.push_back( N-1,  0 ,   c); M.push_back( N-1, N-2,  c ); M.push_back( N-1, N-1,  c4); }
	else{           M.push_back( 0  , 0  , 5*c); M.push_back( 0  , 1  ,   c);
		            M.push_back( N-1, N-2,   c); M.push_back( N-1, N-1, 5*c);}
	double *Qx = new double[N], *Px = new double[N];
	double *Qy = new double[N], *Py = new double[N];
	double *Qz = new double[N], *Pz = new double[N];
	for( int i=0; i<N; ++i){Px[i] = CPs[i][0];
							Py[i] = CPs[i][1];
							Pz[i] = CPs[i][2];}

	M.createFieldForSolveLinearSystem( isClosed ? (N * 3 ) : (N * 3 - 2) );
	
	M.umfPack_Prepare();
	M.umfPack_Solve_forPrecomp( Px, Qx);
	M.umfPack_Solve_forPrecomp( Py, Qy);
	M.umfPack_Solve_forPrecomp( Pz, Qz);
	M.umfPack_Release();


	// 0 <= t <= N-1
	double stepD = (N-1) / (double)(NumOfPoint -1); //３点なら２等分
	TVector3 pos;
	target.resize( NumOfPoint );
	if( isClosed ){
		for( int ii=0; ii<NumOfPoint; ++ii)
		{
			double t = ii * stepD;
			pos.Set( 0,0,0 );
			for( int j = -1; j < N + 1; ++j)
			{
				if( fabs( t - j ) > 2 ) continue; //t_bSpline3は [-2, 2]の外ではゼロ
				int index = ( j < 0  ) ? N-1   : ( j > N-1) ? 0 : j ;
				double g = t_bSpline3(t - j);
				pos[0] += g * Qx[index];
				pos[1] += g * Qy[index];
				pos[2] += g * Qz[index];
			}
			target[ii] = pos;
		}
	}else{
		for( int ii=0; ii<NumOfPoint; ++ii)
		{
			double t = ii * stepD;
			pos.Set( 0,0,0 );
			for( int j = -1; j < N + 1; ++j)
			{
				if( fabs( t - j ) > 2 ) continue; //t_bSpline3は [-2, 2]の外ではゼロ
				int index = ( j < 0  ) ? 0   : ( j > N-1) ? N-1 : j ;
				double g = t_bSpline3(t - j);
				pos[0] += g * Qx[index];
				pos[1] += g * Qy[index];
				pos[2] += g * Qz[index];
			}
			target[ii] = pos;
		}
	}

	delete[] Qx; delete[] Qy; delete[] Qz;
	delete[] Px; delete[] Py; delete[] Pz;
}





/////////////////////////////////////////////////////////////////////////////////////////
double	Tstroke3D::Length() const{
	if( size() < 2 ) return 0;
	double d = 0;
	for( int i = 1; i < (int) m_points.size(); ++i ) d += t_distance( m_points[i], m_points[i-1] );
	if( m_isClosed ) d += t_distance( m_points.front(), m_points.back() );
	return d;
}

double	Tstroke3D::Length(int sIndex,int eIndex) const
{
	if( sIndex < 0                      ) sIndex = 0;
	if( eIndex >= (int) m_points.size() ) eIndex = (int)m_points.size() -1;
	
	double d = 0;
	for( int i = sIndex + 1; i <= eIndex; ++i ) 
		d += t_distance( m_points[i], m_points[i-1] );
	return d;
}



void Tstroke3D::devideEquals(int n)
{
	vector<TVector3> result;
	t_stroke_devideEquals( n, m_points, result );
	m_points = result;
}


void Tstroke3D::Smoothing(void)
{
	vector<TVector3> result(m_points.size());
	result[  0            ] = m_points.front();
	for (int i = 1; i < (int)m_points.size() - 1; i++) result[i].SetAdditionWithCoef(	0.25, m_points[i-1],
																						0.5 , m_points[ i ], 
																						0.25, m_points[i+1]);
	result[result.size()-1] = m_points.back ();
	m_points = result;
}


//Smoothing stroke using B-Spline algorithm
//B-Splineの制御点として、現在のストロークを利用し、スムージングをかける
void Tstroke3D::SmoothingUsingBSpline()
{
	vector<TVector3> result;
	t_stroke_Smooth_BSpline3( m_points, result,(int)m_points.size());
	m_points = result;
}


double Tstroke3D::distanceFromEyeRayToPoints(const TVector3 &eyeP, 
											 const TVector3 &rayDir, 
											 int        &nearestPointIdx) const
{
	if( m_points.size() == 0 ) return DBL_MAX;

	double minDist  = DBL_MAX;
	nearestPointIdx = -1;

	TVector3 ndir = rayDir;
	ndir.Normalize_Self();

	for( int i = 0; i < (int) m_points.size(); ++i )
	{
		double d = t_distPointToLine( m_points[i], eyeP, ndir );
		if( d < minDist ){ minDist = d; nearestPointIdx  = i; }
	}
	return minDist;
}

double Tstroke3D::distanceFromEyeRayToCurve (const TVector3 &eyeP, 
	                                         const TVector3 &rayDir,
											 int  &nearestPointIdx) const
{	
	if( m_points.size() <= 1 ) return 0;
	double minDist = distanceFromEyeRayToPoints( eyeP, rayDir, nearestPointIdx);

	for( int i = 1; i < (int) m_points.size(); ++i )
	{
		bool isOnSegm = false;
		double d = t_distLineToLineSegment_sq(eyeP, eyeP + rayDir, m_points[i], m_points[i-1], isOnSegm ); 
		if( isOnSegm )
			minDist = min( sqrt( d ), minDist);
		
	}

	if( isClosed() ){
		bool isOnSegm = false;
		double d = t_distLineToLineSegment_sq(eyeP, eyeP + rayDir, m_points.front(), m_points.back(), isOnSegm ); 
		if( isOnSegm ){
			minDist = min( sqrt( d ), minDist);
		}
	}
	return minDist;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Tstroke3D_tube///////////////////////////////////////////////////////////////////////////////////////////////
Tstroke3D_tube::Tstroke3D_tube(int tubeRes,double radius,bool isClosed)
: Tstroke3D(isClosed)
{
	m_tubeRes = tubeRes;
	m_radius  = radius ;
}

Tstroke3D_tube::Tstroke3D_tube(const vector<TVector3> &ps ,int tubeRes,double radius,bool isClosed)
: Tstroke3D(ps, isClosed)
{
	m_tubeRes = tubeRes;
	m_radius  = radius ;
	createTubeFromAxis();
}

Tstroke3D_tube::Tstroke3D_tube(const Tstroke3D_tube &src)
: Tstroke3D( src.m_points, src.isClosed() )
{
	m_tubeRes = src.m_tubeRes;
	m_radius  = src.m_radius ;
	createTubeFromAxis();
}



#ifdef USE_TOGL
void Tstroke3D_tube::drawTube() const 
{
	glBegin( GL_TRIANGLE_STRIP );
	for( int i = 0; i < (int)m_tubeVerts   .size() -1; ++i)
	for( int j = 0; j < (int)m_tubeVerts[i].size()   ; ++j){
		glNormal3dv( m_tubeNorms[ i ][j].data ); glVertex3dv( m_tubeVerts[ i ][j].data );
		glNormal3dv( m_tubeNorms[i+1][j].data ); glVertex3dv( m_tubeVerts[i+1][j].data );
	}
	glEnd();
}
#endif
	

//Create generalized Cylinder from the axis
//both for initial creation and reflesh the shape
//以下の実装ではglobalなx軸方向を固定して、localな座標系を考慮する
#if 0

void Tstroke3D_tube::createTubeFromAxis()
{
	if( m_points.size() == 0 ) return;

	//prepare fields
	if( m_tubeVerts        .size() == 0               || 
		m_tubeVerts        .size() != m_points.size() ||
		m_tubeVerts.front().size() != m_tubeRes       )
	{
		m_tubeVerts.resize( m_points.size(), vector<TVector3>( m_tubeRes ) );
		m_tubeNorms.resize( m_points.size(), vector<TVector3>( m_tubeRes ) );
	}
	
	TVector3 gCenter;
	TVector3 globalN, n;//fitする平面の法線(近似)
	for( int i = 0; i < (int) m_points.size(); ++i) gCenter += m_points[i];
	for( int i = 1; i < (int) m_points.size(); ++i)
	{
		TI_Set_V1subtV2_outmult_V3subtV4( m_points[i], gCenter, m_points[i-1], gCenter, n);
		if( n*globalN <0 ) globalN -= n;
		else               globalN += n;
	}
	globalN.Normalize_Self();



	TMatrix9  rotM;
	TVector3 localX, localY, v;

	for( int i = 0; i < (int) m_points.size(); ++i)
	{
		localY.Set(0,0,0);
		if( i != 0                        ) localY.AddSubtract( m_points[ i ], m_points[i-1] );
		if( i != (int)m_points.size() - 1 ) localY.AddSubtract( m_points[i+1], m_points[ i ] );
		localX = globalN;
		
		if( !TI_getRotationMatFromLocalXYdir( localX, localY, rotM ))
		{
		}
		
		for( int j = 0; j < m_tubeRes; ++j)
		{
			v.Set( m_radius * cos( 2 * M_PI / (m_tubeRes -1) * j ), 0,
				   m_radius * sin( 2 * M_PI / (m_tubeRes -1) * j ));
			
			TI_MatMultVec( v, rotM );
			m_tubeVerts[i][j].SetAddition( m_points[i], v);
			m_tubeNorms[i][j] = v;
			m_tubeNorms[i][j].Normalize_Self();
		}
	}
}

#else


/*
closed strokeの場合は, 縦のtriangle stripの列をひとつ多くする
*/

void Tstroke3D_tube::createTubeFromAxis()
{
	if( m_points.size() < 2 ) return;
	
	const int vRes = isClosed() ? (int)m_points.size()+1 : (int)m_points.size();//vertical resolution of surface

	//prepare fields
	if( m_tubeVerts.size() != vRes  || (m_tubeVerts.size() != 0 && m_tubeVerts.front().size() != m_tubeRes) )
	{
		m_tubeVerts.clear();
		m_tubeNorms.clear();
		m_tubeVerts.resize( vRes, vector<TVector3>( m_tubeRes ) );
		m_tubeNorms.resize( vRes, vector<TVector3>( m_tubeRes ) );
	}
	
	TVector3 axis, localY, preLocalX, newLocalX, v;
	TMatrix16 rotM, rotM_Y;


		             localY.SetSubtract( m_points[ 1 ], m_points[0]      );
	if( isClosed() ) localY.AddSubtract( m_points[ 0 ], m_points.back()  );
	t_getRotMat_YToV1( localY, axis, rotM );

	if( isClosed() ) for( int i = 0; i < vRes; ++i)
	{
		if( i==0){
			localY.SetSubtract( m_points[ 1 ], m_points[ 0 ]   );
			localY.AddSubtract( m_points[ 0 ], m_points.back() );
			t_getRotMat_YToV1( localY, axis, rotM );
		} else {
			localY.Set(0,0,0);
			if     ( i == (int) m_points.size()    ){ localY.AddSubtract( m_points[ 1 ]   , m_points[ 0 ]   ); 
				                                      localY.AddSubtract( m_points[ 0 ]   , m_points.back() ); }
			else if( i == (int) m_points.size() -1 ){ localY.AddSubtract( m_points[ 0 ]   , m_points.back() ); 
			                                          localY.AddSubtract( m_points.back() , m_points[i-1]   ); }
			else                                    { localY.AddSubtract( m_points[i+1]   , m_points[ i ]   );
													  localY.AddSubtract( m_points[ i ]   , m_points[i-1]   ); }

			preLocalX.Set( rotM.data[0], rotM.data[1], rotM.data[2]); // preLocalX = rotM * (0,1,0)
			if( !t_getRotMat_YToV1( localY, axis, rotM ) && localY.data[1] < 0) rotM.SetAsRotateX( M_PI );
			newLocalX.Set( rotM.data[0], rotM.data[1], rotM.data[2]); // preLocalX = rotM * (0,1,0)

			//localXを合わせる
			preLocalX.Add_CoefMultVec( -(preLocalX * localY), localY ); 
			preLocalX.Normalize_Self();

			rotM_Y.SetAsRotateY( t_getAngle_FixedAxis( newLocalX, preLocalX, localY) );
			rotM = rotM * rotM_Y;
		}	

		for( int j = 0; j < m_tubeRes; ++j)
		{
			v.Set( m_radius * cos( 2 * M_PI / (m_tubeRes -1) * j ), 0,
				   m_radius * sin( 2 * M_PI / (m_tubeRes -1) * j ));
			t_MatMultVec( v, rotM );
			m_tubeVerts[i][j].SetAddition((i==m_points.size() ) ?  m_points.front():m_points[i], v);
			m_tubeNorms[i][j] = v; 
			m_tubeNorms[i][j].Normalize_Self();
		}
	}
	else for( int i = 0; i < vRes; ++i)
	{
		if( i==0){
			localY.AddSubtract( m_points[ 1 ], m_points[0] );
			t_getRotMat_YToV1( localY, axis, rotM );
		} else {
			localY.Set(0,0,0);
   												localY.AddSubtract( m_points[ i ], m_points[i-1] );
			if( i != (int)m_points.size() - 1 ) localY.AddSubtract( m_points[i+1], m_points[ i ] );

			preLocalX.Set( rotM.data[0], rotM.data[1], rotM.data[2]); // preLocalX = rotM * (0,1,0)
			if( !t_getRotMat_YToV1( localY, axis, rotM ) && localY.data[1] < 0) rotM.SetAsRotateX( M_PI );
			newLocalX.Set( rotM.data[0], rotM.data[1], rotM.data[2]); // preLocalX = rotM * (0,1,0)

			//localXを合わせる
			preLocalX.Add_CoefMultVec( -(preLocalX * localY), localY );
			preLocalX.Normalize_Self();

			rotM_Y.SetAsRotateY( t_getAngle_FixedAxis( newLocalX, preLocalX, localY) );
			rotM = rotM * rotM_Y;
		}	

		for( int j = 0; j < m_tubeRes; ++j)
		{
			v.Set( m_radius * cos( 2 * M_PI / (m_tubeRes -1) * j ), 0,
				   m_radius * sin( 2 * M_PI / (m_tubeRes -1) * j ));
			t_MatMultVec( v, rotM );
			m_tubeVerts[i][j].SetAddition( m_points[i], v);
			m_tubeNorms[i][j] = v; 
			m_tubeNorms[i][j].Normalize_Self();
		}
	}
}

#endif



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Tstroke3D_tube///////////////////////////////////////////////////////////////////////////////////////////////


Tstroke3D_fiber::Tstroke3D_fiber()
{
	m_ROI       .clear();
	m_localRot  .clear();
	m_pointsOrig.clear();

	m_bPointFix.clear() ;
}


Tstroke3D_fiber::Tstroke3D_fiber(const vector<TVector3> &points ,int circleRes,double radius,bool isClosed) 
: Tstroke3D_tube( points, circleRes, radius, isClosed )
{
	fprintf( stderr, "const\n");

	//m_points等は初期化済み
	m_pointsOrig = m_points;
	m_ROI     .resize( m_points.size(), 0 );
	m_localRot.resize( m_points.size()    );
	m_bPointFix.resize(m_points.size(), 0 ); 
}

Tstroke3D_fiber::Tstroke3D_fiber(const Tstroke3D_fiber  &src)
: Tstroke3D_tube( src.m_points, src.getTubeResolution(), src.getTubeRadius(), src.isClosed())
{
	fprintf( stderr, " copy const\n");
	//m_points等は初期化済み
	m_pointsOrig = m_points;
	m_ROI     .resize( m_points.size(), 0 );
	m_localRot.resize( m_points.size()    );
	m_bPointFix = src.m_bPointFix; 
}


/*
leftは自分
left = rightとするだけでなく　*thisを返すことで、a=(b=c)を可能にする
*/
const Tstroke3D_fiber& Tstroke3D_fiber::operator=(const Tstroke3D_fiber& right)
{
	m_points   = right.m_points;

	setClosed		 ( right.isClosed()          );
	setTubeRadius    ( right.getTubeRadius()     );
	setTubeResolution( right.getTubeResolution() );
	createTubeFromAxis();
	
	m_pointsOrig = m_points;
	m_ROI     .resize( m_points.size(), 0 );
	m_localRot.resize( m_points.size()    );
	m_bPointFix = right.m_bPointFix; 

	return *this;
}



void Tstroke3D_fiber::SmoothingPtByPt( int ptId )
{
	if( ptId <= 0 || (int) m_points.size() -1 <= ptId ) return;
	m_bPointFix.resize(m_points.size(), 0); //外から直接push_backしたときにはbPointFixが短い状態になる
	if( m_bPointFix[ ptId ] ) return;

	TVector3 p;
	p.SetAdditionWithCoef( 0.333, m_points[ptId-1], 
	                       0.333 , m_points[ptId  ],
	                       0.333, m_points[ptId+1] );
	m_points[ptId  ] = p;
}


void Tstroke3D_fiber::SmoothingPtByPt_strong( int ptId )
{
	if( size() < 3 || ptId < 0 || (int) m_points.size() -1 < ptId ) return;
	m_bPointFix.resize(m_points.size(), 0); //外から直接push_backしたときにはbPointFixが短い状態になる
	if( m_bPointFix[ ptId ] ) return;

	if( ptId != 0 && (int) m_points.size() -1 != ptId )
	{
		m_points[ptId].SetAdditionWithCoef( 0.5, m_points[ptId-1], 0.5, m_points[ptId+1] );
	}
	else if( isClosed() && ptId == 0 )
	{
		m_points[ptId] .SetAdditionWithCoef( 0.5, m_points.back(), 0.5, m_points[ptId+1] );
	}
	else if( isClosed() && ptId == (int) m_points.size() -1  )
	{
		m_points[ptId] .SetAdditionWithCoef( 0.5, m_points.front(), 0.5, m_points[ptId-1] );
	}
}


void Tstroke3D_fiber::SmoothingPtByPt_keepLengthRatio( int ptId )
{
	if( ptId <= 0 || (int) m_points.size() -1 <= ptId ) return;
	m_bPointFix.resize(m_points.size(), 0); //外から直接push_backしたときにはbPointFixが短い状態になる
	if( m_bPointFix[ ptId ] ) return;

	double lpre = t_distance( m_points[ptId  ], m_points[ptId-1]);
	double lnex = t_distance( m_points[ptId  ], m_points[ptId+1]);
	TVector3 p;
	p.SetAdditionWithCoef( lnex/(lnex+lpre), m_points[ptId-1], 
	                       lpre/(lnex+lpre), m_points[ptId+1] );//内分点
	m_points[ptId] += p;
	m_points[ptId] *= 0.5;
}





	
void Tstroke3D_fiber::prepareDeform( int grabedPointIdx, double initialRoiRadius)
{
	fprintf( stderr, "start fiber deformation!!!!\n");
	m_pointsOrig = m_points;
	m_ROI     .clear(); m_ROI     .resize( m_points.size(), 0);
	m_localRot.clear(); m_localRot.resize( m_points.size()   );
	m_ROI[ grabedPointIdx ] = true;

	m_bPointFix.resize(m_points.size(), 0); //外から直接push_backしたときにはbPointFixが短い状態になる
	if( isClosed() ) updateROIforClosedFiber( grabedPointIdx, initialRoiRadius );
	else             updateROIforOpenFiber  ( grabedPointIdx, initialRoiRadius);
}

void Tstroke3D_fiber::prepareDeform( int grabedPointIdx)
{
	fprintf( stderr, "start fiber deformation!!!!\n");
	m_pointsOrig = m_points;
	m_ROI     .clear(); m_ROI     .resize( m_points.size(), 0);
	m_localRot.clear(); m_localRot.resize( m_points.size()   );
	m_ROI[ grabedPointIdx ] = true;
}

void Tstroke3D_fiber::finishDeform ()
{
	for( int i = 0; i < (int) m_points.size(); ++i) m_pointsOrig[i] = m_points[i];
	createTubeFromAxis();
	fprintf( stderr, "finish fiber deformation!!!!\n");
}

void Tstroke3D_fiber::moveVertex   ( int grabedPointIdx, const TVector3 &newPos){
	deformationFiber1( grabedPointIdx, newPos );
}




const int c_r[3][3] = 
{
	{  0,  2,  1},
	{  2,  0,  0},
	{  1,  0,  0}
};

const int pm_r[3][3] = 
{
	{ 0,  1, -1},
	{-1,  0,  1},
	{ 1, -1,  0}
};


double Tstroke3D_fiber::w1 = 20;
double Tstroke3D_fiber::w2 = 100;
double Tstroke3D_fiber::w3 = 20;
double Tstroke3D_fiber::w4 = 20;
double Tstroke3D_fiber::ROI_LengthRate = 0.5;

/*
Curve deformation of "the fiber mesh"
http://www-ui.is.s.u-tokyo.ac.jp/~takeo/papers/fibermesh.pdf

arg min  Σ|| L(vi) - ri*Ri*δi ||^2  +  Σ|| vi - v'i ||^2   +   Σ|| ri*Ri - rj*Rj ||^2   + Σ||ri*Ri - R'i||^2 
vi, ri    all                            i∈c1                      i∈E                       i∈c2
  
           pointNum*3                 + constNum*3            +   edgeNum*9                 + constNum2 * 9

	numOfUnknown = 3 * pointNum + 3 * pointNum
*/
void Tstroke3D_fiber::deformationFiber( int grabedPointIdx, const TVector3 &point)
{

	m_bPointFix.resize(m_points.size(), 0); //外から直接push_backしたときにはbPointFixが短い状態になる

	if( isClosed() ) updateROIforClosedFiber( grabedPointIdx, point );
	else             updateROIforOpenFiber  ( grabedPointIdx, point );


	//check ROI, if all vertices are controled by the user the system just apply Rigid trans
	const int pSize = (int) m_pointsOrig.size();
	bool allVertsAreActive = true;
	for( int i = 0; i < pSize; ++i) if( !m_ROI[i] ) allVertsAreActive = false;

	if( allVertsAreActive )
	{
		TVector3 rigidTrans(point - m_pointsOrig[grabedPointIdx] );
		for( int i = 0; i < pSize; ++i) m_points[i].SetAddition( rigidTrans, m_pointsOrig[i] );
		createTubeFromAxis();
		return;
	}

	int constIndx1 = -1, constIndx2 = -1;
	{
		int pivot = 0;
		for( ; pivot < pSize; ++pivot) if( !m_ROI[ pivot ] )
			if( m_ROI[ ( pivot ==    0    ) ? pSize-1 : pivot -1 ] || 
				m_ROI[ ( pivot == pSize- 1) ?    0    : pivot +1 ] ) { constIndx1 = pivot; pivot++; break; }
		for( ; pivot < pSize; ++pivot) if( !m_ROI[ pivot ] )
			if( m_ROI[ ( pivot ==   0     ) ? pSize -1 : pivot -1 ] || 
				m_ROI[ ( pivot == pSize- 1) ?    0     : pivot +1 ] ){ constIndx2 = pivot; pivot++; break; }
	}

	//Calculate At/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//At size = 3 * n_point (1) + 9 * n_edge + n_constP * 3 + n_constRot * 9 + constraint//////////////////////////////////
	int numRotConsts = (constIndx1 == -1 || constIndx2 == -1 ) ? 0 :
					   (constIndx1 == constIndx2             ) ? 1 : 2 ;
	int numPosConsts = 0;
	for( int i = 0; i < pSize; ++i) if( !m_ROI[i] ) ++numPosConsts;

	const int nRows_At =  6 * pSize;           //縦の長さ
	//横の長さ
	const int nCols_At =  isClosed() ?         
		                  3 * pSize            + 
		                  9 * pSize            +   //edgeが一本多い
						  3 + 3 * numPosConsts + 
						  9     * numRotConsts      
						  :
		                  3 * pSize            + //laplacian constraint          |L(vi) - ri*Ri*δi|
		                  9 * ( pSize -1 )     + //neighbor rotation constraint  |ri*Ri -   rj*Rj  | 
						  3 + 3 * numPosConsts + //positional constraint         |  vi  -     vi'  |
						  9 *     numRotConsts ; //rotation constraint           |ri*Ri -     Ri'  |横の長さ
	
	
	TSparseMatrix At( nRows_At, nCols_At);
	double *b   = new double[nCols_At];
	double *Atb = new double[nRows_At];
	int eqIndex = 0;

	//1 laplacian constraints |L(vi) - ri*Ri*δi| /////////////////////////////////////////////
	{
		TVector3 abc, delta;
		int idxPiv, idxPre;
		for( int i = 0; i < pSize; ++i){
			if( isClosed() ){
				idxPiv = i;
				idxPre = (i==0) ? pSize-1 : i-1;
			}else{
				idxPiv = (i==0) ? 0 :  i ;
				idxPre = (i==0) ? 1 : i-1;
			}

			delta.SetSubtract( m_pointsOrig[idxPiv], m_pointsOrig[idxPre] );
			t_MatMultVec( m_localRot[i], delta, abc);

			for( int xyz = 0; xyz < 3; ++xyz)
			{			
				At.push_back( idxPiv * 6 + xyz, eqIndex,    w1 );
				At.push_back( idxPre * 6 + xyz, eqIndex,   -w1 );
			
				if( pm_r[xyz][0] != 0) At.push_back( idxPiv * 6 + 3, eqIndex,  - w1 * pm_r[xyz][0] * abc.data[ c_r[xyz][0] ]  );
				if( pm_r[xyz][1] != 0) At.push_back( idxPiv * 6 + 4, eqIndex,  - w1 * pm_r[xyz][1] * abc.data[ c_r[xyz][1] ]  );
				if( pm_r[xyz][2] != 0) At.push_back( idxPiv * 6 + 5, eqIndex,  - w1 * pm_r[xyz][2] * abc.data[ c_r[xyz][2] ]  );
				b[ eqIndex ] = w1 * abc.data[xyz];
				++eqIndex;
			}
		}	
	}

	//2  Rotation Σ|| ri*Ri - rj*Rj ||^2  ////////////////////////////////////////////////
	{
		int edgeID = isClosed() ? 0 : 1;
		int p0idx, p1idx;
		for( ; edgeID < pSize; ++edgeID)
		{
			p0idx = edgeID;
			p1idx = ( edgeID == 0 ) ? pSize -1 : edgeID - 1; //isClosedでなければ edgeID==0とはならない

			const TMatrix9 &R0 = m_localRot[ p0idx ]; 
			const TMatrix9 &R1 = m_localRot[ p1idx ];

			for( int ii = 0; ii < 3; ++ii )
			for( int jj = 0; jj < 3; ++jj )
			{
				if( pm_r[ii][0] != 0) At.push_back( p0idx * 6 + 3, eqIndex,  w2 * (  pm_r[ii][0] * R0.data[ jj*3 + c_r[ii][0] ] )  );
				if( pm_r[ii][1] != 0) At.push_back( p0idx * 6 + 4, eqIndex,  w2 * (  pm_r[ii][1] * R0.data[ jj*3 + c_r[ii][1] ] ) );
				if( pm_r[ii][2] != 0) At.push_back( p0idx * 6 + 5, eqIndex,  w2 * (  pm_r[ii][2] * R0.data[ jj*3 + c_r[ii][2] ] ) );

				if( pm_r[ii][0] != 0) At.push_back( p1idx * 6 + 3, eqIndex,  w2 * ( -pm_r[ii][0] * R1.data[ jj*3 + c_r[ii][0] ] )  );
				if( pm_r[ii][1] != 0) At.push_back( p1idx * 6 + 4, eqIndex,  w2 * ( -pm_r[ii][1] * R1.data[ jj*3 + c_r[ii][1] ] )  );
				if( pm_r[ii][2] != 0) At.push_back( p1idx * 6 + 5, eqIndex,  w2 * ( -pm_r[ii][2] * R1.data[ jj*3 + c_r[ii][2] ] )  );
				b[eqIndex] = w2 * ( - R0.data[ jj*3 + ii ] + R1.data[ jj*3 + ii ] ); 
				++eqIndex;
			}
		}
	}

	//3 Position Constraint   3equation for grabed and m equations for not moved points///////////////////////////////////////////////
	for( int xyz = 0; xyz < 3; ++xyz ){
		At.push_back( grabedPointIdx * 6 + xyz, eqIndex,  w3 ); 
		b[eqIndex] = w3 * point.data[xyz]; 
		++eqIndex; 
	}
	for( int cPidx = 0; cPidx < pSize; ++cPidx) if( !m_ROI[cPidx] )
	{
		for( int xyz = 0; xyz < 3; ++xyz )
		{
			At.push_back( cPidx * 6 + xyz, eqIndex, w3 ); 
			b[eqIndex] = w3 * m_pointsOrig[ cPidx ].data[xyz]; 
			++eqIndex; 
		}
	}

	//4  |ri Ri - Ri'|~2  num of equation = 9 since R is 3 * 3 matrix/////////////////////////////////////////////
	for( int C2 = 0; C2 < numRotConsts; ++C2)
	{
		int cPidx = (C2==0)? constIndx1:constIndx2 ;

		const TMatrix9 &R = m_localRot[ cPidx ]; 

		for( int ii = 0; ii < 3; ++ii )
		for( int jj = 0; jj < 3; ++jj )
		{
			if( pm_r[ii][0] != 0) At.push_back( cPidx * 6 + 3,  eqIndex, w4 * ( pm_r[ii][0] * R.data[ jj*3 + c_r[ii][0] ] ) );
			if( pm_r[ii][1] != 0) At.push_back( cPidx * 6 + 4,  eqIndex, w4 * ( pm_r[ii][1] * R.data[ jj*3 + c_r[ii][1] ] ) );
			if( pm_r[ii][2] != 0) At.push_back( cPidx * 6 + 5,  eqIndex, w4 * ( pm_r[ii][2] * R.data[ jj*3 + c_r[ii][2] ] ) );
			b[eqIndex] = w4 * ( - R.data[ jj*3 + ii ] + ( (ii==jj)? 1 : 0 ) ); 
			++eqIndex;
		}
	}

	//calc AtA and Atb//////////////////////////////
	TSparseMatrix AtA( nRows_At, nRows_At );
	int nonZero = At.calcMmultMt( AtA );
	At.multVec( nCols_At, b, nRows_At, Atb );

	double*   result = new double[ nRows_At    ];

	memset( result, 0, sizeof( double) * nRows_At );
	AtA.createFieldForSolveLinearSystem(nonZero);
	AtA.solveLinearEquationUmfPack( Atb, result);

	//get deformed vertex and update Rotations
	TMatrix9 rot;
	for( int i = 0; i < pSize; ++i)
	{
		int piv = i * 6;
		m_points[i].Set( result[ piv+0 ], result[ piv+1 ], result[ piv+2 ] );

		rot.Set(          1      , -result[piv + 5],   result[piv + 4],
			      result[piv + 5],         1       ,  -result[piv + 3],
			     -result[piv + 4],  result[piv + 3],          1       );
		t_extractRotationTermFromMatrix( rot * m_localRot[i], m_localRot[i]);
	}

	createTubeFromAxis();

	delete[] result;
	delete[] Atb;
	delete[] b;
}






/*
毎回roi分だけしかmatrixを生成しない実装
*/
void Tstroke3D_fiber::deformationFiber1( int grabedPointIdx, const TVector3 &point)
{
	m_bPointFix.resize(m_points.size(), 0); //外から直接push_backしたときにはbPointFixが短い状態になる
	if( isClosed() ) updateROIforClosedFiber( grabedPointIdx, point );
	else             updateROIforOpenFiber  ( grabedPointIdx, point );
	//check ROI, if all vertices are controled by the user the system just apply Rigid trans
	int roiPointSize = 0;
	for( int i = 0; i < (int)m_points.size(); ++i) if( m_ROI[i] ) ++roiPointSize;


	if(      roiPointSize == 0 ) return;
	else if( roiPointSize == 1 ){
		m_points[ grabedPointIdx ] = point;
		createTubeFromAxis();
		return;
	}
	else if( roiPointSize == m_points.size() ){
		TVector3 rigidTrans( point - m_pointsOrig[grabedPointIdx] );
		for( int i = 0; i < (int) m_points.size(); ++i) m_points[i].SetAddition( rigidTrans, m_pointsOrig[i] );
		createTubeFromAxis();
		return;
	}


	//all activeでないため、端点を固定したdeformationになる
	vector< int > vmap_FULtoROI( m_points.size(), -1  );
	vector< int > vmap_ROItoFUL( roiPointSize   , -1  ); 

	if( m_ROI.back() && m_ROI.front() ) //front - backをまたぐroi region
	{
		int roiStartIdx = -1;
		for( int i = (int)m_points.size()-2; i >= 0; --i) if( !m_ROI[i] ) { roiStartIdx = i+1; break; }
		
		int idx = 0;
		for( int i = roiStartIdx; i < (int)m_points.size(); ++i){
			vmap_FULtoROI[ i   ] = idx;
			vmap_ROItoFUL[ idx ] = i  ;
			++idx;
		}
		for( int i = 0; m_ROI[i]; ++i) {
			vmap_FULtoROI[ i   ] = idx;
			vmap_ROItoFUL[ idx ] = i ;
			++idx;
		}
	}else{
		int idx = 0;
		for( int i = 0; i < (int) m_points.size(); ++i) if( m_ROI[i] ){
			vmap_FULtoROI[ i   ] = idx;
			vmap_ROItoFUL[ idx ] = i ;
			idx++;
		}
	}


	//このroi pointsの両端点を固定した変形を行う
	const int pSize = roiPointSize  ;

	//Calculate At/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//At size = 3 * n_point (1) + 9 * n_edge + n_constP * 3 + n_constRot * 9 + constraint//////////////////////////////////
	const int nRows_At =  6 * pSize;   
	const int nCols_At =  3 *   pSize      + //laplacian constraint          |L(vi) - ri*Ri*δi|
		                  9 * ( pSize -1 ) + //neighbor rotation constraint  |ri*Ri -   rj*Rj  | 
						  3 * ( 1 + 2    ) + //positional constraint         |  vi  -     vi'  | //両端点 + grab = 2点
						  9 *       2      ; //rotation constraint           |ri*Ri -     Ri'  |横の長さ 両端点

	TSparseMatrix At( nRows_At, nCols_At);
	double *b      = new double[nCols_At];
	double *Atb    = new double[nRows_At];
	double* result = new double[nRows_At];

	int eqIndex = 0;


	//1 laplacian constraints |L(vi) - ri*Ri*δi| /////////////////////////////////////////////
	{
		TVector3 abc, delta;
		for( int i = 0; i < pSize; ++i)
		{
			int idxPiv = (i==0) ? 0 :  i ;
			int idxPre = (i==0) ? 1 : i-1;
			
			delta.SetSubtract( m_pointsOrig[ vmap_ROItoFUL[idxPiv] ], m_pointsOrig[ vmap_ROItoFUL[idxPre] ] );
			t_MatMultVec( m_localRot[ vmap_ROItoFUL[ i ] ], delta, abc);

			for( int xyz = 0; xyz < 3; ++xyz)
			{			
				At.push_back( idxPiv * 6 + xyz, eqIndex,    w1 );
				At.push_back( idxPre * 6 + xyz, eqIndex,   -w1 );
			
				if( pm_r[xyz][0] != 0) At.push_back( idxPiv * 6 + 3, eqIndex,  - w1 * pm_r[xyz][0] * abc.data[ c_r[xyz][0] ]  );
				if( pm_r[xyz][1] != 0) At.push_back( idxPiv * 6 + 4, eqIndex,  - w1 * pm_r[xyz][1] * abc.data[ c_r[xyz][1] ]  );
				if( pm_r[xyz][2] != 0) At.push_back( idxPiv * 6 + 5, eqIndex,  - w1 * pm_r[xyz][2] * abc.data[ c_r[xyz][2] ]  );
				b[ eqIndex ] = w1 * abc.data[xyz];
				++eqIndex;
			}
		}	
	}
	//fprintf( stderr, " eqIndex %d / %d \n", eqIndex, pSize * 3);


	//2  Rotation Σ|| ri*Ri - rj*Rj ||^2  ////////////////////////////////////////////////
	{
		for( int edgeID = 1; edgeID < pSize; ++edgeID)
		{
			int p0idx =   edgeID    ;
			int p1idx =   edgeID - 1; 

			const TMatrix9 &R0 = m_localRot[ vmap_ROItoFUL[p0idx] ]; 
			const TMatrix9 &R1 = m_localRot[ vmap_ROItoFUL[p1idx] ];

			for( int ii = 0; ii < 3; ++ii )
			for( int jj = 0; jj < 3; ++jj )
			{
				if( pm_r[ii][0] != 0) At.push_back( p0idx * 6 + 3, eqIndex,  w2 * (  pm_r[ii][0] * R0.data[ jj*3 + c_r[ii][0] ] )  );
				if( pm_r[ii][1] != 0) At.push_back( p0idx * 6 + 4, eqIndex,  w2 * (  pm_r[ii][1] * R0.data[ jj*3 + c_r[ii][1] ] ) );
				if( pm_r[ii][2] != 0) At.push_back( p0idx * 6 + 5, eqIndex,  w2 * (  pm_r[ii][2] * R0.data[ jj*3 + c_r[ii][2] ] ) );

				if( pm_r[ii][0] != 0) At.push_back( p1idx * 6 + 3, eqIndex,  w2 * ( -pm_r[ii][0] * R1.data[ jj*3 + c_r[ii][0] ] )  );
				if( pm_r[ii][1] != 0) At.push_back( p1idx * 6 + 4, eqIndex,  w2 * ( -pm_r[ii][1] * R1.data[ jj*3 + c_r[ii][1] ] )  );
				if( pm_r[ii][2] != 0) At.push_back( p1idx * 6 + 5, eqIndex,  w2 * ( -pm_r[ii][2] * R1.data[ jj*3 + c_r[ii][2] ] )  );
				b[eqIndex] = w2 * ( - R0.data[ jj*3 + ii ] + R1.data[ jj*3 + ii ] ); 
				++eqIndex;
			}
		}
	}
	//fprintf( stderr, " eqIndex %d / %d \n", eqIndex, pSize * 3 + 9 * (pSize-1) );

	//3 Position Constraint   3equation for grabed and m equations for not moved points///////////////////////////////////////////////
	//grab + 両端点
	int constIdx = vmap_FULtoROI[ grabedPointIdx ];
	At.push_back( constIdx*6 + 0, eqIndex, w3 ); b[eqIndex] = w3 * point.data[0];  ++eqIndex;
	At.push_back( constIdx*6 + 1, eqIndex, w3 ); b[eqIndex] = w3 * point.data[1];  ++eqIndex;
	At.push_back( constIdx*6 + 2, eqIndex, w3 ); b[eqIndex] = w3 * point.data[2];  ++eqIndex;
		 
	constIdx = vmap_ROItoFUL[0];
	At.push_back( 0*6 + 0, eqIndex, w3 ); b[eqIndex] = w3 * m_pointsOrig[constIdx].data[0];  ++eqIndex;
	At.push_back( 0*6 + 1, eqIndex, w3 ); b[eqIndex] = w3 * m_pointsOrig[constIdx].data[1];  ++eqIndex;
	At.push_back( 0*6 + 2, eqIndex, w3 ); b[eqIndex] = w3 * m_pointsOrig[constIdx].data[2];  ++eqIndex;
	
	constIdx = vmap_ROItoFUL[ pSize-1 ];
	At.push_back((pSize-1)*6 + 0, eqIndex, w3 ); b[eqIndex] = w3 * m_pointsOrig[constIdx].data[0];  ++eqIndex;
	At.push_back((pSize-1)*6 + 1, eqIndex, w3 ); b[eqIndex] = w3 * m_pointsOrig[constIdx].data[1];  ++eqIndex;
	At.push_back((pSize-1)*6 + 2, eqIndex, w3 ); b[eqIndex] = w3 * m_pointsOrig[constIdx].data[2];  ++eqIndex;

	//fprintf( stderr, " eqIndex %d\n", eqIndex);

	//4  |ri Ri - Ri'|~2  num of equation = 9 since R is 3 * 3 matrix/////////////////////////////////////////////
	for( int C2 = 0; C2 < 2; ++C2)
	{
		int cPidx = (C2==0) ? 0: pSize-1 ;

		const TMatrix9 &R = m_localRot[ vmap_ROItoFUL[ cPidx ] ]; 

		for( int ii = 0; ii < 3; ++ii )
		for( int jj = 0; jj < 3; ++jj )
		{
			if( pm_r[ii][0] != 0) At.push_back( cPidx * 6 + 3,  eqIndex, w4 * ( pm_r[ii][0] * R.data[ jj*3 + c_r[ii][0] ] ) );
			if( pm_r[ii][1] != 0) At.push_back( cPidx * 6 + 4,  eqIndex, w4 * ( pm_r[ii][1] * R.data[ jj*3 + c_r[ii][1] ] ) );
			if( pm_r[ii][2] != 0) At.push_back( cPidx * 6 + 5,  eqIndex, w4 * ( pm_r[ii][2] * R.data[ jj*3 + c_r[ii][2] ] ) );
			b[eqIndex] = w4 * ( - R.data[ jj*3 + ii ] + ( (ii==jj)? 1 : 0 ) ); 
			++eqIndex;
		}
	}

	//fprintf( stderr, " eqIndex %d\n", eqIndex);


	//calc AtA and Atb//////////////////////////////
	TSparseMatrix AtA( nRows_At, nRows_At );
	int nonZero = At.calcMmultMt( AtA );
	AtA.createFieldForSolveLinearSystem(nonZero);

	At.multVec( nCols_At, b, nRows_At, Atb );


	memset( result, 0, sizeof( double) * nRows_At );
	AtA.solveLinearEquationUmfPack( Atb, result);



	//get deformed vertex and update Rotations
	TMatrix9 rot;
	for( int i = 0; i < pSize; ++i)
	{
		TVector3 &p = m_points  [ vmap_ROItoFUL[i] ];
		TMatrix9  &R = m_localRot[ vmap_ROItoFUL[i] ];
		int piv = i * 6;
		p.Set( result[ piv+0 ], result[ piv+1 ], result[ piv+2 ] );

		rot.Set(          1      , -result[piv + 5],   result[piv + 4],
			      result[piv + 5],         1       ,  -result[piv + 3],
			     -result[piv + 4],  result[piv + 3],          1       );


		t_extractRotationTermFromMatrix( rot * R, R);
	}

	createTubeFromAxis();


	delete[] result;
	delete[] Atb;
	delete[] b;
}








//ROIはfix pointを超えて成長しない
void Tstroke3D_fiber::updateROIforOpenFiber(int grabedPointIdx, const TVector3 &newPos)
{
	updateROIforOpenFiber( grabedPointIdx, ROI_LengthRate * t_distance( newPos, m_pointsOrig[grabedPointIdx]) );
}




void Tstroke3D_fiber::updateROIforClosedFiber  (int grabedPointIdx, const TVector3 &newPos)
{
	updateROIforClosedFiber(grabedPointIdx, ROI_LengthRate * t_distance( newPos, m_pointsOrig[grabedPointIdx] ));
}





void Tstroke3D_fiber::updateROIforOpenFiber  (int grabedPointIdx, double roiLength)
{
	m_ROI[ grabedPointIdx ] = true;
	
	const int pSize = (int) m_pointsOrig.size();
	double sumOfLength = 0;
	for(int pivot = grabedPointIdx + 1 ; pivot < pSize && !m_bPointFix[pivot]; ++pivot)
	{
		sumOfLength += t_distance( m_pointsOrig[pivot] , m_pointsOrig[pivot-1] );
		if( sumOfLength > roiLength ) break;
		m_ROI[ pivot ] = true;
	}

	sumOfLength = 0;
	for(int pivot = grabedPointIdx -1 ; pivot >= 0 && !m_bPointFix[pivot]; --pivot)
	{
		sumOfLength += t_distance( m_pointsOrig[pivot] , m_pointsOrig[pivot+1] );
		if( sumOfLength > roiLength ) break;
		m_ROI[ pivot ] = true;
	}
}
void Tstroke3D_fiber::updateROIforClosedFiber(int grabedPointIdx, double roiLength)
{
	m_ROI[ grabedPointIdx ] = true;
	
	const int pSize = (int) m_points.size();
	double sumOfLength = 0;

	//calc ROI for circular points
	int    pivot  = grabedPointIdx + 1;
	while( pivot != grabedPointIdx )
	{
		if( pivot > pSize-1 ) pivot = 0; 
		if( m_bPointFix[pivot] ) break;

		sumOfLength += (pivot == 0)? t_distance( m_pointsOrig.front(), m_pointsOrig.back()   ):
			                         t_distance( m_pointsOrig[pivot] , m_pointsOrig[pivot-1] );
		if( sumOfLength > roiLength ) break;
		m_ROI[ pivot ] = true;
		++pivot;
	}

	pivot       = grabedPointIdx -1;
	sumOfLength = 0;
	while( pivot != grabedPointIdx )
	{
		if( pivot < 0 ) pivot = pSize - 1;
		if( m_bPointFix[pivot] ) break;

		sumOfLength += (pivot == pSize - 1)? t_distance( m_pointsOrig.front(), m_pointsOrig.back()   ):
			                                 t_distance( m_pointsOrig[pivot] , m_pointsOrig[pivot+1] );
		if( sumOfLength > roiLength ) break; 
		m_ROI[ pivot ] = true;
		--pivot;
	}
}
