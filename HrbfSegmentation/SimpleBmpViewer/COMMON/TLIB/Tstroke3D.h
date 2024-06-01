#pragma once

#include "tmath.h"



class Tstroke3D
{
private:
	bool m_isClosed;
public:
	vector<TVector3> m_points;

	virtual ~Tstroke3D(void){}
	Tstroke3D(                                  bool m_isClosed=false):m_isClosed(    m_isClosed){ m_points.clear()       ; }
	Tstroke3D(const Tstroke3D &src                                   ):m_isClosed(src.m_isClosed){ m_points = src.m_points; }
	Tstroke3D(const vector<TVector3> &points, bool m_isClosed=false):m_isClosed(    m_isClosed){ m_points =       points; }

	double	Length(                     ) const;
	double	Length(int sIndex,int eIndex) const;
	void	devideEquals(int n);
	void	Smoothing(void);
	void    Smoothing( int k ){for(int i=0; i<k; ++i) Smoothing();}
	void	SmoothingUsingBSpline();

	inline bool isClosed() const {return m_isClosed;}
	inline void setClosed(bool tf){m_isClosed = tf; }
	inline void push_back( const TVector3 &p){ m_points.push_back( p ) ;}
	inline TVector3& back (){ return m_points.back ();}
	inline TVector3& front(){ return m_points.front();}
	inline int         size ()const{ return (int) m_points.size();}

	inline virtual void clear    ()   { m_points.clear()        ;}
	inline void resize(size_t size){ m_points.resize( size);}
	inline void invertStroke(){ std::reverse( m_points.begin(), m_points.end() ); }
	inline void Trace() const { for(int i = 0; i <(int) m_points.size(); ++i) m_points[i].Trace(); }
	inline void mulMatrixToStroke(const TMatrix16 &M){
		for(int i = 0;i < (int)m_points.size(); ++i) t_MatMultVec( m_points[i], M );                 
	}

	double distanceFromEyeRayToPoints(const TVector3 &eyeP, const TVector3 &rayDir, int &nearestPointIdx) const;
	double distanceFromEyeRayToCurve (const TVector3 &eyeP, const TVector3 &rayDir, int &nearestPointIdx) const;


	
/*
	double distanceFromEyeRayToStroke(const TVector3 &eyeP, 
		                              const TVector3 &rayDir, 
											double     &distFromEye, 
											TVector3 &nearestPoint );

*/
#ifdef USE_TOGL
	void drawLineStrip( float width) const;
#endif


	inline       TVector3& operator[](int idx)      { return m_points[idx]; }
	inline const TVector3& operator[](int idx)const { return m_points[idx]; }
};



//Tstroke with Generalized cylinder(tube)
class Tstroke3D_tube : public Tstroke3D
{
private:
	vector< vector< TVector3 > > m_tubeVerts;
	vector< vector< TVector3 > > m_tubeNorms;
	int    m_tubeRes;
	double m_radius;

public:
	virtual ~Tstroke3D_tube(){}
	Tstroke3D_tube(                            int tubeRes=6,double radius=0.1,bool isClosed=false);
	Tstroke3D_tube(const vector<TVector3> &ps ,int tubeRes=6,double radius=0.1,bool isClosed=false);
	Tstroke3D_tube(const Tstroke3D_tube &src);

	inline void clear(){ 
		m_points   .clear();
		m_tubeVerts.clear();
		m_tubeNorms.clear();
	}


#ifdef USE_TOGL
	void drawTube()const ;
#endif
	//both for initial creation and reflesh the shape
	void createTubeFromAxis();
	inline double getTubeRadius    ()        const { return m_radius    ; }
	inline int    getTubeResolution()        const { return m_tubeRes   ; }
	inline void   setTubeRadius    ( double radius){ m_radius  = radius ; }
	inline void   setTubeResolution( int   tubeRes){ m_tubeRes = tubeRes; }
};




//Tstroke with curved fiber of "the fiber mesh"
class Tstroke3D_fiber : public Tstroke3D_tube
{
	vector< TVector3 > m_pointsOrig;//original set of points
	vector< int      > m_ROI;
	vector< TMatrix9 > m_localRot;

	// fixされているかどうか, 外からm_pointsに直接push_backした場合には，m_pointsとの長さが合わない可能性があるので、このfieldにアクセスするfunctionは注意が必要
	vector< short > m_bPointFix; 

public:
	static double w1 ;
	static double w2 ;
	static double w3 ;
	static double w4 ;
	static double ROI_LengthRate;

public:
	virtual ~Tstroke3D_fiber(){}
	Tstroke3D_fiber();
	Tstroke3D_fiber(const Tstroke3D_fiber     &src);
	Tstroke3D_fiber(const vector<TVector3> &points ,int circleRes=6,double radius=0.1,bool isClosed=false);
	const Tstroke3D_fiber& operator=(const Tstroke3D_fiber& right);



	//deformation interface
	void prepareDeform( int grabedPointIdx);
	void prepareDeform( int grabedPointIdx, double initialRoiRadius);
	void moveVertex   ( int grabedPointIdx, const TVector3 &newPos);
	void finishDeform ();


	void addFixPoint   ( int idx ){ if( m_bPointFix.size() != m_points.size() ) m_bPointFix.resize(m_points.size(), 0);        m_bPointFix[idx] = 1; }
	void removeFixPoint( int idx ){ if( m_bPointFix.size() != m_points.size() ) m_bPointFix.resize(m_points.size(), 0);        m_bPointFix[idx] = 0; }
	inline short isFixed(int idx ){ if( m_bPointFix.size() != m_points.size() ) m_bPointFix.resize(m_points.size(), 0); return m_bPointFix[idx]; }
	void removeAllFixPoint()   { m_bPointFix.clear(); m_bPointFix.resize( size(), 0 ); }



	void SmoothingPtByPt( int ptId );
	void SmoothingPtByPt_strong( int ptId );
	void SmoothingPtByPt_keepLengthRatio( int ptId );

	inline bool isInCurrentROI( int i)
	{
		if( i<0 || (int)m_ROI.size()-1 <i) return false;
		return m_ROI[i]?true:false;
	}

private:
	void updateROIforOpenFiber  (int grabedPointIdx, const TVector3 &newPos);
	void updateROIforClosedFiber(int grabedPointIdx, const TVector3 &newPos);
	void updateROIforOpenFiber  (int grabedPointIdx, double roiLength);
	void updateROIforClosedFiber(int grabedPointIdx, double roiLength);

	void deformationFiber       (int grabedPointIdx, const TVector3 &newPos);
	void deformationFiber1      (int grabedPointIdx, const TVector3 &newPos);
};










///////////////////////////////////////////////////////////////

//3-degree B-Spline weighting function 
inline double t_bSpline3(double t){
	if(t < 0  ) t = -t;
	if     (t < 1.0) return 0.5 *t*t*t - t*t  + 4.0/6.0;
	else if(t < 2.0){
		double k = t - 2.0;
		return  -k*k*k / 6.0;
	}else
		return 0;
}

//set b spline point p(t) in pos
//t should be 0 <= t <= N-1
inline void t_bSpline3_pos(const vector<TVector3> &CPs, double t , TVector3 &pos)
{	
	const int N = (int)CPs.size();

	pos.Set( 0,0,0 );
	//t_bSpline3は [-2, 2]の外ではゼロ
	for( int i = -2; i < N + 2; ++i)
	{
		if( fabs( t - i ) > 2 ) continue;

		int index = ( i < 0  ) ? 0   :
					( i > N-1) ? N-1 : i ;

		pos.Add_CoefMultVec( t_bSpline3(t - i), CPs[index] ); 
	}
}

//pointsを制御点としたBsplineを生成し、n個の点をサンプリングしてtargetに挿入
inline void t_stroke_Smooth_BSpline3( const vector<TVector3> &points, vector<TVector3> &target , int NumOfPoint)
{
	const int N = (int) points.size();
	if(N < 2 || NumOfPoint < 2){ fprintf( stderr, "strange input (@bspline smoothing)\n"); target = points; return;}
	
	// -1 <= t <= N-1 + 1
	double stepD = (N-1 + 2 ) / (double)(NumOfPoint -1); //３点なら２等分
	
	target.resize( NumOfPoint );
	for( int i=0; i<NumOfPoint; ++i)
		t_bSpline3_pos(points, i * stepD - 1, target[i]);
}

void t_stroke_Smooth_BSpline3Interpolation(const vector<TVector3> &points, vector<TVector3> &target, bool isClosed, int NumOfPoint);
