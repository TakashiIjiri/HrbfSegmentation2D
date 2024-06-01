#pragma once
#include "TOGLImage.h"
#include "TTriangleMesh.h"

class RBF_PosCP{
public:
	double   val_c, val_s, val_j;
	double p[5];//xyrgb‚Ì5ŽŸŒ³
	RBF_PosCP( const TVector2 &_p, double r, double g, double b, double _valC, double _valS, double _valJ)
	{
		val_c = _valC;
		val_s = _valS;
		val_j = _valJ;
		p[0] = _p[0];
		p[1] = _p[1];
		p[2] = r;
		p[3] = g;
		p[4] = b;
	}
};

class RBF_GradCP{
public:
	double dir[3]; // dir(xyr) 3ŽŸŒ³
	double p  [5]; //xyrgb‚Ì5ŽŸŒ³
	int    pixX,pixY;
	RBF_GradCP( const TVector2 &_p, double r, double g, double b, TVector3 &grad, int x,int y){
		p[0] = _p[0];
		p[1] = _p[1];
		p[2] = r;
		p[3] = g;
		p[4] = b;
		dir[0]  = grad[0];
		dir[1]  = grad[1];
		dir[2]  = grad[2];
		pixX = x;
		pixY = y;
	}
};




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                constraints 
//                pos    pos + ¤c    pos + ¤s   pos + ¤j   
//	     color    ›        ›           ~           ~      
// space scale    ›        ›           ›           ~       
//       joint    ›@      ›           ›           › 
//
//Solver --> leastSquare or HRBF


class RBFManager
{
	RBFManager(void){}
public:
	~RBFManager(void){}

	inline static RBFManager* getInst( ){ static RBFManager p; return &p; }
	//static double m_scalePitch ; 
	//static double m_colorPitch ;
	static double m_approxCoef ;
	static int    m_polynomMode;
	static int    m_basisFuncMode;

	static void  RBFs_color     ( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs,								     double cPitch, double sPitch, float *scholarField);
	static void  RBFs_scale     ( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs,								     double cPitch, double sPitch, float *scholarField);
	static void  RBFs_joint     ( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs,								     double cPitch, double sPitch, float *scholarField, TTriangleMesh *zeroMesh);
	static void  RBFs_color_graC( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField);
	static void  RBFs_scale_graS( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField);
	static void  RBFs_joint_graC( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField);
	static void  RBFs_joint_graJ( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField);

	static void HRBFs_color_graC( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField);
	static void HRBFs_scale_graS( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField,                    TTriangleMesh *zeroMesh);
	static void HRBFs_joint_graJ( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField, double colDirCoef, TOGL3DImage4 *volume, TTriangleMesh *zeroMesh);

};



inline void rbfKernel1D( const double &d0, double &val, double &g0, double &h0 )
{
	const int &m=RBFManager::getInst()->m_basisFuncMode;

	double r   = fabs( d0 );
	double A, B;
	if(       m == 1){ val = r; //x
					   A   = (r != 0 ) ?  1/r       : 10000000; 
					   B   = (r != 0 ) ? -1/(r*r*r) : 10000000; 
	}else if( m == 2){ val = r*r*r;//x^3
					   A   = 3 * r;
					   B   = (r == 0) ? 0 : 3/r;
	}else if( m == 3){ val = r*r*r*r*r ; //x^5 
					   A   =  5 * r*r*r;
					   B   = 15 * r    ;
	}else if( m == 4){ val = (r==0) ? 0 : r*r*log(r)    ;//x^2 log(x)
					   A   = (r==0) ? 0 :   2*log(r) + 1;
					   B   = (r==0) ? 2 :   2 / (r*r)   ;

	}else if( m == 5){ val = (r==0) ? 0 : r*r*r*r*log(r);//x^4 log(x)
					   A   = (r==0) ? 0 : r*r*( 4*log(r) + 1 );
					   B   = (r==0) ? 0 : 8*log(r) + 6; 
					   
	}else {            val = (r>1)  ? 0 : (1-r)*(1-r)*(1-r)*(1-r)*( 4*r + 1);//(4r+1)(1-x)^4
					   A   =              -20 * (1-r)*(1-r)*(1-r);
					   B   = (r==0) ? 0 : +60 * (1-r)*(1-r) / r  ; 	              
	}
	g0 = A * d0       ;
	h0 = B * d0*d0 + A;
}


inline void rbfKernel2D( const double &d0, const double &d1, 
	                     double &val, 
						 double &g0 , double &g1 ,
						 double &h00, double &h01, double &h11)
{
	const int &m=RBFManager::getInst()->m_basisFuncMode;

	double r   = fabs( sqrt( d0*d0 + d1*d1) );
	double A, B;
	if(       m == 1){ val = r; //x
					   A   = (r != 0 ) ?  1/r       : 10000000; 
					   B   = (r != 0 ) ? -1/(r*r*r) : 10000000; 
	}else if( m == 2){ val = r*r*r;//x^3
					   A   = 3 * r;
					   B   = (r == 0) ? 0 : 3/r;
	}else if( m == 3){ val = r*r*r*r*r ; //x^5 
					   A   =  5 * r*r*r;
					   B   = 15 * r    ;
	}else if( m == 4){ val = (r==0) ? 0 : r*r*log(r)    ;//x^2 log(x)
					   A   = (r==0) ? 0 :   2*log(r) + 1;
					   B   = (r==0) ? 2 :   2 / (r*r)   ;

	}else if( m == 5){ val = (r==0) ? 0 : r*r*r*r*log(r);//x^4 log(x)
					   A   = (r==0) ? 0 : r*r*( 4*log(r) + 1 );
					   B   = (r==0) ? 0 : 8*log(r) + 6; 
					   
	}else {            val = (r>1)  ? 0 : (1-r)*(1-r)*(1-r)*(1-r)*( 4*r + 1);//(4r+1)(1-x)^4
					   A   =              -20 * (1-r)*(1-r)*(1-r);
					   B   = (r==0) ? 0 : +60 * (1-r)*(1-r) / r  ; 	              
	}

	g0  = A * d0;
	g1  = A * d1;
	h00 = B * d0*d0 + A;  h01 = B * d0*d1    ;
	                      h11 = B * d1*d1 + A;
}


inline void rbfKernel3D( const double &d0, const double &d1, const double &d2, 
	                     double &val, 
						 double &g0 , double &g1 , double &g2 ,
						 double &h00, double &h01, double &h02,
						              double &h11, double &h12,
									               double &h22)
{
	const int &m=RBFManager::getInst()->m_basisFuncMode;

	double r   = fabs( sqrt(d0*d0 + d1*d1 + d2*d2) );
	double A, B;
	if(       m == 1){ val = r; //x
					   A   = (r != 0 ) ?  1/r       : 10000000; 
					   B   = (r != 0 ) ? -1/(r*r*r) : 10000000; 
	}else if( m == 2){ val = r*r*r;//x^3
					   A   = 3 * r;
					   B   = (r == 0) ? 0 : 3/r;
	}else if( m == 3){ val = r*r*r*r*r ; //x^5 
					   A   =  5 * r*r*r;
					   B   = 15 * r    ;
	}else if( m == 4){ val = (r==0) ? 0 : r*r*log(r)    ;//x^2 log(x)
					   A   = (r==0) ? 0 :   2*log(r) + 1;
					   B   = (r==0) ? 2 :   2 / (r*r)   ;

	}else if( m == 5){ val = (r==0) ? 0 : r*r*r*r*log(r);//x^4 log(x)
					   A   = (r==0) ? 0 : r*r*( 4*log(r) + 1 );
					   B   = (r==0) ? 0 : 8*log(r) + 6; 
					   
	}else {            val = (r>1)  ? 0 : (1-r)*(1-r)*(1-r)*(1-r)*( 4*r + 1);//(4r+1)(1-x)^4
					   A   =              -20 * (1-r)*(1-r)*(1-r);
					   B   = (r==0) ? 0 : +60 * (1-r)*(1-r) / r  ; 	              
	}

	g0  = A * d0;
	g1  = A * d1;
	g2  = A * d2;
	h00 = B * d0*d0 + A;  h01 = B * d0*d1    ;  h02 = B * d0*d2    ;
	                      h11 = B * d1*d1 + A;  h12 = B * d1*d2    ;
						                        h22 = B * d2*d2 + A;

}