#include "StdAfx.h"
#include "RBFManager.h"
#include "tsparsematrix.h"
#include "TCore.h"

double RBFManager::m_approxCoef    = 1;
int    RBFManager::m_polynomMode   = 1;
int    RBFManager::m_basisFuncMode = 1;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     constraints                
//                pos   pos+∇c   pos+∇s  pos+∇j
//        color   △                --      --
// space  scale   △                △      --
//        joint   △                △      △
//RBF in COLOR space////////////////////////////////////////////////////////////////////////////////////////////////////
void RBFManager::RBFs_color( const TOGL2DImage &img, const vector<RBF_PosCP > &posCPs, double cPitch, double sPitch, float *scholarField)
{
	fprintf( stderr, "RBFs_color\n");
	const int W  = img.m_width, H  = img.m_height        ;
	const int Nc = (int)posCPs.size()  ;
	const int Np = (m_polynomMode < 0 ) ? 0: 
		           (m_polynomMode ==0 ) ? 1: //f(x) = Σφ(|x-ci|^2) + 0次(1) 
		           (m_polynomMode ==1 ) ? 2: //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) 
		                                  3; //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) + 2次( dimN + dimN-1 +...+1)
	const int N = Nc + Np;
	double *w = new double[ N ]; memset( w, 0, sizeof( double ) * N);
	double *b = new double[ N ]; memset( b, 0, sizeof( double ) * N);

	//compute RBF matrix////////////////////////////////////////////////////////////////// 
	TDenseMatrix  dMat(N,N);

	double val, g0, h0;
	//値制約
	for( int i=0; i<Nc; ++i)
	{
		const double *pi = posCPs[i].p;
		for( int j=i; j<Nc; ++j)
		{
			const double *pj = posCPs[j].p; 
			rbfKernel1D( pi[2]-pj[2], val, g0, h0);

			if( i==j) dMat[i][j]              = val + m_approxCoef;
			else      dMat[i][j] = dMat[j][i] = val               ;
		}
		if( m_polynomMode >= 0 ) { dMat[i][Nc+0 ] = dMat[ Nc+0 ][i] = 1              ; }
		if( m_polynomMode >= 1 ) { dMat[i][Nc+1 ] = dMat[ Nc+1 ][i] = pi[ 2 ]        ; }
		if( m_polynomMode >= 2 ) { dMat[i][Nc+2 ] = dMat[ Nc+2 ][i] = pi[ 2 ]*pi[ 2 ]; }
		b[i] = posCPs[i].val_c;
	}

	dMat.solveLinearSystem_LU( b, w );

	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{	
		double p = cPitch * img.m_RGBA[ (x+y*W)*4+0 ] / 255.0; 	

		double v=0;
		for( int i=0; i<Nc; ++i){
			rbfKernel1D( p-posCPs[i].p[2], val, g0, h0);
			v += w[i] * val;
		}
		
		if( m_polynomMode >= 0 ){ v += w[Nc+0 ]      ; }
		if( m_polynomMode >= 1 ){ v += w[Nc+1 ] * p  ; }
		if( m_polynomMode >= 2 ){ v += w[Nc+2 ] * p*p; }

		scholarField[ x+ y*W ] = (float) v;
	}
	delete[] b;
	delete[] w;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//RBF in SCALE space////////////////////////////////////////////////////////////////////////////////////////////////////
void RBFManager::RBFs_scale( const TOGL2DImage &img, const vector<RBF_PosCP > &posCPs,double cPitch, double sPitch,  float *scholarField)
{
	fprintf( stderr, "RBFs_scale\n");
	const int W  = img.m_width, H  = img.m_height        ;
	const int Nc = (int)posCPs.size()  ;
	const int Np = (m_polynomMode <  0 ) ? 0: 
		           (m_polynomMode == 0 ) ? 1:     //f(x) = Σφ(|x-ci|^2) + 0次(1) 
		           (m_polynomMode == 1 ) ? 1+2:   //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) 
		                               1+2+3; //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) + 2次( dimN + dimN-1 +...+1)
	const int N = Nc + Np;

	double *w = new double[ N ]; memset( w, 0, sizeof( double ) * N);
	double *b = new double[ N ]; memset( b, 0, sizeof( double ) * N);

	//compute RBF matrix////////////////////////////////////////////////////////////////// 
	TDenseMatrix  dMat(N,N);

	double g0, g1, h00, h01, h11, val;
	//値制約
	for( int i=0; i<Nc; ++i)
	{
		const double *pi = posCPs[i].p;
		for( int j=i; j<Nc; ++j)
		{
			const double *pj = posCPs[j].p;
			rbfKernel2D( pi[0]-pj[0], pi[1]-pj[1], val, g0,g1,h00,h01,h11);

			if( i==j) dMat[i][j]              = val + m_approxCoef;
			else      dMat[i][j] = dMat[j][i] = val ;
		}
		if( m_polynomMode >= 0 ) {                                   dMat[i][Nc+0   ] = dMat[ Nc+0   ][i] = 1          ; }
		if( m_polynomMode >= 1 ) { for( int k=0    ; k<2; ++k      ) dMat[i][Nc+1+k ] = dMat[ Nc+1+k ][i] = pi[k]      ; }
		if( m_polynomMode >= 2 ) { for( int k=0,n=0; n<2; ++n      )
								   for( int m=n    ; m<2; ++m, ++k ) dMat[i][Nc+3+k ] = dMat[ Nc+3+k ][i] = pi[n]*pi[m]; }
		b[i] = posCPs[i].val_s;
	}

	if(isCtrKeyOn()) dMat.solveLinearSystem_umfpack( b, w );
	else             dMat.solveLinearSystem_LU( b, w );

	double p[2];
	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{	
		p[0] = (x+0.5) * sPitch;
		p[1] = (y+0.5) * sPitch;

		double v=0;
		for( int i=0; i<Nc; ++i) {
			rbfKernel2D( p[0]-posCPs[i].p[0], 
				         p[1]-posCPs[i].p[1], val, g0,g1,h00,h01,h11);
			v += w[i] * val;
		}
		
		if( m_polynomMode >= 0 ){                                   v += w[Nc+0   ]            ; }
		if( m_polynomMode >= 1 ){ for( int k=0    ; k<2; ++k      ) v += w[Nc+1+k ] * p[k]     ; }
		if( m_polynomMode >= 2 ){ for( int k=0,n=0; n<2; ++n      )
								for( int m=n    ; m<2; ++m, ++k ) v += w[Nc+3+k ] * p[n]*p[m]; }

		scholarField[ x+ y*W ] = (float) v;
	}

	delete[] b;
	delete[] w;
}


//RBF in JOINT space////////////////////////////////////////////////////////////////////////////////////////////////////
void RBFManager::RBFs_joint(const TOGL2DImage &img, const vector<RBF_PosCP > &posCPs, double cPitch, double sPitch, float *scholarField, TTriangleMesh *zeroMesh)
{
	fprintf( stderr, "RBFs_joint\n");
	const int W  = img.m_width         ;
	const int H  = img.m_height        ;
	const int Nc = (int)posCPs.size()  ;
	const int Np = (m_polynomMode< 0 ) ? 0: 
		           (m_polynomMode==0 ) ? 1:     //f(x) = Σφ(|x-ci|^2) + 0次(1) 
		           (m_polynomMode==1 ) ? 1+3:   //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) 
		                               1+3+6;   //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) + 2次( dimN + dimN-1 +...+1)
	const int N = Nc + Np;

	double *w = new double[ N ]; memset( w, 0, sizeof( double ) * N);
	double *b = new double[ N ]; memset( b, 0, sizeof( double ) * N);

	//compute RBF matrix////////////////////////////////////////////////////////////////// 
	TDenseMatEx  dMat(N);

	double val, g0,g1,g2, h00,h01,h02,h11,h12,h22;
	//値制約
	for( int i=0; i<Nc; ++i)
	{
		const double *pi = posCPs[i].p;
		for( int j=i; j<Nc; ++j)
		{
			rbfKernel3D( pi[0]-posCPs[j].p[0],
				         pi[1]-posCPs[j].p[1],
						 pi[2]-posCPs[j].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);
			if( i==j) dMat[i][j]              = val + m_approxCoef;
			else      dMat[i][j] = dMat[j][i] = val;
		}
		if( m_polynomMode >= 0 ) {                                   dMat[i][Nc+0   ] = dMat[ Nc+0   ][i] = 1          ; }
		if( m_polynomMode >= 1 ) { for( int k=0    ; k<3; ++k      ) dMat[i][Nc+1+k ] = dMat[ Nc+1+k ][i] = pi[k]      ; }
		if( m_polynomMode >= 2 ) { for( int k=0,n=0; n<3; ++n      )
								   for( int m=n    ; m<3; ++m, ++k ) dMat[i][Nc+4+k ] = dMat[ Nc+4+k ][i] = pi[n]*pi[m]; }
		b[i] = posCPs[i].val_j;
	}
	dMat.Trace();
	dMat.LU_FactAndSolve( 0, b, w, 0);
	for( int i=0; i<N; ++i) fprintf( stderr, "%d %f    ", i, w[i]);

	double p[3];
	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{	
		p[0] = sPitch * (x+0.5);
		p[1] = sPitch * (y+0.5);
		p[2] = cPitch * img.m_RGBA[ (x+ y*W ) * 4] / 255.0;

		double v=0;
		for( int i=0; i<Nc; ++i) {
			rbfKernel3D( p[0]-posCPs[i].p[0],
				         p[1]-posCPs[i].p[1],
						 p[2]-posCPs[i].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);
			v += w[i] * val;
		}
		
		if( m_polynomMode >= 0 ){                                   v += w[Nc+0   ]            ; }
		if( m_polynomMode >= 1 ){ for( int k=0    ; k<3; ++k      ) v += w[Nc+1+k ] * p[k]     ; }
		if( m_polynomMode >= 2 ){ for( int k=0,n=0; n<3; ++n      )
								  for( int m=n    ; m<3; ++m, ++k ) v += w[Nc+4+k ] * p[n]*p[m]; }

		scholarField[ x+ y*W ] = (float) v;
	}



	
	if( zeroMesh != 0){
		//zero mesh生成
		const int    D = W;
		const double zPitch = 1 * cPitch / D;
		float *volume = new float[W*H*D];
		for( int z=0; z<D; ++z)
		for( int y=0; y<H; ++y)
		for( int x=0; x<W; ++x)
		{	
			double p[3];
			p[0] = (x+0.5) * sPitch;
			p[1] = (y+0.5) * sPitch;
			p[2] = (z+0.5) * zPitch;

			double v=0;

			for( int i=0; i<Nc; ++i){
				rbfKernel3D( p[0]-posCPs[i].p[0],
					         p[1]-posCPs[i].p[1],
							 p[2]-posCPs[i].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);	
				v += w[ i] * val;
			}

			if( m_polynomMode >= 0 ){                                   v += w[Nc+0  ]            ;}
			if( m_polynomMode >= 1 ){ for( int k=0    ; k<3; ++k      ) v += w[Nc+1+k] * p[k]     ;}
			if( m_polynomMode >= 2 ){ for( int k=0,n=0; n<3; ++n      )
									  for( int m=n    ; m<3; ++m, ++k ) v += w[Nc+4+k] * p[n]*p[m];}

			volume[ x+ y*W + z*W*H] = (float) v;
		}

		t_MartchingCubes_floatVolume( W,H,D, volume, 0, TVector3(W*sPitch, H*sPitch,  cPitch), *zeroMesh );
		delete[] volume;
	}





	delete[] b;
	delete[] w;
}


void  RBFManager::RBFs_color_graC( const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField)
{
	fprintf( stderr, "RBFs_scale_graS\n");
	const int W  = img.m_width         ;
	const int H  = img.m_height        ;
	const int Nc = (int)posCPs.size()  ;
	const int Ng = (int)graCPs.size()  ;
	const int Np = (m_polynomMode==0) ? 1 :    //f(x) = Σφ(|x-ci|^2) + 0次(1) 
		           (m_polynomMode==1) ? 2 :    //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) 
		           (m_polynomMode==2) ? 3 : 0; //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) + 2次( dimN + dimN-1 +...+1)

	double *w = new double[ Nc + Np      ]; memset( w, 0, sizeof( double ) * (Nc + Np      ) );
	double *b = new double[ Nc + Np + Ng ]; memset( b, 0, sizeof( double ) * (Nc + Np +  Ng) );

	TDenseMatrix  dMat(Nc + Np + Ng, Nc + Np);

	double val, g0, h0;
	//値制約
	for( int i=0; i<Nc; ++i)
	{
		const double *pi = posCPs[i].p;
		for( int j=i; j<Nc; ++j)
		{
			const double *pj = posCPs[j].p; 
			rbfKernel1D( pi[2]-pj[2], val, g0, h0);
			if( i==j) dMat[i][j]              = val + m_approxCoef;
			else      dMat[i][j] = dMat[j][i] = val;
		}
		if( m_polynomMode >= 0 ) { dMat[i][Nc+0 ] = dMat[ Nc+0 ][i] = 1              ; }
		if( m_polynomMode >= 1 ) { dMat[i][Nc+1 ] = dMat[ Nc+1 ][i] = pi[ 2 ]        ; }
		if( m_polynomMode >= 2 ) { dMat[i][Nc+2 ] = dMat[ Nc+2 ][i] = pi[ 2 ]*pi[ 2 ]; }
		b[i] = posCPs[i].val_c;
	}

	//grad制約
	for( int i=0,s=(int)graCPs.size(); i<s; ++i)
	{
		const double *pi = graCPs[i].p;
		const int    idx = Nc + Np + i;

		for( int j=0; j<Nc; ++j)
		{
			rbfKernel1D( pi[2]-posCPs[j].p[2], val, g0, h0);
			dMat[idx][j] = g0;
		}
		if( m_polynomMode >= 1 ) { dMat[ idx ][ Nc+1 ] = 1      ; }
		if( m_polynomMode >= 2 ) { dMat[ idx ][ Nc+2 ] = 2*pi[2]; }
		b[ idx  ] = graCPs[i].dir[2] > 0 ? 1 : -1;
	}

	dMat.solveLinearSystem_LU( b, w ); 

	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{	
		double c = img.m_RGBA[ 4*( x + y*W) ] * cPitch / 255.0;

		double v=0;

		for( int i=0; i<Nc; ++i){
			rbfKernel1D( c-posCPs[i].p[2], val, g0,h0);
			v += w[i] * val;
		}

		if( m_polynomMode >= 0 ){ v += w[Nc+0]        ;}
		if( m_polynomMode >= 1 ){ v += w[Nc+1] * c    ;}
		if( m_polynomMode >= 2 ){ v += w[Nc+2] * c * c;}

		scholarField[ x+ y*W ] = (float) v;
	}
	delete[] b;
	delete[] w;
}


//RBF in SCALE space with gradoConst////////////////////////////////////////////////////////////////////////////////////////////////////
void RBFManager::RBFs_scale_graS( const TOGL2DImage &img, const vector<RBF_PosCP > &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField)
{	
	fprintf( stderr, "RBFs_scale_graS\n");
	const int W  = img.m_width         ;
	const int H  = img.m_height        ;
	const int Nc = (int)posCPs.size()  ;
	const int Ng = (int)graCPs.size()  ;
	const int Np = (m_polynomMode==0) ? 1     :    //f(x) = Σφ(|x-ci|^2) + 0次(1) 
		           (m_polynomMode==1) ? 1+2   :    //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) 
		           (m_polynomMode==2) ? 1+2+3 : 0; //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) + 2次( dimN + dimN-1 +...+1)

	double *w = new double[ Nc + Np        ]; memset( w, 0, sizeof( double ) * (Nc + Np        ) );
	double *b = new double[ Nc + Np + 2*Ng ]; memset( b, 0, sizeof( double ) * (Nc + Np + 2* Ng) );

	TDenseMatrix  dMat(Nc + Np + 2*Ng, Nc + Np);
	double val, g0,g1, h00,h01,h11;
	//値制約///
	for( int i=0; i<Nc; ++i)
	{
		const double *pi = posCPs[i].p;
		for( int j=i; j<Nc; ++j)
		{
			rbfKernel2D( pi[0]-posCPs[j].p[0], 
				         pi[1]-posCPs[j].p[1], val,g0,g1, h00, h01, h11);

			if( i==j) dMat[i][j]              = val+ m_approxCoef;
			else      dMat[i][j] = dMat[j][i] = val ;
		}
		if( m_polynomMode >= 0 ) {                                   dMat[i][Nc+0   ] = dMat[Nc+0   ][i] = 1          ;}
		if( m_polynomMode >= 1 ) { for( int k=0    ; k<2; ++k      ) dMat[i][Nc+1+k ] = dMat[Nc+1+k ][i] = pi[k]      ;}
		if( m_polynomMode >= 2 ) { for( int k=0,n=0; n<2; ++n      )
								 for( int m=n    ; m<2; ++m, ++k ) dMat[i][Nc+3+k ] = dMat[Nc+3+k ][i] = pi[n]*pi[m];}
		b[i] = posCPs[i].val_s;
	}

	//grad制約
	for( int i=0,s=(int)graCPs.size(); i<s; ++i)
	{
		const double *pi = graCPs[i].p;
		const int    idx = Nc + Np + i*2;

		for( int j=0; j<Nc; ++j)
		{
			rbfKernel2D( pi[0]-posCPs[j].p[0], 
				         pi[1]-posCPs[j].p[1], val, g0, g1, h00,h01,h11);
			dMat[idx+0][j] = g0;
			dMat[idx+1][j] = g1;
		}
		if( m_polynomMode >= 1 ) { dMat[idx+0][ Nc+1 ] = 1;
			                       dMat[idx+1][ Nc+2 ] = 1;}
		if( m_polynomMode >= 2 ) { dMat[idx+0][ Nc+3 ] = 2*pi[0];   dMat[idx+1][ Nc+4] =   pi[0];
				                   dMat[idx+0][ Nc+4 ] =   pi[1];   dMat[idx+1][ Nc+5] = 2*pi[1]; }
		double l = sqrt( graCPs[i].dir[0] * graCPs[i].dir[0] + graCPs[i].dir[1] * graCPs[i].dir[1]);
		b[idx  ] = graCPs[i].dir[0] / l;
		b[idx+1] = graCPs[i].dir[1] / l;
	}

	dMat.solveLinearSystem_LU( b, w ); 

	double p[2];
	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{	
		p[0] = (x+0.5) * sPitch;
		p[1] = (y+0.5) * sPitch;

		double v=0;

		for( int i=0; i<Nc; ++i)
		{
			rbfKernel2D( p[0]-posCPs[i].p[0],
				         p[1]-posCPs[i].p[1], val, g0,g1,h00,h01,h11);
			v += w[i] * val;
		}
		if( m_polynomMode >= 0 ){                                 v += w[Nc+0  ]       ;}
		if( m_polynomMode >= 1 ){ for( int k=0    ; k<2; ++k    ) v += w[Nc+1+k] * p[k];}
		if( m_polynomMode >= 2 ){ for( int k=0,n=0; n<2; ++n    )
								for( int m=n    ; m<2; ++m,++k) v += w[Nc+3+k] * p[n]*p[m];}

		scholarField[ x+ y*W ] = (float) v;
	}
	delete[] b;
	delete[] w;
}


void RBFManager::RBFs_joint_graJ( const TOGL2DImage &img, const vector<RBF_PosCP > &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField)
{
	fprintf( stderr, "RBFs_joint_graJ\n");
	const int W  = img.m_width         ;
	const int H  = img.m_height        ;
	const int Nc = (int)posCPs.size()  ;
	const int Ng = (int)graCPs.size()  ;
	const int Np = (m_polynomMode==0) ? 1    :   //f(x) = Σφ(|x-ci|^2) + 0次(1) 
		           (m_polynomMode==1) ? 1+3  :   //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) 
		           (m_polynomMode==2) ? 1+3+6:0; //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) + 2次( dimN + dimN-1 +...+1)

	double *w = new double[ Nc + Np        ]; memset( w, 0, sizeof( double ) * (Nc + Np        ) );
	double *b = new double[ Nc + Np + 3*Ng ]; memset( b, 0, sizeof( double ) * (Nc + Np + 3* Ng) );

	TDenseMatrix  dMat(Nc + Np + 3*Ng, Nc + Np);

	double val, g0,g1,g2,  h00,h01,h02, h11,h12, h22;
	//値制約///
	for( int i=0; i<Nc; ++i)
	{
		const double *pi = posCPs[i].p;
		for( int j=i; j<Nc; ++j)
		{
			rbfKernel3D( pi[0]-posCPs[j].p[0],
				         pi[1]-posCPs[j].p[1],
						 pi[2]-posCPs[j].p[2], val, g0,g1,g2, h00,h01,h02, h11,h12, h22);

			if( i==j) dMat[i][j]              = val + m_approxCoef;
			else      dMat[i][j] = dMat[j][i] = val               ;
		}
		if( m_polynomMode >= 0 ) {                                   dMat[i][Nc+0   ] = dMat[Nc+0   ][i] = 1          ;}
		if( m_polynomMode >= 1 ) { for( int k=0    ; k<3; ++k      ) dMat[i][Nc+1+k ] = dMat[Nc+1+k ][i] = pi[k]      ;}
		if( m_polynomMode >= 2 ) { for( int k=0,n=0; n<3; ++n      )
								 for( int m=n    ; m<3; ++m, ++k ) dMat[i][Nc+4+k ] = dMat[Nc+4+k ][i] = pi[n]*pi[m];}
		b[i] = posCPs[i].val_j;
	}

	//grad制約
	for( int i=0,s=(int)graCPs.size(); i<s; ++i)
	{
		const double *pi = graCPs[i].p;
		const int    idx = Nc + Np + i*3;

		for( int j=0; j<Nc; ++j)
		{
			rbfKernel3D( pi[0]-posCPs[j].p[0], 
				         pi[1]-posCPs[j].p[1],
						 pi[2]-posCPs[j].p[2], val, g0, g1, g2,  h00,h01,h02, h11,h12, h22);
			dMat[idx+0][j] = g0;
			dMat[idx+1][j] = g1;
			dMat[idx+2][j] = g2;
		}
		//a0 + a1 *x + a2 *y + a3 *i + a4 *xx + a5* xy + a6 *xi + a7 *yy +a8 * yi +a9* ii
		if( m_polynomMode >= 1 ) { dMat[idx][ Nc+1 ] = 1;       dMat[idx+1][ Nc+2] = 1 ;      dMat[idx+2][ Nc+3 ] = 1;       }
		if( m_polynomMode >= 2 ) { dMat[idx][ Nc+4 ] = 2*pi[0]; dMat[idx+1][ Nc+5] =   pi[0]; dMat[idx+2][ Nc+6 ] =   pi[0];
				                   dMat[idx][ Nc+5 ] =   pi[1]; dMat[idx+1][ Nc+7] = 2*pi[1]; dMat[idx+2][ Nc+8 ] =   pi[1];
		                           dMat[idx][ Nc+6 ] =   pi[2]; dMat[idx+1][ Nc+8] =   pi[2]; dMat[idx+2][ Nc+9 ] =2* pi[2];}
		b[idx  ] = graCPs[i].dir[0];
		b[idx+1] = graCPs[i].dir[1];
		b[idx+2] = graCPs[i].dir[2];
	}

	dMat.solveLinearSystem_QR( b, w);

	double p[3];
	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{	
		p[0] = sPitch * (x+0.5);
		p[1] = sPitch * (y+0.5);
		p[2] = cPitch * img.m_RGBA[ 4*(x + y *W) ] / 255.0;

		double v=0;
		for( int i=0; i<Nc; ++i)
		{
			rbfKernel3D( p[0]-posCPs[i].p[0], 
				         p[1]-posCPs[i].p[1],
						 p[2]-posCPs[i].p[2], val, g0, g1, g2,  h00,h01,h02, h11,h12, h22);
			v += w[i] * val;
		}
		if( m_polynomMode >= 0 ){                                 v += w[Nc+0  ]       ;}
		if( m_polynomMode >= 1 ){ for( int k=0    ; k<3; ++k    ) v += w[Nc+1+k] * p[k];}
		if( m_polynomMode >= 2 ){ for( int k=0,n=0; n<3; ++n    )
								  for( int m=n    ; m<3; ++m,++k) v += w[Nc+4+k] * p[n]*p[m];}

		scholarField[ x+ y*W ] = (float) v;
	}

	delete[] b;
	delete[] w;
}



//todo 



//!!HRBF!!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!!HRBF!!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!!HRBF!!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!!HRBF!!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     constraints                
//                pos   pos+∇c   pos+∇s  pos+∇j
//        color   --               --      --
// space  scale   --                 △      --
//        joint   --                 △      △

//!!HRBF!!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!!HRBF!!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//f(x)/////////////////////////////////////////
//           A x  = b
//
//     Φ 　 -∇φ  Px       αi    c(値)
//     ∇φ　 - Hφ ∇P      βj =  n(勾配)
//      Pt     ∇Pt  0       pi     0
//scale spaceでHRBF計算 ▽s
void RBFManager::HRBFs_scale_graS( const TOGL2DImage &img, const vector<RBF_PosCP > &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField, TTriangleMesh *zeroMesh)
{
	fprintf( stderr, "HRBFs_scale_graS\n");
	const int W  = img.m_width         ;
	const int H  = img.m_height        ;

	const int Nc = (int) posCPs.size() ;
	const int Ng = (int) graCPs.size();
	const int Np = (m_polynomMode==0 ) ? 1         :   //f(x) = Σφ(|x-ci|^2) + 0次(1) 
		           (m_polynomMode==1 ) ? 1 + 2     :   //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) 
		           (m_polynomMode==2 ) ? 1 + 2 + 3 :0; //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) + 2次( dimN + dimN-1 +...+1)
	
	const int N    = Nc + Ng * 2 + Np;
	const int Pidx = Nc + Ng * 2;

	TDenseMatrix  dMat(N, N);
	double *w = new double[ N ]; memset( w, 0, sizeof( double ) * ( N ) );
	double *b = new double[ N ]; memset( b, 0, sizeof( double ) * ( N ) );

	fprintf( stderr, "aaaaaaaa%d %d %d\n", N, Nc, Ng);

	double g0,g1, h00,h01,h11, val;
	//値制約 f(xi ) = ci
	for( int i=0; i<Nc; ++i)
	{
		const double *pi = posCPs[i].p;
		
		for( int j=0; j<Nc; ++j)
		{
			rbfKernel2D( pi[0]-posCPs[j].p[0], 
				         pi[1]-posCPs[j].p[1], val, g0,g1,h00,h01, h11);
			dMat[i][j] = (i==j)? (val + m_approxCoef): val;
		}
		for( int j=0; j<Ng; ++j)
		{
			rbfKernel2D( pi[0]-graCPs[j].p[0], 
				         pi[1]-graCPs[j].p[1], val, g0,g1,h00,h01, h11);
			dMat[i][Nc + j*2    ] = -g0;
			dMat[i][Nc + j*2 + 1] = -g1;
		}

		//P(x) : f(xi ) = ci
		if( m_polynomMode >= 0 ) {                                   dMat[i][Pidx+0   ] = dMat[Pidx+0   ][i] = 1          ;}
		if( m_polynomMode >= 1 ) { for( int k=0    ; k<2; ++k      ) dMat[i][Pidx+1+k ] = dMat[Pidx+1+k ][i] = pi[k]      ;}
		if( m_polynomMode >= 2 ) { for( int k=0,n=0; n<2; ++n      )
								   for( int m=n    ; m<2; ++m, ++k ) dMat[i][Pidx+3+k ] = dMat[Pidx+3+k ][i] = pi[n]*pi[m];}
		b[ i ] = posCPs[i].val_s;
	}

	
	//勾配制約 ▽f(xi) = ni
	for( int i=0; i<Ng; ++i)
	{
		const double *pi = graCPs[i].p;
		const int    idx = Nc + 2*i   ;

		for( int j=0; j<Nc; ++j)
		{
			rbfKernel2D( pi[0]-posCPs[j].p[0], pi[1]-posCPs[j].p[1], val, g0,g1, h00, h01, h11);
			dMat[idx  ][ j ] = g0;   
			dMat[idx+1][ j ] = g1;  
		}

		for( int j=0; j<Ng; ++j)
		{
			rbfKernel2D( pi[0]-graCPs[j].p[0], 
				         pi[1]-graCPs[j].p[1], val, g0,g1, h00, h01, h11);
			dMat[idx  ][Nc+2*j  ] = -h00;   dMat[idx  ][Nc+2*j+1] = -h01; 
			dMat[idx+1][Nc+2*j  ] = -h01;   dMat[idx+1][Nc+2*j+1] = -h11;
		}

		//P = a0 + a1 x + a2 y + a3 xx + a4 xy + + a5 yy
		if( m_polynomMode >= 1 ) {  dMat[idx+0][ Pidx+1 ] = dMat[ Pidx+1][idx+0] = 1;         dMat[idx+1][ Pidx+2 ] = dMat[ Pidx+2][idx+1] = 1;}
		if( m_polynomMode >= 2 ) {  dMat[idx+0][ Pidx+ 3] = dMat[ Pidx+3][idx+0] = 2*pi[0];   dMat[idx+1][ Pidx+4 ] = dMat[ Pidx+4][idx+1] =   pi[0];
							        dMat[idx+0][ Pidx+ 4] = dMat[ Pidx+4][idx+0] =   pi[1];   dMat[idx+1][ Pidx+5 ] = dMat[ Pidx+5][idx+1] = 2*pi[1];}

		b[idx  ] = graCPs[i].dir[0];
		b[idx+1] = graCPs[i].dir[1];
	}


	dMat.solveLinearSystem_LU( b, w );



	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{	
		double p[2];
		p[0] = sPitch * (x+0.5) ;
		p[1] = sPitch * (y+0.5) ;

		double v=0;

		for( int i=0; i<Nc; ++i){
			rbfKernel2D( p[0]-posCPs[i].p[0], p[1]-posCPs[i].p[1], val, g0,g1, h00, h01, h11);
			v += w[ i] * val;
		}
		for( int i=0; i<Ng; ++i){
			rbfKernel2D( p[0]-graCPs[i].p[0], p[1]-graCPs[i].p[1], val, g0,g1, h00, h01, h11);
			v += - w[Nc+2*i  ] * g0 ;
			v += - w[Nc+2*i+1] * g1 ;
		}

		if( m_polynomMode >= 0 ){                                   v += w[Pidx+0  ]            ;}
		if( m_polynomMode >= 1 ){ for( int k=0    ; k<2; ++k      ) v += w[Pidx+1+k] * p[k]     ;}
		if( m_polynomMode >= 2 ){ for( int k=0,n=0; n<2; ++n      )
								  for( int m=n    ; m<2; ++m, ++k ) v += w[Pidx+3+k] * p[n]*p[m];}

		scholarField[ x+ y*W ] = (float) v;
	}



	if( zeroMesh != 0 ){
		//zero mesh生成
		const int D = W;
		float *volume = new float[W*H*D];
		for( int z=0; z<D; ++z)
		for( int y=0; y<H; ++y)
		for( int x=0; x<W; ++x)
		{	
			double p[2];
			p[0] = sPitch * (x+0.5) ;
			p[1] = sPitch * (y+0.5) ;

			double v=0;

			for( int i=0; i<Nc; ++i){
				rbfKernel2D( p[0]-posCPs[i].p[0], p[1]-posCPs[i].p[1], val, g0,g1, h00, h01, h11);
				v +=  w[i] * val;
			}
			for( int i=0; i<Ng; ++i)
			{
				rbfKernel2D( p[0]-graCPs[i].p[0], p[1]-graCPs[i].p[1], val, g0,g1, h00, h01, h11);
				v += - w[Nc+2*i  ] * g0 ;
				v += - w[Nc+2*i+1] * g1 ;
			}

			if( m_polynomMode >= 0 ){                                   v += w[Pidx+0  ]            ;}
			if( m_polynomMode >= 1 ){ for( int k=0    ; k<2; ++k      ) v += w[Pidx+1+k] * p[k]     ;}
			if( m_polynomMode >= 2 ){ for( int k=0,n=0; n<2; ++n      )
									  for( int m=n    ; m<2; ++m, ++k ) v += w[Pidx+3+k] * p[n]*p[m];}
			volume[ x+ y*W + z*W*H] = (float) v;
		}

		t_MartchingCubes_floatVolume( W,H,D, volume, 0, TVector3(W*sPitch, H*sPitch, cPitch), *zeroMesh );
		delete[] volume;
	}


	delete[] b;
	delete[] w;
}

void RBFManager::HRBFs_color_graC( const TOGL2DImage &img, const vector<RBF_PosCP > &posCPs, const vector<RBF_GradCP> &graCPs,double cPitch, double sPitch,  float *scholarField)
{
	fprintf( stderr, "HRBFs_color_graC\n");
	const int W  = img.m_width         ;
	const int H  = img.m_height        ;

	const int Nc = (int) posCPs.size() ;
	const int Ng = (int) graCPs.size();
	const int Np = (m_polynomMode==0 ) ? 1 :    //f(x) = Σφ(|x-ci|^2) + 0次(1) 
		           (m_polynomMode==1 ) ? 2 :    //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) 
		           (m_polynomMode==2 ) ? 3 : 0; //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) + 2次( dimN + dimN-1 +...+1)
	
	const int N    = Nc + Ng  + Np;
	const int Pidx = Nc + Ng ;

	TDenseMatrix  dMat(N, N);
	double *w = new double[ N ]; memset( w, 0, sizeof( double ) * ( N ) );
	double *b = new double[ N ]; memset( b, 0, sizeof( double ) * ( N ) );

	double g0, h0, val;
	//値制約 f(xi ) = ci
	for( int i=0; i<Nc; ++i)
	{		
		const double c = posCPs[i].p[2];
		for( int j=i; j<Nc; ++j)
		{
			const double *pj = posCPs[j].p; 
			rbfKernel1D( c-pj[2], val, g0,h0);

			if( i==j) dMat[i][j]              = val + m_approxCoef;
			else      dMat[i][j] = dMat[j][i] = val;
		}
		for( int j=0; j<Ng; ++j)
		{
			rbfKernel1D( c - graCPs[j].p[2], val, g0,h0);
			dMat[i][ Nc + j ] = -g0;
		}
		if( m_polynomMode >= 0 ) { dMat[i][Pidx+0 ] = dMat[ Pidx+0 ][i] = 1  ; }
		if( m_polynomMode >= 1 ) { dMat[i][Pidx+1 ] = dMat[ Pidx+1 ][i] = c  ; }
		if( m_polynomMode >= 2 ) { dMat[i][Pidx+2 ] = dMat[ Pidx+2 ][i] = c*c; }
		b[i] = posCPs[i].val_c;
	}

	//勾配制約 ▽f(xi) = ni
	for( int i=0; i<Ng; ++i)
	{
		const double c   = graCPs[i].p[2];
		const int    idx = Nc + i        ;
		for( int j=0; j<Nc; ++j)
		{
			rbfKernel1D(c-posCPs[j].p[2], val, g0, h0);
			dMat[idx  ][ j ] = g0;   
		}

		for( int j=0; j<Ng; ++j)
		{
			rbfKernel1D( c - graCPs[j].p[2], val, g0, h0);
			dMat[ idx ][ Nc+j ] = -h0; 
		}
		//P = a0 + a1 x + a2 y + a3 xx + a4 xy + + a5 yy
		if( m_polynomMode >= 1 ) {  dMat[idx][ Pidx+ 1] = dMat[ Pidx+1][idx] = 1  ; }
		if( m_polynomMode >= 2 ) {  dMat[idx][ Pidx+ 2] = dMat[ Pidx+2][idx] = 2*c; }  

		b[idx  ] = graCPs[i].dir[2]>0?1:-1;
	}
	//dMat.Trace();
	dMat.solveLinearSystem_LU( b, w );

	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{	
		double c =  cPitch * img.m_RGBA[ 4*(x + y*W) ] / 255.0;
		
		double v=0;
		for( int i=0; i<Nc; ++i){
			rbfKernel1D( c-posCPs[i].p[2], val, g0, h0);
			v += w[ i] * val;
		}

		for( int i=0; i<Ng; ++i)
		{
			rbfKernel1D( c-graCPs[i].p[2], val, g0, h0);
			v += - w[ Nc+i ] * g0 ;
		}

		if( m_polynomMode >= 0 ){  v += w[ Pidx  ]      ;}
		if( m_polynomMode >= 1 ){  v += w[ Pidx+1] * c  ;}
		if( m_polynomMode >= 2 ){  v += w[ Pidx+2] * c*c;}

		scholarField[ x+ y*W ] = (float) v;
	}

	//debug
	for( int cc=0;cc<255;++cc){
		double c = cPitch * cc;
		double v=0;
		for( int i=0; i<Nc; ++i){
			rbfKernel1D( c-posCPs[i].p[2], val, g0,h0);
			v += w[ i] * val;
		}

		for( int i=0; i<Ng; ++i)
		{
			rbfKernel1D( c-graCPs[i].p[2], val, g0,h0);
			v += - w[ Nc+i ] * g0 ;
		}

		if( m_polynomMode >= 0 ){  v += w[ Pidx  ]      ;}
		if( m_polynomMode >= 1 ){  v += w[ Pidx+1] * c  ;}
		if( m_polynomMode >= 2 ){  v += w[ Pidx+2] * c*c;}
		fprintf( stderr, "%f\n", v);
	}

	for( int i=0; i<Nc; ++i){ fprintf( stderr, "%d %f\n", i, posCPs[i].p[2]);}


	delete[] b;
	delete[] w;
}



//joint spaceでHRBF計算 gradientも3Dを考慮 (▽j)
void RBFManager::HRBFs_joint_graJ(  const TOGL2DImage &img, const vector<RBF_PosCP> &posCPs, const vector<RBF_GradCP> &graCPs, double cPitch, double sPitch, float *scholarField, 
	double colDirCoef, TOGL3DImage4 *volume4, TTriangleMesh *zeroMesh)
{
	fprintf( stderr, "HRBFs_joint_graJ\n");
	const int W  = img.m_width         ;
	const int H  = img.m_height        ;
	const int Nc = (int) posCPs.size() ;
	const int Ng = (int) graCPs.size();
	const int Np = (m_polynomMode==0 ) ? 1         :    //f(x) = Σφ(|x-ci|^2) + 0次(1) 
		           (m_polynomMode==1 ) ? 1 + 3     :    //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) 
		           (m_polynomMode==2 ) ? 1 + 3 + 6 : 0; //f(x) = Σφ(|x-ci|^2) + 0次(1) + 1次(dimN) + 2次( dimN + dimN-1 +...+1)
	const int N    = Nc + Ng * 3 + Np;
	const int Pidx = Nc + Ng * 3;

	TDenseMatEx dMat( N );
	double *w = new double[ N ]; memset( w, 0, sizeof( double ) * ( N ) );
	double *b = new double[ N ]; memset( b, 0, sizeof( double ) * ( N ) );

	//値制約 f(xi ) = ci
	double val, g0, g1, g2, h00,h01,h02,h11,h12,h22;
	for( int i=0; i<Nc; ++i)
	{		
		const double *pi = posCPs[i].p;
		for( int j=0; j<Nc; ++j){
			rbfKernel3D( pi[0]-posCPs[j].p[0],pi[1]-posCPs[j].p[1],pi[2]-posCPs[j].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);
			dMat[i][j] = (i==j)? (val + m_approxCoef): val;
		}

		for( int j=0; j<Ng; ++j)
		{
			rbfKernel3D( pi[0]-graCPs[j].p[0], pi[1]-graCPs[j].p[1],pi[2]-graCPs[j].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);
			dMat[i][Nc + j*3    ] = -g0;
			dMat[i][Nc + j*3 + 1] = -g1;
			dMat[i][Nc + j*3 + 2] = -g2;
		}
		//P(x) : f(xi ) = ci
		if( m_polynomMode >= 0 ) {                                   dMat[i][Pidx+0   ] = dMat[Pidx+0   ][i] = 1          ;}
		if( m_polynomMode >= 1 ) { for( int k=0    ; k<3; ++k      ) dMat[i][Pidx+1+k ] = dMat[Pidx+1+k ][i] = pi[k]      ;}
		if( m_polynomMode >= 2 ) { for( int k=0,n=0; n<3; ++n      )
								   for( int m=n    ; m<3; ++m, ++k ) dMat[i][Pidx+4+k ] = dMat[Pidx+4+k ][i] = pi[n]*pi[m];}
		b[ i ] = posCPs[i].val_j;
	}

	//勾配制約 ▽f(xi) = ni
	for( int i=0; i<Ng; ++i)
	{
		const double *pi = graCPs[i].p;
		const int    idx = Nc + 3*i   ;

		for( int j=0; j<Nc; ++j){
			rbfKernel3D( pi[0]-posCPs[j].p[0], pi[1]-posCPs[j].p[1], pi[2]-posCPs[j].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);
			dMat[idx  ][ j ] = g0;   
			dMat[idx+1][ j ] = g1;  
			dMat[idx+2][ j ] = g2;  
		}

		for( int j=0; j<Ng; ++j){
			rbfKernel3D( pi[0]-graCPs[j].p[0], pi[1]-graCPs[j].p[1], pi[2]-graCPs[j].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);
			dMat[idx  ][Nc+3*j] = -h00;   dMat[idx  ][Nc+3*j+1] = -h01;   dMat[idx  ][Nc+3*j+2] = -h02; 
			dMat[idx+1][Nc+3*j] = -h01;   dMat[idx+1][Nc+3*j+1] = -h11;   dMat[idx+1][Nc+3*j+2] = -h12;
			dMat[idx+2][Nc+3*j] = -h02;   dMat[idx+2][Nc+3*j+1] = -h12;   dMat[idx+2][Nc+3*j+2] = -h22;

			if( i == j ) {
				dMat[idx  ][Nc+3*j  ] += m_approxCoef*10;   
				dMat[idx+1][Nc+3*j+1] += m_approxCoef*10;
				dMat[idx+2][Nc+3*j+2] += m_approxCoef*10;
			}
		}

		//P = a0 + a1 x + a2 y + a3 i + a4 xx + + a5 xy + a6 xi + a7 yy + a8 yi + a9 ii
		if( m_polynomMode >= 1 ) {dMat[idx  ][ Pidx+1 ] = dMat[ Pidx+1][idx  ] = 1;         
		                          dMat[idx+1][ Pidx+2 ] = dMat[ Pidx+2][idx+1] = 1;
		                          dMat[idx+2][ Pidx+3 ] = dMat[ Pidx+3][idx+2] = 1; }
		if( m_polynomMode >= 2 ) {dMat[idx  ][ Pidx+4 ] = dMat[ Pidx+4][idx  ] = 2*pi[0];   
							      dMat[idx  ][ Pidx+5 ] = dMat[ Pidx+5][idx  ] =   pi[1];   
								  dMat[idx  ][ Pidx+6 ] = dMat[ Pidx+6][idx  ] =   pi[2];   
								  dMat[idx+1][ Pidx+5 ] = dMat[ Pidx+5][idx+1] =   pi[0];
								  dMat[idx+1][ Pidx+7 ] = dMat[ Pidx+7][idx+1] = 2*pi[1];
								  dMat[idx+1][ Pidx+8 ] = dMat[ Pidx+8][idx+1] =   pi[2];
								  dMat[idx+2][ Pidx+6 ] = dMat[ Pidx+6][idx+2] =   pi[0];
								  dMat[idx+2][ Pidx+8 ] = dMat[ Pidx+8][idx+2] =   pi[1];
								  dMat[idx+2][ Pidx+9 ] = dMat[ Pidx+9][idx+2] = 2*pi[2]; }
		b[idx  ] = graCPs[i].dir[0];
		b[idx+1] = graCPs[i].dir[1];
		b[idx+2] = graCPs[i].dir[2];
	}

	dMat.LU_factorization_Full(0);
	dMat.LU_SolveLinearSystem( b,w);

	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{	
		double p[3];
		p[0] = (x+0.5) * sPitch;
		p[1] = (y+0.5) * sPitch;
		p[2] =  cPitch * img.m_RGBA[ (x+W*y)*4 ] /255.0;

		double v=0;

		for( int i=0; i<Nc; ++i){
			rbfKernel3D(p[0]-posCPs[i].p[0],  p[1]-posCPs[i].p[1],  p[2]-posCPs[i].p[2], val,  g0,g1,g2, h00,h01,h02,h11,h12,h22);
			v += w[ i] * val;
		}
		for( int i=0; i<Ng; ++i){
			rbfKernel3D( p[0]-graCPs[i].p[0],  p[1]-graCPs[i].p[1],  p[2]-graCPs[i].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);	
			v += - w[Nc+3*i  ] * g0 ;
			v += - w[Nc+3*i+1] * g1 ;
			v += - w[Nc+3*i+2] * g2 ;
		}

		if( m_polynomMode >= 0 ){                                   v += w[Pidx+0  ]            ;}
		if( m_polynomMode >= 1 ){ for( int k=0    ; k<3; ++k      ) v += w[Pidx+1+k] * p[k]     ;}
		if( m_polynomMode >= 2 ){ for( int k=0,n=0; n<3; ++n      )
								  for( int m=n    ; m<3; ++m, ++k ) v += w[Pidx+4+k] * p[n]*p[m];}

		scholarField[ x+ y*W ] = (float) v;
	}


	if( W*H < 20000 && zeroMesh != 0 && volume4 != 0)
	{
		//zero mesh生成
		const int    D = W;
		volume4->allocateImage( W, H, D, &TCore::getInst()->m_ogl);
		const double zPitch = 1 * cPitch / D;
		float *field3D = new float[W*H*D];
		for( int z=0; z<D; ++z)
		for( int y=0; y<H; ++y)
		for( int x=0; x<W; ++x)
		{	
			double p[3];
			p[0] = (x+0.5) * sPitch;
			p[1] = (y+0.5) * sPitch;
			p[2] = (z+0.5) * zPitch;

			double v=0;

			for( int i=0; i<Nc; ++i){
				rbfKernel3D( p[0]-posCPs[i].p[0], p[1]-posCPs[i].p[1], p[2]-posCPs[i].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);	
				v += w[ i] * val;
			}
			for( int i=0; i<Ng; ++i)
			{
				rbfKernel3D( p[0]-graCPs[i].p[0], p[1]-graCPs[i].p[1], p[2]-graCPs[i].p[2], val, g0,g1,g2, h00,h01,h02,h11,h12,h22);	
				v += - w[Nc+3*i  ] * g0 ;
				v += - w[Nc+3*i+1] * g1 ;
				v += - w[Nc+3*i+2] * g2 ;
			}

			if( m_polynomMode >= 0 ){                                   v += w[Pidx+0  ]            ;}
			if( m_polynomMode >= 1 ){ for( int k=0    ; k<3; ++k      ) v += w[Pidx+1+k] * p[k]     ;}
			if( m_polynomMode >= 2 ){ for( int k=0,n=0; n<3; ++n      )
									  for( int m=n    ; m<3; ++m, ++k ) v += w[Pidx+4+k] * p[n]*p[m];}

			field3D[ x+ y*W + z*W*H] = (float) v;

			byte rr,gg,bb; heuColor( 8*v, rr, gg, bb);
			volume4->m_RGBA[ 4 * (x+ y*W + z*W*H) + 0] = rr;
			volume4->m_RGBA[ 4 * (x+ y*W + z*W*H) + 1] = gg;
			volume4->m_RGBA[ 4 * (x+ y*W + z*W*H) + 2] = bb;
		}

		t_MartchingCubes_floatVolume( W,H,D, field3D, 0, TVector3(W*sPitch, H*sPitch,  cPitch), *zeroMesh );
		delete[] field3D;
	}
	delete[] b;
	delete[] w;
}
