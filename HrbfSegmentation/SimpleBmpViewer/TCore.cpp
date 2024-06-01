#include "StdAfx.h"
#include "TCore.h"
#include "tsparsematrix.h"
#include "Tutil.h"

#include "SimpleBmpViewer.h"


#ifndef TGCUT_STATIC_LINK
#define TGCUT_STATIC_LINK
#endif 

#ifndef FOREBACK_MAX 
#define FOREBACK_MAX 1000000000.0
#endif 

#include "./TGCutDll2010/TGcutStdHead.h" //include after defining the "TGCUT_STATIC_LINK"必要あり
#include <list>
#include "TRBF3D.h"


/*---------------------------------------------------
画像は[0,1]×[0,1]に配置される
pixel pitchは 1 / max(W, H)

画素値は m_colorPitch * [0,1]となる

m_colorCoef = 0.6としておく
これはpitchというよりはcolorのscaleingであることに注意

0.6は大きいように見えるが、対象領域のサイズとも比例する値になってしまっているためやむをえない

図5では、ほぼ全体を覆うものを抽出しているが
[0,1][0,1]のものを抽出するさいには　0.6

[0,0.5][0,0.5]のものを抽出するさいには　0.3程度で同等の効果が得られ
[0,0.1][0,0.1]のものを抽出するさいには　0.06程度で同等の効果が得られることに注意


---------------------------------------------------*/

//



static double maxImageSize = 1.0  ;


//byte *sobelX, byte sobelY should be allocated before
//color [0,255]-->[0,1]
static void sobel( const byte *rgba, const int W, const int H, double *sobelX, double *sobelY)
{
	double filterX[3][3]={	{-1,  0,  1}, {-2,  0,  2}, {-1,  0,  1}};
	double filterY[3][3]={  {-1, -2, -1}, { 0,  0,  0}, { 1,  2,  1}};

	for(int y = 0, idx=0; y < H ; y++)
	for(int x = 0       ; x < W ; x++, ++idx)
	{
		double vx = 0, vy = 0 ;		
		for( int wy = -1; wy < 2; ++wy)
		for( int wx = -1; wx < 2; ++wx)
		{
			int pixId = 4*( idx + wx + wy * W );
			double pixVal;
			if( x + wx < 0 || x + wx > W-1 || y + wy < 0 || y + wy > H-1 ) pixVal = 0;
			else pixVal	= (rgba[ pixId + 0] + rgba[ pixId + 1] + rgba[ pixId + 2] ) * 0.33333 / 255.0;
			vx += filterX[wy+1][wx+1] * pixVal;
			vy += filterY[wy+1][wx+1] * pixVal;
		}
		sobelX[ idx ] = fabs( vx ) ;
		sobelY[ idx ] = fabs( vy ) ;
	}
}

TCore::~TCore(void){
}

TCore::TCore(void)
{
	m_bAutoDrawMode = false;
	m_theta    = 0    ;
	m_phai     = 10   ;
	m_colPitch = 0.6;

	//RBF parametes////////////////////////////////
	RBFManager::getInst()->m_basisFuncMode = 2           ;
	RBFManager::getInst()->m_polynomMode   = 0           ;
	RBFManager::getInst()->m_approxCoef    = 0           ;

	//solver and constraint mode///////////////////
	m_solverMode =  HRBF_JOINT_GraJ ;
	m_bPosConstF = false; m_bPosConstB = false; m_bPosConstM = true;
	m_bGraConstF = false; m_bGraConstB = false; m_bGraConstM = true;


	m_ogl.SetClearColor( 0.1, 0.1, 0.1, 0.5);
	m_ogl.setInitEyePosition( TVector3(0,0,5), TVector3(0,1,0), TVector3( 0,0,0));
	m_ogl.setOrthoModeViewSize(  2.0);

	//load image/////////////////////////////////////////////////////////////////////////
	CString         filter("bmp Files (*.bmp;*.bmp)|*.bmp; *.bmp||");
	CFileDialog     selDlg(TRUE, NULL, NULL, OFN_HIDEREADONLY, filter);
	bool bInverted;
	if (selDlg.DoModal() == IDOK) m_imgOriginal.allocateFromFile( selDlg.GetPathName(), bInverted, &m_ogl );
	else exit(1);
	const int W = m_imgOriginal.m_width ;
	const int H = m_imgOriginal.m_height;

	if( W > H ) m_pixPitch = maxImageSize / W;
	else        m_pixPitch = maxImageSize / H;
	m_imgSize.Set( W*m_pixPitch, H*m_pixPitch);
	
	m_imgOriginal.flipImageInY();                  m_imgOriginal.m_DoInterpolation = true ;
	m_imgSmooth.allocateImage( m_imgOriginal, 0 ); m_imgSmooth  .m_DoInterpolation = false;
	m_imgVis   .allocateImage( m_imgOriginal, 0 ); m_imgVis     .m_DoInterpolation = true ;
	//for( int i=0; i< 4; ++i) m_imgOrig  .gaussianFilter33();
	for( int i=0; i< 5; ++i) m_imgSmooth.gaussianFilter33();


	m_dlg.Create( IDD_DIALOG_PARAM );
	m_dlg.ShowWindow( SW_SHOW );

	m_UVs.resize( W*H);
	vector<TVector3>  Vs( W*H);
	vector<TTriangle> Ps; Ps.reserve( W*H*2);

	for( int y = 0; y<H-1; ++y)
	for( int x = 0; x<W-1; ++x){
		int idx = x + y*W;
		Ps.push_back( TTriangle( idx, idx+1   ,idx+1+W) );
		Ps.push_back( TTriangle( idx, idx+1+W ,idx+  W) );
	}

	m_imgManif.initFromVertsPolys( Vs, Ps);
	m_mesh     = m_imgManif;
	for( int y = 0; y<H; ++y)
	for( int x = 0; x<W; ++x){
		double xp = (x+0.5)*m_pixPitch;
		double yp = (y+0.5)*m_pixPitch;
		int idx = x + y*W;
		m_UVs[idx].Set( (x+0.5)/(W-1.0), (y+0.5)/(H-1.0));
		m_imgManif.m_verts[idx].Set( xp,yp, m_imgSmooth.m_RGBA[ idx*4 ] / 255.0 * m_colPitch);
	}
	m_imgManif.updateNormal();
	m_HRBF_colorDirCoef = 0.5;
}




///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
static float ambiZ[4] = {0,.7f,.2f,.5f};
static float diffZ[4] = {0,1,0.2f ,.8f};
static float shinZ[1] = {1 };
static float ImgManifAmbi[4] = {.5f, .5f, .5f,.8f};
static float ImgManifDiff[4] = {.8f, .8f, .8f,.8f};


static const void setOglMatrix( const double rx, const double ry, const double phai, const double theta)
{
	glTranslated( 0,  ry*0.5, 0);
	glRotated   (-theta, 1,0,0); 
	glTranslated( ry*0.5,0,0);
	glRotated   ( phai , 0,0,1); 
	glTranslated( -ry*0.5,0,0);
	glTranslated( 0, -ry*0.5, 0);
}

static void drawImage(  TOGL2DImage &img, double sizeX, double sizeY)
{
	img.bind(0);
	glDisable( GL_LIGHTING );
	glEnable( GL_TEXTURE_2D);
	glColor3d(1,1,1);
	glBegin( GL_QUADS );
		glTexCoord2d(0,0); glVertex3d( 0,0,0);
		glTexCoord2d(1,0); glVertex3d( sizeX,0,0);
		glTexCoord2d(1,1); glVertex3d( sizeX,sizeY,0);
		glTexCoord2d(0,1); glVertex3d( 0,sizeY,0);
	glEnd();
	glDisable( GL_TEXTURE_2D);
}

static void drawLines(const vector<pair<TVector2, TVector2>> &points)
{
	glTranslated( 0,0,0.1);

	glDisable( GL_LIGHTING );
	glLineWidth( 5 );

	glColor3d( 1,1,0);
	glBegin( GL_LINES );
	for( int i=0, s=(int)points.size(); i<s; ++i){ glVertex2dv( points[i].first .data ); glVertex2dv( points[i].second.data );}
	glEnd();

	glTranslated( 0,0,0.1);
	
	if( false) {
		glLineWidth( 1 );
		glColor3d( 0,0,1);

		glBegin( GL_LINES );
		for( int i=0, s=(int)points.size(); i<s; ++i){ glVertex2dv( points[i].first .data ); glVertex2dv( points[i].second.data ); }
		glEnd();
	}
	glTranslated( 0,0, -0.1);
	glTranslated( 0,0, -0.1);
}


//userの指定したhandleを描画
static void draw_HandleIn2D(vector< TBoundConst > CPs_bound, vector< TInOutConst > CPs_inOut )
{
	glDisable( GL_LIGHTING );
	glLineWidth(4);
	glPointSize( 12 );
	glColor3d(0,0,1);

	glTranslated( 0,0,1);
	glBegin( GL_LINES );
	for( int i=0; i<(int) CPs_bound.size(); ++i)
	{
		glColor3d( 0,0,1); glVertex2dv( CPs_bound[i].m_pos.data );
		glColor3d( 1,0,1); glVertex2d( CPs_bound[i].m_pos[0] + 0.05 * CPs_bound[i].m_dir[0], 
			                           CPs_bound[i].m_pos[1] + 0.05 * CPs_bound[i].m_dir[1]);
	}
	glEnd();

	glBegin( GL_POINTS);
		glColor3d( 1,1,0); for( int i=0; i<(int)CPs_bound.size(); ++i)  glVertex2dv( CPs_bound[i].m_pos.data );
		glColor3d( 1,0,1); for( int i=0; i<(int)CPs_inOut.size(); ++i)  glVertex2dv( CPs_inOut[i].m_pos.data );
	glEnd();
	glTranslated( 0,0,-1);
}


static void drawRect( double rx, double ry)
{
	glEnable( GL_BLEND );
	glDisable( GL_LIGHTING );
	glColor4d(0,0.9,0.9,0.3);
	glBegin( GL_TRIANGLES );
		glVertex3d( 0 , 0,-0.01); glVertex3d( rx, 0,-0.01); glVertex3d( rx,ry,-0.01); 
		glVertex3d( 0 , 0,-0.01); glVertex3d( rx,ry,-0.01); glVertex3d( 0 ,ry,-0.01); 
	glEnd();
	glDisable( GL_BLEND );
}


static void drawSurface( const TTriangleMesh &mesh, bool bTransparent, float diff[4], float ambi[4] )
{
	float shin [1] = {128};
	float spec [4] = {1,1,1,1};
	glEnable(GL_LIGHTING );
	glMaterialfv( GL_FRONT, GL_SHININESS, shin );
	glMaterialfv( GL_FRONT, GL_SPECULAR,  spec );
	glMaterialfv( GL_FRONT, GL_DIFFUSE,   diff );
	glMaterialfv( GL_FRONT, GL_AMBIENT,   ambi );
	if( bTransparent ) glEnable( GL_BLEND );
		mesh.drawPolygons();
	if( bTransparent ) glDisable( GL_BLEND );
}
static void drawSurface_texture( const TTriangleMesh &mesh, TOGL2DImage &img, vector<TVector2> &UVs, bool bTransparent)
{
	float shin [1] = {128};
	float spec [4] = {1,1,1,1};
	float ambi [4] = {.7f,.7f,.7f,1};
	float diff [4] = {1,1,1,1};
	glEnable(GL_LIGHTING );
	img.bind(0);

	glMaterialfv( GL_FRONT, GL_SHININESS, shin );
	glMaterialfv( GL_FRONT, GL_SPECULAR,  spec );
	glMaterialfv( GL_FRONT, GL_DIFFUSE,   diff );
	glMaterialfv( GL_FRONT, GL_AMBIENT,   ambi );
	if( bTransparent ) glEnable( GL_BLEND );
	glEnable( GL_TEXTURE_2D );
		glBegin(GL_TRIANGLES );
		for( int i = 0; i < mesh.getPnum(); ++i){
			int *p = mesh.m_polys[i].idx;	
			glNormal3dv( mesh.m_v_norm[ p[0]].data ); glTexCoord2dv( UVs[p[0]].data ); glVertex3dv( mesh.m_verts[ p[0] ].data );
			glNormal3dv( mesh.m_v_norm[ p[1]].data ); glTexCoord2dv( UVs[p[1]].data ); glVertex3dv( mesh.m_verts[ p[1] ].data );
			glNormal3dv( mesh.m_v_norm[ p[2]].data ); glTexCoord2dv( UVs[p[2]].data ); glVertex3dv( mesh.m_verts[ p[2] ].data );
		}
		glEnd();

	glDisable( GL_TEXTURE_2D );
	if( bTransparent ) glDisable( GL_BLEND );
}

static void drawSurface_texture3D( const TTriangleMesh &mesh, TOGL3DImage4 &img, double colorPitch, bool bTransparent)
{
	float shin [1] = {128};
	float spec [4] = {1,1,1,1};
	float diff [4] = {1,1,1,1};
	float ambi [4] = {.7f,.7f,.7f,1};

	glMaterialfv( GL_FRONT, GL_SHININESS, shin );
	glMaterialfv( GL_FRONT, GL_SPECULAR,  spec );
	glMaterialfv( GL_FRONT, GL_DIFFUSE,   diff );
	glMaterialfv( GL_FRONT, GL_AMBIENT,   ambi );
	glEnable(GL_LIGHTING );
	img.bind(0);

	if( bTransparent ) glEnable( GL_BLEND );
	glEnable( GL_TEXTURE_3D );
		glBegin(GL_TRIANGLES );
		for( int i = 0; i < mesh.getPnum(); ++i)
		{
			int *p = mesh.m_polys[i].idx;
			const TVector3 &p1 = mesh.m_verts[ p[0] ];
			const TVector3 &p2 = mesh.m_verts[ p[1] ];
			const TVector3 &p3 = mesh.m_verts[ p[2] ];
			TVector3 t1 = p1; t1.data[2] /= colorPitch;
			TVector3 t2 = p2; t2.data[2] /= colorPitch;
			TVector3 t3 = p3; t3.data[2] /= colorPitch;

			glNormal3dv( mesh.m_v_norm[ p[0]].data ); glTexCoord3dv( t1.data ); glVertex3dv( p1.data );
			glNormal3dv( mesh.m_v_norm[ p[1]].data ); glTexCoord3dv( t2.data ); glVertex3dv( p2.data );
			glNormal3dv( mesh.m_v_norm[ p[2]].data ); glTexCoord3dv( t3.data ); glVertex3dv( p3.data );
		}
		glEnd();
	glDisable( GL_TEXTURE_3D );
	if( bTransparent ) glDisable( GL_BLEND );
}

static void draw_CPs(const vector<RBF_PosCP > &posCPs, 
	                 const vector<RBF_GradCP> &graCPs)
{
	float aR[4] = {.8f,.1f,.1f,1}, dR[4] = {.8f,.1f,.1f,1};
	float aG[4] = {.1f,.8f,.1f,1}, dG[4] = {.1f,.8f,.1f,1};
	float aB[4] = {.1f,.1f,.8f,1}, dB[4] = {.1f,.1f,.8f,1};

	glEnable( GL_LIGHTING );
	for( int i=0; i<posCPs.size(); ++i){
		if     ( posCPs[i].val_j < 0) { glMaterialfv( GL_FRONT, GL_DIFFUSE, dB );glMaterialfv( GL_FRONT, GL_AMBIENT, aB );} 
		else if( posCPs[i].val_j > 0) { glMaterialfv( GL_FRONT, GL_DIFFUSE, dR );glMaterialfv( GL_FRONT, GL_AMBIENT, aR );} 
		else                          { glMaterialfv( GL_FRONT, GL_DIFFUSE, dG );glMaterialfv( GL_FRONT, GL_AMBIENT, aG );} 
		t_drawSphere_norm( 0.01, TVector3( posCPs[i].p[0], posCPs[i].p[1], posCPs[i].p[2] ) );
	}

	glColor3d(1,0,0);
	glDisable( GL_LIGHTING );
	glLineWidth( 7 );
	glBegin( GL_LINES );
	for( int i=0; i<graCPs.size(); ++i)
	{
		const double *p = graCPs[i].p;
		const double *d = graCPs[i].dir;
		glVertex3d( p[0]         , p[1]         , p[2]          );
		glVertex3d( p[0]+0.2*d[0], p[1]+0.2*d[1], p[2]+0.2*d[2] );
	}
	glEnd();
}



void TCore::drawScene()
{
	fprintf( stderr, "draw scene!!\n");
	glEnable( GL_LIGHT1);
	glEnable( GL_LIGHT0);
	glDisable( GL_CULL_FACE );
	const double rx = m_imgSize[0];
	const double ry = m_imgSize[1];

	//draw image original image----------------------------------------------------------------
	glPushMatrix();
		drawLines ( m_boundLines          );
		drawImage ( m_imgOriginal, rx, ry );
		if( !isCtrKeyOn() ) draw_HandleIn2D( m_CPs_bound, m_CPs_inOut); 
	glPopMatrix();

	//draw distance field in 2D----------------------------------------------------------------
	glPushMatrix();
		glTranslated( -rx*2.0,0,0);
		drawImage( m_imgVis, rx, ry );
	glPopMatrix();

	//draw Image Manifold ----------------------------------------------------------------
	glPushMatrix();
		glTranslated( rx * 2.0, 0, 0);
		setOglMatrix( rx, ry, m_phai, m_theta );
		drawRect(rx,ry);
		drawSurface( m_imgManif, false, ImgManifDiff, ImgManifAmbi);			
		//drawSurface_texture3D( m_visPlane, m_imgVolume, false );
		if( isCtrKeyOn() ) drawSurface( m_zeroMesh, true, diffZ, ambiZ);			
		if( isAltKeyOn() || m_bAutoDrawMode ) draw_CPs( m_posCPs, m_gradCPs );
	glPopMatrix();

	//draw Image Manifold with scalarfield vis------------------------------------------------
	glPushMatrix();
		glTranslated( rx * 4.0, 0, 0);
		setOglMatrix( rx, ry, m_phai, m_theta );
		drawRect(rx,ry);
		if( isSpaceKeyOn() ) drawSurface_texture( m_imgManif, m_imgVis, m_UVs, false );
		else                 drawSurface_texture3D( m_imgManif, m_imgVolume, m_colPitch, false );

		if( isCtrKeyOn() ) drawSurface( m_zeroMesh, true, diffZ, ambiZ);	
		if( isAltKeyOn() ) draw_CPs( m_posCPs, m_gradCPs );

	glPopMatrix();

	//draw Image Manifold with scalarfield vis------------------------------------------------
	glPushMatrix();
		glTranslated( rx * 6.0, 0, 0);
		setOglMatrix( rx, ry, m_phai, m_theta );
		drawRect(rx,ry);
		glEnable( GL_LIGHTING );	
		drawSurface( m_mesh, false, diffZ, ambiZ );
		drawLines( m_boundLines );
	glPopMatrix();	

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TCore::addBoundConst( TVector2 pos,  TVector2 dir)
{
	const int W = m_imgOriginal.m_width ;
	const int H = m_imgOriginal.m_height;

	if( ! t_isInWindow2D( TVector3(0,0,0), TVector3(m_imgSize[0], m_imgSize[1],0) , TVector3( pos[0], pos[1],0), 0) ) return;
	if( dir.Length() < m_pixPitch * 2 ) return;
	dir.Normalize_Self();
	
	//ADD constraint points///////////////////////////////////////////////////////
	m_CPs_bound.push_back( TBoundConst( pos , dir ) );

	updateSegmentation();
}


void TCore::addInOutConst( vector< TVector2 > &ps,  bool tf     )
{	
	const int W = m_imgOriginal.m_width ;
	const int H = m_imgOriginal.m_height;

	for( int i=0; i<(int)ps.size(); ++i){
		if( ! t_isInWindow2D( TVector3(0,0,0), TVector3(m_imgSize[0], m_imgSize[1],0), TVector3( ps[i][0], ps[i][1],0), 0) ) continue;
		m_CPs_inOut.push_back( TInOutConst( ps[i], tf )) ;
	}
	updateSegmentation();
}


void TCore::OnKeyDown( char nChar )
{
	if( nChar == 49 /*1*/) {RBFManager::m_basisFuncMode = 1; updateSegmentation();m_ogl.RedrawWindow();}
	if( nChar == 50 /*2*/) {RBFManager::m_basisFuncMode = 2; updateSegmentation();m_ogl.RedrawWindow();}
	if( nChar == 51 /*3*/) {RBFManager::m_basisFuncMode = 3; updateSegmentation();m_ogl.RedrawWindow();}
	if( nChar == 52 /*4*/) {RBFManager::m_basisFuncMode = 4; updateSegmentation();m_ogl.RedrawWindow();}
	if( nChar == 53 /*5*/) {RBFManager::m_basisFuncMode = 5; updateSegmentation();m_ogl.RedrawWindow();}
	if( nChar == 54 /*6*/) {RBFManager::m_basisFuncMode = 6; updateSegmentation();m_ogl.RedrawWindow();}
	if( nChar == 55 /*7*/) 
	{
		/*
		m_bAutoDrawMode = true;
		int N          = 100;
		double thetaStep = 1.5;
		for( int i=0; i<N; ++i)
		{
			m_theta-= thetaStep;
			TCore::getInst()->m_ogl.RedrawWindow();
			Sleep( 100 );
		}
		for( int i=0; i<N; ++i)
		{
			m_theta += thetaStep;
			TCore::getInst()->m_ogl.RedrawWindow();
			Sleep( 100 );
		}
		m_bAutoDrawMode = false;
		*/
	}
	if( nChar == 56 /*8*/) {
		/*
		m_bAutoDrawMode = true;
		vector<TVector3 > Vs(4);
		vector<TTriangle> Ps(2);
		Ps[0].Set(0,1,3);
		Ps[1].Set(0,3,2);
		Vs[0].Set(0,1,0.7); Vs[1].Set(1,1,0.7);
		Vs[2].Set(0,1,0); Vs[3].Set(1,1,0);
		m_visPlane.initFromVertsPolys(Vs, Ps );
		int    N = 100;
		double step = -1.0 / N;
		for( int i=0; i<N; ++i)
		{
			m_visPlane.m_verts[0].data[1] += step;
			m_visPlane.m_verts[1].data[1] += step;
			m_visPlane.m_verts[2].data[1] += step;
			m_visPlane.m_verts[3].data[1] += step;
			TCore::getInst()->m_ogl.RedrawWindow();
			Sleep( 100 );
		}
		m_visPlane.Clear();	
		m_bAutoDrawMode = false;
		*/
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



static void t_regionGrowing(const int W, const int H, const vector<TPixelInfo> &_seeds, const float *voxVal, byte *voxFlg, int postGrowthNum)
{
	//region growing///////////////////////////////////////////////////////////////////////////////
	memset( voxFlg, 0, sizeof( byte  ) * W*H ); // 0:non visited, 1 fronteer, 2 fixed as fore, 3 fixed as back
	
	list< TPixelInfo > fronteer;
	for( int i=0; i<(int)_seeds.size(); ++i) fronteer.push_back( _seeds[i] );
	TPixelInfo pivPix(0,0,0);

	while( !fronteer.empty() )
	{
		pivPix = fronteer.back(); fronteer.pop_back();
		const int &i = pivPix.idx, &x = pivPix.x, &y = pivPix.y;
			
		//+の領域のgrow
		if( 0 < x   && voxFlg[i- 1] == 0){
			if( voxVal[i-1] < 0){ voxFlg[i -1] = 3;                                                  }
			else                { voxFlg[i -1] = 2; fronteer.push_back( TPixelInfo( x-1, y, i-1 ) ); }
		}
		if( x <W-1  && voxFlg[i+ 1] == 0){
			if( voxVal[i+1] < 0){ voxFlg[i +1] = 3;                                                  }
			else                { voxFlg[i +1] = 2; fronteer.push_back( TPixelInfo( x+1, y, i+1 ) ); }
		}
		if( 0 < y   && voxFlg[i- W] == 0){
			if( voxVal[i-W] < 0){ voxFlg[i -W] = 3;                                                    }
			else                { voxFlg[i -W] = 2; fronteer.push_back( TPixelInfo( x  , y-1, i-W ) ); }
		}
		if( y < H-1 && voxFlg[i+ W] == 0){
			if( voxVal[i+W] < 0){ voxFlg[i +W] = 3;                                                    }
			else                { voxFlg[i +W] = 2; fronteer.push_back( TPixelInfo( x  , y+1, i+W ) ); }
		}
	}
	for( int i=0; i<W*H; ++i) if( voxFlg[i] == 0 )voxFlg[i] =3;  

	//boundaryをpostGrowthNum分だけ広げる
	for( int kk=0; kk<postGrowthNum; ++kk)
	{
		for( int y = 1; y < H-1; ++y)
		for( int x = 1; x < W-1; ++x){
			const int idx = x + W*y;
			if( voxFlg[idx] == 1 || voxFlg[idx] == 3) continue;
			if( voxFlg[idx+1] == 3 ) voxFlg[idx+1] = 1;
			if( voxFlg[idx-1] == 3 ) voxFlg[idx-1] = 1;
			if( voxFlg[idx+W] == 3 ) voxFlg[idx+W] = 1;
			if( voxFlg[idx-W] == 3 ) voxFlg[idx-W] = 1;
		}

		for( int i=0; i<W*H; ++i) if( voxFlg[i] == 1 )voxFlg[i] =2;  
	}
}



void TCore::updateSegmentation()
{
	fprintf( stderr, "update segmentation here !");
	const int    W    = m_imgSmooth.m_width ;
	const int    H    = m_imgSmooth.m_height;
	const byte*  rgba = m_imgSmooth.m_RGBA;

	vector< TPixelInfo > regionCandidate;

	m_posCPs .clear();
	m_gradCPs.clear();

	//位置/勾配制約作成////////////////////////////////////////////////////////////////
	const double colC    = m_colPitch ;//User parameter 1 color Pitch
	const double pixR    = m_pixPitch ;//User parameter 2 scale Pitch
	const double posOfst = 0.5          ;//User parameter 3 offset 
	double blendA        = 1.0          ;//User parameter 4 (なくてもいいかも)

	TVector2 gradMid, gradInt, gradExt;
	TVector3 midDir, intDir, extDir;

	for( int i=0; i<(int) m_CPs_bound.size(); ++i)
	{
		m_CPs_bound[i].m_dir.Normalize_Self();

		const TVector2 cpD  = m_CPs_bound[i].m_dir;
		const TVector2 cpP  = m_CPs_bound[i].m_pos + posOfst * pixR * cpD;
		const TVector2 pInt =   cpP                + posOfst * pixR * cpD;   
		const TVector2 pExt =   cpP                - posOfst * pixR * cpD;  
		fprintf( stderr, "\ncps %d\n", i);
		cpP .Trace();
		pInt.Trace();
		pExt.Trace();

		int  cpX = (int)( cpP[0]/pixR),    cpY = (int)(  cpP[1]/pixR),     cpIdx = ( cpX +  cpY*W) * 4;
		int intX = (int)(pInt[0]/pixR),   intY = (int)( pInt[1]/pixR),    intIdx = (intX + intY*W) * 4;
		int extX = (int)(pExt[0]/pixR),   extY = (int)( pExt[1]/pixR),    extIdx = (extX + extY*W) * 4;

		//image gradient のサンプリング (境界が谷付近の時，少しだけ前景に近いところでsamplingしたほうがいい結果が出る)
		m_imgSmooth.getGrad5(  cpX,  cpY, gradMid );  gradMid *= colC / pixR / 255.0 ;
		m_imgSmooth.getGrad5( intX, intY, gradInt );  gradInt *= colC / pixR / 255.0;  
		m_imgSmooth.getGrad5( extX, extY, gradExt );  gradExt *= colC / pixR / 255.0; 

		double len = gradMid.Length();
		if( gradMid*cpD < 0 ) midDir.Set( blendA * len* cpD[0] - (1-blendA) * gradMid[0], blendA * len* cpD[1] - (1-blendA) * gradMid[1], 0 );
		else                  midDir.Set( blendA * len* cpD[0] + (1-blendA) * gradMid[0], blendA * len* cpD[1] + (1-blendA) * gradMid[1], 0 );
		midDir.Normalize_Self();
		midDir[2] = m_HRBF_colorDirCoef * (midDir[0]*gradMid[0] + midDir[1]*gradMid[1]);
		midDir.Normalize_Self();
		
		//Image manifold r(x,y) =(x,y,I(x,y))上の点 r(xM,yM)における，接平面上の点は, a rx + b ry = a*(1,0,Ix) + b*(0,1,Iy)で与えられる		
		intDir.Set( midDir[0], midDir[1],  m_HRBF_colorDirCoef * (midDir[0]*gradInt[0] + midDir[1]*gradInt[1] )); 
		extDir.Set( midDir[0], midDir[1],  m_HRBF_colorDirCoef * (midDir[0]*gradExt[0] + midDir[1]*gradExt[1] ));
		intDir.Normalize_Self();
		extDir.Normalize_Self();

		//位置制約の挿入 (todo joint spaceなら距離もjoint spaceで計算する)
		if( m_bPosConstM ){
			m_posCPs.push_back( RBF_PosCP(  cpP, colC*rgba[ cpIdx]/255.0, colC*rgba[ cpIdx+1]/255.0, colC*rgba[ cpIdx+2]/255.0, 0.0 , 0.0, 0.0       ));
		}
		if( m_bPosConstF ){
			double ds = t_distance(pInt,cpP);
			double dc = fabs( colC *rgba[intIdx]/255.0 - colC *rgba[cpIdx]/255.0 ) ;
			double dj = sqrt( ds * ds + dc*dc);
			if( true ){
				TVector3 POS = TVector3(cpP[0], cpP[1], colC*rgba[ cpIdx]/255.0) + posOfst * pixR* midDir;
				m_posCPs.push_back( RBF_PosCP( TVector2(POS[0],POS[1]), POS[2], POS[2], POS[2], dc, posOfst* pixR, posOfst* pixR));
			}else {
				m_posCPs.push_back( RBF_PosCP( pInt, colC*rgba[intIdx]/255.0, colC*rgba[intIdx+1]/255.0, colC*rgba[intIdx+2]/255.0, dc, ds, dj));
			}
		}
		if( m_bPosConstB ){
			double ds = t_distance(pExt,cpP);
			double dc = fabs( colC *rgba[extIdx]/255.0 - colC *rgba[cpIdx]/255.0 ) ;
			double dj = sqrt( ds * ds + dc*dc);

			if( true ){
				TVector3 POS = TVector3(cpP[0], cpP[1], colC*rgba[ cpIdx]/255.0) - posOfst* pixR * midDir;
				m_posCPs.push_back( RBF_PosCP( TVector2(POS[0],POS[1]), POS[2], POS[2], POS[2], -dc, -posOfst* pixR, -posOfst* pixR));
			}else {
				m_posCPs.push_back( RBF_PosCP( pExt, colC*rgba[extIdx]/255.0, colC*rgba[extIdx+1]/255.0, colC*rgba[extIdx+2]/255.0,-dc,-ds,-dj ));
			}
		}
		if( m_bGraConstM ) m_gradCPs.push_back( RBF_GradCP( cpP, colC*rgba[ cpIdx]/255.0, colC*rgba[ cpIdx+1]/255.0, colC*rgba[ cpIdx+2]/255.0, midDir,  cpX,  cpY  ) );
		if( m_bGraConstF ) m_gradCPs.push_back( RBF_GradCP(pInt, colC*rgba[intIdx]/255.0, colC*rgba[intIdx+1]/255.0, colC*rgba[intIdx+2]/255.0, intDir, intX, intY  ) );
		if( m_bGraConstB ) m_gradCPs.push_back( RBF_GradCP(pExt, colC*rgba[extIdx]/255.0, colC*rgba[extIdx+1]/255.0, colC*rgba[extIdx+2]/255.0, extDir, extX, extY  ) );

		regionCandidate.push_back( TPixelInfo( intX, intY, intIdx/4 ) );
	}
	fprintf( stderr, "const num Pos%d Dir%d\n", m_posCPs.size(), m_gradCPs.size() );
	//位置のみの制約
	const int preCPsize = (int)m_posCPs.size();
	for( int i=0; i<m_CPs_inOut.size(); ++i)
	{
		TVector2 &cpP = m_CPs_inOut[i].m_pos;
		double   sign = m_CPs_inOut[i].m_tf ? 1:-1;
		int cpIdx = ((int)( cpP[0] / pixR) + (int)( cpP[1] / pixR)*W) * 4;

		double r = rgba[ cpIdx + 0] * colC;
		double g = rgba[ cpIdx + 1] * colC;
		double b = rgba[ cpIdx + 2] * colC;

		double minDist= DBL_MAX;
		int    nearIdx;
		for( int k=0; k<preCPsize; ++k){
			double *p = m_gradCPs[k].p;
			double d = (p[0]-cpP[0])*(p[0]-cpP[0]) + (p[1]-cpP[1])*(p[1]-cpP[1]) + (p[2]-r)*(p[2]-r);
			if( d < minDist ){
				minDist = d;
				nearIdx = k;
			}
		}

		minDist = (preCPsize == 0) ? 1 : sqrt( minDist );

		double dc = fabs( m_posCPs[nearIdx].p[2] - r);
		double ds = sqrt( (m_posCPs[nearIdx].p[0] - cpP[0])*(m_posCPs[nearIdx].p[0] - cpP[0]) + 
			              (m_posCPs[nearIdx].p[1] - cpP[1])*(m_posCPs[nearIdx].p[1] - cpP[1]) ) ;
		double dj = sqrt( dc*dc + ds*ds );
		m_posCPs.push_back( RBF_PosCP( cpP, r,g,b, sign*dc, sign*ds, sign*dj ) );
	}

	//compute RBF  Joint/Scale/Color/S+J/S+C//////////////////////////////////////////////////////////////
	float* field  = new float[W*H]; memset( field , 0, sizeof( float ) * W*H );
	byte*  voxFlg = new byte [W*H]; memset( voxFlg, 0, sizeof( byte  ) * W*H );

	if     ( m_solverMode == RBF_COLOR      ) RBFManager:: RBFs_color     ( m_imgSmooth, m_posCPs,            colC, pixR, field );
	else if( m_solverMode == RBF_SCALE      ) RBFManager:: RBFs_scale     ( m_imgSmooth, m_posCPs,            colC, pixR, field );
	else if( m_solverMode == RBF_JOINT      ) RBFManager:: RBFs_joint     ( m_imgSmooth, m_posCPs,            colC, pixR, field, &m_zeroMesh );
	else if( m_solverMode == RBF_COLOR_GraC ) RBFManager:: RBFs_color_graC( m_imgSmooth, m_posCPs, m_gradCPs, colC, pixR, field );
	else if( m_solverMode == RBF_SCALE_GraS ) RBFManager:: RBFs_scale_graS( m_imgSmooth, m_posCPs, m_gradCPs, colC, pixR, field );
	else if( m_solverMode == RBF_JOINT_GraJ ) RBFManager:: RBFs_joint_graJ( m_imgSmooth, m_posCPs, m_gradCPs, colC, pixR, field );
	else if( m_solverMode ==HRBF_COLOR_GraC ) RBFManager::HRBFs_color_graC( m_imgSmooth, m_posCPs, m_gradCPs, colC, pixR, field );
	else if( m_solverMode ==HRBF_SCALE_GraS ) RBFManager::HRBFs_scale_graS( m_imgSmooth, m_posCPs, m_gradCPs, colC, pixR, field, &m_zeroMesh  );
	else if( m_solverMode ==HRBF_JOINT_GraJ ) RBFManager::HRBFs_joint_graJ( m_imgSmooth, m_posCPs, m_gradCPs, colC, pixR, field, m_HRBF_colorDirCoef, &m_imgVolume, &m_zeroMesh );

	// region growingで内部の領域を取得
	const double r = m_pixPitch;
	m_imgVis.allocateImage( W, H, &m_ogl );


	//image 計算
	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x){
		int i= x + y*W;
		byte br,bg,bb; 
		heuColor( 8.0 * field[i], br,bg,bb ); 
		m_imgVis.setPix( x,y, br,bg, bb, 128 );
	}

	//boundary 計算
	m_boundLines.clear(); 
	if( true ){
		t_regionGrowing( W, H, regionCandidate, field, voxFlg, (int)( posOfst ) );
		for( int y=0; y<H-1; ++y)
		for( int x=0; x<W-1; ++x){
			int i = x + y * W;
			//   ↓ちょっとトリッキーですね 2 * 3ならfore * backってこと
			if( voxFlg[i] * voxFlg[i+1] == 6 ) m_boundLines.push_back( make_pair( TVector2((x+1)*r, y*r), TVector2((x+1)*r,(y+1)*r) ) );
			if( voxFlg[i] * voxFlg[i+W] == 6 ) m_boundLines.push_back( make_pair( TVector2( x*r,(y+1)*r), TVector2((x+1)*r,(y+1)*r) ) );
		}
	}else{
		for( int y=0; y<H-1; ++y)
		for( int x=0; x<W-1; ++x){
			int i = x + y * W;
			//   ↓ちょっとトリッキーですね 2 * 3ならfore * backってこと
			if( field[i] * field[i+1] <= 0 ) m_boundLines.push_back( make_pair( TVector2((x+1)*r, y*r), TVector2((x+1)*r,(y+1)*r) ) );
			if( field[i] * field[i+W] <= 0 ) m_boundLines.push_back( make_pair( TVector2( x*r,(y+1)*r), TVector2((x+1)*r,(y+1)*r) ) );
		}
	}

	double coef = 1.0;// / max( fieldMax[0], fieldMax[1] ) * 30.0;
	for( int y=0; y<H; ++y)
	for( int x=0; x<W; ++x)
	{
		int i= x + y*W;
		m_mesh.m_verts[i].Set( (x+0.5)*pixR, (y+0.5)*pixR, max(-5, min( 5, coef * field[i]) ));
	}
	m_mesh.updateNormal();

	delete[] voxFlg;
	delete[] field;

	vector<TVector3> Vs;
	vector<TTriangle> Ps;

	//meshのいらない部分取り除く[px,1-px] * [py,1-py] * [0.01, colPitch-0.01] のそとにあるものは取り除く

	for( int i=0; i<m_zeroMesh.getPnum(); ++i)
	{
		const TVector3 &x0 = m_zeroMesh.m_verts[ m_zeroMesh.m_polys[i].idx[0] ];
		const TVector3 &x1 = m_zeroMesh.m_verts[ m_zeroMesh.m_polys[i].idx[1] ];
		const TVector3 &x2 = m_zeroMesh.m_verts[ m_zeroMesh.m_polys[i].idx[2] ];
		if( m_pixPitch >= x0[0] || x0[0] >= m_imgSize[0] - m_pixPitch ) continue;  
		if( m_pixPitch >= x1[0] || x1[0] >= m_imgSize[0] - m_pixPitch ) continue;  
		if( m_pixPitch >= x2[0] || x2[0] >= m_imgSize[0] - m_pixPitch ) continue;
		if( m_pixPitch >= x0[1] || x0[1] >= m_imgSize[1] - m_pixPitch ) continue;  
		if( m_pixPitch >= x1[1] || x1[1] >= m_imgSize[1] - m_pixPitch ) continue;  
		if( m_pixPitch >= x2[1] || x2[1] >= m_imgSize[1] - m_pixPitch ) continue;  

		if( 0.01       >= x0[2] || x0[2] >= m_colPitch-0.01 ) continue;  
		if( 0.01       >= x1[2] || x1[2] >= m_colPitch-0.01 ) continue;  
		if( 0.01       >= x2[2] || x2[2] >= m_colPitch-0.01 ) continue; 
		Ps.push_back( m_zeroMesh.m_polys[i] );
	}
	for( int i=0; i<m_zeroMesh.getVnum(); ++i) Vs.push_back( m_zeroMesh.m_verts[i] );

	m_zeroMesh.initFromVertsPolys( Vs, Ps );
	fprintf( stderr, "finish update segmentation!!\n\n\n");

	m_ogl.RedrawWindow();
}

//vis作る ok
//surface作る
//全体の大きさとうをまとめる










