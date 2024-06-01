#ifndef __TOGL_H_INCLUDED__
#define __TOGL_H_INCLUDED__

/*-----------------------------------------------
TOGL for OpenGL on MFC 
Written by Takashi Ijiri @ riken / 2011/10/6

call OnCreate( ) 
     OnSize  ( )
	 OnDraw_Begin()
	 OnDraw_End()
	 OnDestroy()   from ****View.h
------------------------------------------------*/
#pragma warning (disable:4786)
#pragma warning (disable:4996)

#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")

#ifndef GLEW_STATIC
#define GLEW_STATIC
#endif

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cmath>
#include "tmath.h"



class TOGLEyeParam
{
public:
	TVector3 m_pos, m_focus, m_upDir;

	TOGLEyeParam(){}
	TOGLEyeParam( const TOGLEyeParam& src ){
		m_pos   = src.m_pos;
		m_focus = src.m_focus;
		m_upDir = src.m_upDir;
	}
	inline void set(  TVector3 _pos,  TVector3 _focus, TVector3 _upDir ){
		m_pos   = _pos; 
		m_focus = _focus; 
		m_upDir = _upDir;
	}	
	inline void SetPos  ( double x, double y, double z){ m_pos  .Set(x,y,z);}
	inline void SetFocus( double x, double y, double z){ m_focus.Set(x,y,z);}
	inline void SetUpDir( double x, double y, double z){ m_upDir.Set(x,y,z);}

	inline void callGluLookAt(){
		gluLookAt(  m_pos  .data[0]  , m_pos.data[1], m_pos  .data[2],
					m_focus.data[0], m_focus.data[1], m_focus.data[2],
					m_upDir.data[0], m_upDir.data[1], m_upDir.data[2]);
	}
};




class TOGL
{
	CWnd*  m_pWnd;//trgt MFC view 
	CDC*   m_pDC ;//devide context of trgtMFC view
	HGLRC  m_hRC ;//open gl handle

	bool   m_bIsDrawing ; 
	CPoint m_prevpos    ;
	int    m_cx, m_cy   ;

	enum {
		T_BUTTON_NON, 
		T_BUTTONDOWN_ZOOM,      T_BUTTONDOWN_ROTATE    , T_BUTTONDOWN_ROTYDIR, 
		T_BUTTONDOWN_TRANSLATE, T_BUTTONDOWN_ORTHO_ZOOM, T_BUTTONDOWN_ORTHO_TRANS,
	} m_ButtonState;

	double m_orthoViewSize;
	TOGLEyeParam m_eyeParam;
	TOGLEyeParam m_initEyeParam;
	float    m_clearColor[4];

public:
	inline HGLRC* getHGLRC(){ return &m_hRC;}
	/////////////////////////////////////////////////////////////////////////////////////////
	//CONSTRACTOR////////////////////////////////////////////////////////////////////////////
	TOGL():
		m_pWnd(0), m_pDC(0),
		m_ButtonState(T_BUTTON_NON), m_bIsDrawing(false),m_orthoViewSize(10)
	{
		m_clearColor[0] = m_clearColor[1] = m_clearColor[2] = 1.0f; m_clearColor[3] = 0.5f ;

		m_initEyeParam.SetPos  ( 0, 20, 50 );
		m_initEyeParam.SetFocus( 0, 20, 0  );
		m_initEyeParam.SetUpDir( 0,  1, 0  );
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	//CALLED BY MFC **VIEW CLASS/////////////////////////////////////////////////////////////
public:
	inline void setOrthoModeViewSize( double s){ m_orthoViewSize = s;}
	virtual BOOL OnCreate(CWnd* pWnd)
	{
		m_pWnd = pWnd ;
		m_pDC  = new CClientDC(pWnd);
		if( !m_pDC              )	return FALSE;
		if( !SetupPixelFormat() )	return FALSE;
		
		m_hRC = wglCreateContext( m_pDC->GetSafeHdc() );
		if( !m_hRC                                        )	return FALSE;	
		if( !wglMakeCurrent( m_pDC->GetSafeHdc(), m_hRC ) ) return FALSE;

		SetDefaultProperties();
		glewInit();//OpenGl extension
		return TRUE;
	}

	virtual void OnDestroy()
	{
		MakeOpenGLCurrent() ;
		wglMakeCurrent(0,0);
		wglDeleteContext( m_hRC );
		if( m_pDC )	delete m_pDC;
	}

	virtual void OnSize(int cx, int cy){ if ( cx <= 0 || cy <= 0) return; m_cx = cx; m_cy = cy;}

	virtual void MakeOpenGLCurrent() const { 
		wglMakeCurrent( m_pDC->GetSafeHdc(), m_hRC );
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	//PRE POST RENDERING////////////////////////////////////////////////////////////////////	
	void OnDraw_Begin( double view_near=0.02, double view_far = 300)
	{
		assert(!m_bIsDrawing) ;
		m_bIsDrawing = true ;
		MakeOpenGLCurrent() ;

		glViewport(0, 0, m_cx, m_cy);
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		gluPerspective(45, m_cx / (double)m_cy , view_near, view_far);

		glMatrixMode( GL_MODELVIEW ) ;
		glLoadIdentity();
		m_eyeParam.callGluLookAt();
		
		glClearColor( m_clearColor[0],m_clearColor[1],m_clearColor[2],m_clearColor[3] ) ;
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT );
	}

	void OnDraw_Begin_forOrtho(double view_near = 0, double view_far = 10)
	{
		assert(!m_bIsDrawing) ;
		m_bIsDrawing = true ;
		MakeOpenGLCurrent() ;

		//viewport
		glViewport(0, 0, m_cx, m_cy);
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		double x_y = m_cx/(double)m_cy, y_x = m_cy/(double)m_cx;
		double r = m_orthoViewSize;
		if( m_cx > m_cy ) glOrtho( -x_y * r, x_y * r, -    r,       r, view_near, view_far);
		else		      glOrtho(       -r,       r, -y_x*r, y_x * r, view_near, view_far);
		
		glMatrixMode( GL_MODELVIEW ) ;
		glLoadIdentity();
		m_eyeParam.callGluLookAt();

		glClearColor( m_clearColor[0],m_clearColor[1],m_clearColor[2],m_clearColor[3] ) ;
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT );	
	}

	void OnDraw_Begin_forOrtho_NoScaling( TVector3 &trans, double scale, double view_near = 0, double view_far = 10)
	{
		assert(!m_bIsDrawing) ;
		m_bIsDrawing = true ;
		MakeOpenGLCurrent() ;

		//viewport
		glViewport(0, 0, m_cx, m_cy);
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		glOrtho( -m_cx/2, m_cx/2, -m_cy/2, m_cy/2, view_near, view_far);
		
		glMatrixMode( GL_MODELVIEW ) ;
		glLoadIdentity();
		m_eyeParam.callGluLookAt();
		glTranslated( trans[0], trans[1], trans[2] );
		glScaled    ( scale   , scale   , scale   );

		glClearColor( m_clearColor[0],m_clearColor[1],m_clearColor[2],m_clearColor[3] ) ;
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT );	
	}

	void OnDraw_End()
	{
		glFinish();
		SwapBuffers( m_pDC->GetSafeHdc() );
		wglMakeCurrent(NULL,NULL);
		m_bIsDrawing = false ;
	}

	inline bool IsDrawing           ()const{ return m_bIsDrawing  ; }
	inline void SetMouseCapture     ()const{ m_pWnd->SetCapture (); }
	inline void ReleaseMouseCapture ()const{ ReleaseCapture     (); }
public:

	TVector3  getEyePoint() const{ return m_eyeParam.m_pos  ;}
	TVector3  getEyeFocus() const{ return m_eyeParam.m_focus;}
	TVector3  getEyeYDir () const{ return m_eyeParam.m_upDir;}
	inline void getEyePoint(TVector3 &trgt) const { trgt.Set( m_eyeParam.m_pos); }
	inline void getEyeFocus(TVector3 &trgt) const { trgt.Set( m_eyeParam.m_focus); }
	inline void getEyeYDir (TVector3 &trgt) const { trgt.Set( m_eyeParam.m_upDir); }

	void setOrthoViewSize(double r){ m_orthoViewSize = r;}
	void setInitEyePosition(const TVector3 &pos, const TVector3 &yDir, const TVector3 &focus){ m_initEyeParam.set( pos, focus, yDir);}
	void setEyePoint       (const TVector3 &pos, const TVector3 &yDir, const TVector3 &focus){ m_eyeParam    .set( pos, focus, yDir);}
	void eyePointInit(){
		m_eyeParam = m_initEyeParam;
	}
	void ButtonDownForZoom      (const CPoint& pos){m_prevpos=pos; m_ButtonState = T_BUTTONDOWN_ZOOM       ; m_pWnd->SetCapture();}
	void ButtonDownForRotate    (const CPoint& pos){m_prevpos=pos; m_ButtonState = T_BUTTONDOWN_ROTATE     ; m_pWnd->SetCapture();}
	void ButtonDownForYdirRot   (const CPoint& pos){m_prevpos=pos; m_ButtonState = T_BUTTONDOWN_ROTYDIR    ; m_pWnd->SetCapture();}
	void ButtonDownForTranslate (const CPoint& pos){m_prevpos=pos; m_ButtonState = T_BUTTONDOWN_TRANSLATE  ; m_pWnd->SetCapture();}
	void ButtonDownForOrthoZoom (const CPoint& pos){m_prevpos=pos; m_ButtonState = T_BUTTONDOWN_ORTHO_ZOOM ; m_pWnd->SetCapture();}
	void ButtonDownForOrthoTrans(const CPoint& pos){m_prevpos=pos; m_ButtonState = T_BUTTONDOWN_ORTHO_TRANS; m_pWnd->SetCapture();}

	void ButtonDownForOrtho_ZoomYRotZoom(const CPoint& pos, int wOffset= 30, int hOffset = 30){ m_prevpos=pos; m_pWnd->SetCapture();
		RECT rect ; m_pWnd->GetClientRect(&rect) ; 
		m_ButtonState = ( pos.x < wOffset || rect.right  - wOffset < pos.x ) ? T_BUTTONDOWN_ORTHO_ZOOM : 
			            ( pos.y < hOffset || rect.bottom - hOffset < pos.y ) ? T_BUTTONDOWN_ROTYDIR    : T_BUTTONDOWN_ORTHO_ZOOM;
	}
	void ButtonDownForOrtho_TransYRotZoom(const CPoint& pos, int wOffset= 30, int hOffset = 30){ m_prevpos=pos; m_pWnd->SetCapture();
		RECT rect ; m_pWnd->GetClientRect(&rect) ; 
		m_ButtonState = ( pos.x < wOffset || rect.right  - wOffset < pos.x ) ? T_BUTTONDOWN_ORTHO_ZOOM : 
			            ( pos.y < hOffset || rect.bottom - hOffset < pos.y ) ? T_BUTTONDOWN_ROTYDIR    : T_BUTTONDOWN_ORTHO_TRANS;
	}
	void ButtonDownFor_TransYRotZoom(const CPoint& pos, int wOffset= 30, int hOffset = 30){ m_prevpos=pos; m_pWnd->SetCapture();
		RECT rect ; m_pWnd->GetClientRect(&rect) ; 
		m_ButtonState = ( pos.x < wOffset || rect.right  - wOffset < pos.x ) ? T_BUTTONDOWN_ZOOM   : 
			            ( pos.y < hOffset || rect.bottom - hOffset < pos.y ) ? T_BUTTONDOWN_ROTYDIR: T_BUTTONDOWN_TRANSLATE;
	}
	void ButtonDownFor_RotYRotZoom(const CPoint& pos, int wOffset= 30, int hOffset = 30){ m_prevpos=pos; m_pWnd->SetCapture();
		RECT rect ; m_pWnd->GetClientRect(&rect) ; 
		m_ButtonState = ( pos.x < wOffset || rect.right  - wOffset < pos.x ) ? T_BUTTONDOWN_ZOOM    : 
			            ( pos.y < hOffset || rect.bottom - hOffset < pos.y ) ? T_BUTTONDOWN_ROTYDIR : T_BUTTONDOWN_ROTATE;
	}

	void ButtonUp(){
		m_ButtonState = T_BUTTON_NON;	
		if( m_pWnd == m_pWnd->GetCapture()) ReleaseCapture();
	}

	bool MouseMove(const CPoint& pos)
	{
		TVector3 &eyeP = m_eyeParam.m_pos  ;
		TVector3 &eyeF = m_eyeParam.m_focus;
		TVector3 &eyeU = m_eyeParam.m_upDir;
		if( m_ButtonState == T_BUTTON_NON || m_pWnd != m_pWnd->GetCapture()) return false;

		if( m_ButtonState == T_BUTTONDOWN_ROTATE){			
			//マウスのxyの動きに応じて視点をfocus中心に回転
			double theta = -(pos.x - m_prevpos.x) / 200.0;
			double phai  = -(pos.y - m_prevpos.y) / 200.0;
			eyeRotation( theta, phai );
		}
		else if( m_ButtonState == T_BUTTONDOWN_ZOOM)
		{
				TVector3 eyeRay = eyeF - eyeP; 
				eyeRay.Normalize_Self();
				eyeRay *= (pos.y - m_prevpos.y) * ((eyeP- eyeF).Length() / 80.0);//この分だけ動かす
				
				TVector3 newEyeP = eyeP + eyeRay;
				if( (newEyeP - eyeF).Length() > 0.02) eyeP = newEyeP;
		}
		else if( m_ButtonState == T_BUTTONDOWN_ROTYDIR )
		{
				TVector3 eyeRay = eyeF - eyeP; eyeRay.Normalize_Self();
				TMatrix16 M; M.RotateAlongArbitraryAxis( eyeRay, (pos.x - m_prevpos.x) * 0.005);
				eyeU = M * eyeU;
		}
		else if(m_ButtonState == T_BUTTONDOWN_TRANSLATE){
			TVector3 transXdir = (eyeP - eyeF) ^ eyeU;
			transXdir.Normalize_Self();

			TVector3 t;
			t.SetAdditionWithCoef( (pos.x - m_prevpos.x), transXdir, 
				                   (pos.y - m_prevpos.y), eyeU );
			t *= (eyeP- eyeF).Length() / 900.0 ;
			eyeP += t;
			eyeF += t;
		}
		else if(m_ButtonState == T_BUTTONDOWN_ORTHO_TRANS){
			TVector3 transXdir = (eyeP - eyeF) ^ eyeU;
			transXdir.Normalize_Self();
			double r = 2*m_orthoViewSize;
			double pixelWidth =  (m_cx==0&&m_cy==0)?0.01 : ( m_cx > m_cy ) ? r/ m_cy : r/ m_cx;
			TVector3 t;
			t.SetAdditionWithCoef(  ( pos.x - m_prevpos.x ) * pixelWidth, transXdir, 
				                    ( pos.y - m_prevpos.y ) * pixelWidth, eyeU );
			eyeP += t;
			eyeF += t;
		}
		else if(m_ButtonState == T_BUTTONDOWN_ORTHO_ZOOM) { m_orthoViewSize *= (1 + ( m_prevpos.y-pos.y) * 0.001); if( m_orthoViewSize < 0.01 ) m_orthoViewSize = 0.01;}
		else{}
		m_prevpos = pos;
		return TRUE ;
	}

	void eyeRotation(double theta,double phai)
	{
		TVector3 &eyeP = m_eyeParam.m_pos  ;
		TVector3 &eyeF = m_eyeParam.m_focus;
		TVector3 &eyeU = m_eyeParam.m_upDir;

		TMatrix16 rotThetaMat , rotPhaiMat;
		TVector3 axisT = eyeU;
		TVector3 axisP = (eyeF - eyeP) ^ eyeU;//こっちの方が速い
		
		rotThetaMat.RotateAlongArbitraryAxis(axisT,theta);
		rotPhaiMat .RotateAlongArbitraryAxis(axisP, phai);
		eyeU = rotPhaiMat * rotThetaMat * eyeU;
		eyeP = rotPhaiMat * rotThetaMat * (eyeP - eyeF) + eyeF;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	//PROJECTION//////////////////////////////////////////////////////////////////////////////
	inline void GetCursorRay(int cx, int cy, TVector3 &rayPos, TVector3 &rayDir) const
	{
		if( !m_bIsDrawing ) MakeOpenGLCurrent() ;
		double modelMat[16],projMat[16] ;int vp[4];
		glGetDoublev(GL_MODELVIEW_MATRIX,modelMat) ;
		glGetDoublev(GL_PROJECTION_MATRIX,projMat) ;
		glGetIntegerv(GL_VIEWPORT,vp) ;

		gluUnProject(cx, vp[3] - cy, 0.0, modelMat, projMat, vp, &rayPos[0], &rayPos[1], &rayPos[2]) ;
		gluUnProject(cx, vp[3] - cy, 0.2, modelMat, projMat, vp, &rayDir[0], &rayDir[1], &rayDir[2]) ;
		
		rayDir -= rayPos; 
		rayDir.Normalize_Self();
		if( !m_bIsDrawing ) wglMakeCurrent(NULL,NULL);
	}
	inline void GetCursorRay( const CPoint &point, TVector3 &rayPos, TVector3 &rayDir) const{
		GetCursorRay( point.x, point.y, rayPos, rayDir);
	}

	
	TVector3 unProject_correctY( double cx, double cy, float depth)
	{
		TVector3 pos;
		unProject_correctY( cx,cy,depth,pos);
		return pos;
	}

	void unProject_correctY( double cx, double cy, float depth, TVector3 &pos)
	{
		if( !m_bIsDrawing ) MakeOpenGLCurrent() ;
		double modelMat[16],projMat[16] ;int vp[4];
		glGetDoublev(GL_MODELVIEW_MATRIX,modelMat) ;
		glGetDoublev(GL_PROJECTION_MATRIX,projMat) ;
		glGetIntegerv(GL_VIEWPORT,vp) ;

		gluUnProject(cx, vp[3] - cy, depth, modelMat, projMat, vp, &pos[0], &pos[1], &pos[2]) ;
		
		if( !m_bIsDrawing ) wglMakeCurrent(NULL,NULL);
	}

	void Project( double  objx,double  objy,double  objz, double& winx,double& winy,double& winz ) const
	{
		if( !m_bIsDrawing ) MakeOpenGLCurrent();
		double model[16],proj[16] ;int vp[4] ;
		glGetDoublev(GL_MODELVIEW_MATRIX,model) ;
		glGetDoublev(GL_PROJECTION_MATRIX,proj) ;
		glGetIntegerv(GL_VIEWPORT, vp) ;

		gluProject( objx,objy,objz,model,proj,vp,&winx,&winy,&winz ) ;
		if( !m_bIsDrawing ) wglMakeCurrent(NULL,NULL);
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	//MISCS///////////////////////////////////////////////////////////////////////////////////
	inline void GetClearColor( TVector3 &color ){ color.Set( m_clearColor[0], m_clearColor[1], m_clearColor[2]  );}
	void SetClearColor( double r, double g, double b, double a){
		m_clearColor[0] = (float)r;
		m_clearColor[1] = (float)g;
		m_clearColor[2] = (float)b;
		m_clearColor[3] = (float)a;
	}


	CDC* GetDC(){	 return m_pDC; }
	int GetWidth (){ return m_cx ; }
	int GetHeight(){ return m_cy ; }
	double GetOrthoViewSize(){ return m_orthoViewSize; }

	void RedrawWindow(){ if( m_pWnd ) m_pWnd->RedrawWindow() ; }


	inline double getDistEyeToScreen()
	{
		if( !m_bIsDrawing ) MakeOpenGLCurrent() ;
		double length = 0;
		glMatrixMode( GL_MODELVIEW );
		glPushMatrix  ();
		glLoadIdentity();
		{
			double modelMat[16],projMat[16];
			int vp[4];
			glGetDoublev(GL_MODELVIEW_MATRIX,modelMat) ;
			glGetDoublev(GL_PROJECTION_MATRIX,projMat) ;
			glGetIntegerv(GL_VIEWPORT,vp) ;
			double x,y,z;
			gluUnProject(vp[0]+vp[2]/2, vp[1]+vp[3]/2, 0.0001, modelMat, projMat, vp, &x,&y,&z);
			length = fabs( z );
		}
		glPopMatrix();

		if( !m_bIsDrawing ) wglMakeCurrent(NULL,NULL);
		return length;
	}
	
	inline void t_drawEyeTransGuid(double fPtRad = 0.1, double circleRad = 2.5 , int circleRes = 30)
	{
		TVector3 &eyeP = m_eyeParam.m_pos  ;
		TVector3 &eyeF = m_eyeParam.m_focus;
		TVector3 &eyeU = m_eyeParam.m_upDir;

		if( m_ButtonState != T_BUTTONDOWN_ZOOM      && 
			m_ButtonState != T_BUTTONDOWN_ROTATE    &&
			m_ButtonState != T_BUTTONDOWN_TRANSLATE ) return;
		
		glPushMatrix();
		glTranslated( eyeF.data[0], eyeF.data[1], eyeF.data[2] );
		double r = circleRad;
		double stepTheta = 1.0 / (double) circleRes * 2 * M_PI;

		glDisable  ( GL_LIGHTING );
		glLineWidth( 2 );

		glBegin( GL_LINE_STRIP );
		glColor3d( 1,0,0); for( int i=0; i<circleRes; ++i) glVertex3d( r*cos(stepTheta*i), r*sin(stepTheta*i), 0); glVertex3d( r*cos(0.0), r*sin(0.0), 0);
		glEnd();
		glBegin( GL_LINE_STRIP );
		glColor3d( 0,1,0); for( int i=0; i<circleRes; ++i) glVertex3d( 0, r*cos(stepTheta*i), r*sin(stepTheta*i)); glVertex3d( 0, r*cos(0.0), r*sin(0.0)); 
		glEnd();
		glBegin( GL_LINE_STRIP );
		glColor3d( 0,0,1); for( int i=0; i<circleRes; ++i) glVertex3d( r*sin(stepTheta*i), 0, r*cos(stepTheta*i) ); glVertex3d( r*sin(0.0), 0, r*cos(0.0) );
		glEnd();

		r = fPtRad;
		glColor3d( 1,1,1);
		glBegin( GL_LINE_STRIP );
		for( int i=0; i<circleRes; ++i) glVertex3d( r*cos(stepTheta*i), r*sin(stepTheta*i), 0); glVertex3d( r*cos(0.0), r*sin(0.0), 0);
		for( int i=0; i<circleRes; ++i) glVertex3d( 0, r*cos(stepTheta*i), r*sin(stepTheta*i)); glVertex3d( 0, r*cos(0.0), r*sin(0.0));
		for( int i=0; i<circleRes; ++i) glVertex3d( r*sin(stepTheta*i), 0, r*cos(stepTheta*i)); glVertex3d( r*sin(0.0), 0, r*cos(0.0) );
		glEnd();

		glBegin( GL_LINES);
		glColor3d( 1,0,0); glVertex3d(0,0,0); glVertex3d( 1, 0, 0 );
		glColor3d( 0,1,0); glVertex3d(0,0,0); glVertex3d( 0, 1, 0 );
		glColor3d( 0,0,1); glVertex3d(0,0,0); glVertex3d( 0, 0, 1 );
		glEnd();

		glPopMatrix();
	}


private:
	BOOL SetupPixelFormat(){
		static PIXELFORMATDESCRIPTOR pfd = {
			sizeof(PIXELFORMATDESCRIPTOR),  // size of this pfd
				1,                              // version number
				PFD_DRAW_TO_WINDOW |            // support window
				PFD_SUPPORT_OPENGL |          // support OpenGL
				PFD_DOUBLEBUFFER,             // double buffered
				PFD_TYPE_RGBA,                  // RGBA type
				32,                             // 24-bit color depth
				0, 0, 0, 0, 0, 0,               // color bits ignored
				0,                              // no alpha buffer
				0,                              // shift bit ignored
				0,                              // no accumulation buffer
				0, 0, 0, 0,                     // accum bits ignored
				//        32,                             // 32-bit z-buffer
				16,	// NOTE: better performance with 16-bit z-buffer
				0,                              // no stencil buffer
				0,                              // no auxiliary buffer
				PFD_MAIN_PLANE,                 // main layer
				0,                              // reserved
				0, 0, 0                         // layer masks ignored
		};
		
		int pixelformat = ChoosePixelFormat(m_pDC->GetSafeHdc(), &pfd);
		if ( !pixelformat                                            ) return FALSE;
		if ( !SetPixelFormat(m_pDC->GetSafeHdc(), pixelformat, &pfd) ) return FALSE;
		return TRUE;
	}

	
	void SetDefaultProperties()
	{
		glClearColor( m_clearColor[0],m_clearColor[1],m_clearColor[2],m_clearColor[3] ) ;
		glClearDepth( 1.0f );
		glEnable( GL_DEPTH_TEST );
		glDisable( GL_BLEND );
		//glBlendFunc( GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA );
		glBlendFunc( GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA );
		//glBlendFunc( GL_SRC_ALPHA,GL_DST_ALPHA );
		//glBlendFunc( GL_ONE_MINUS_DST_ALPHA,GL_DST_ALPHA );

		GLfloat lightpos[] = { 1000,1000,-50000,1 };
		GLfloat spec[] = { 0,0,0.05f,1 } ;
		GLfloat diff[] = { 0.5f,0.5f,0.5f,1.0f };
		GLfloat amb [] = { 0.3f,0.3f,0.3f,1.0f };

		GLfloat shininess = 1.5f ;

		glMaterialfv( GL_FRONT, GL_SPECULAR,  spec );
		glMaterialfv( GL_FRONT, GL_DIFFUSE,   diff );
		glMaterialfv( GL_FRONT, GL_AMBIENT,   amb );
		glMaterialfv( GL_FRONT, GL_SHININESS, &shininess );
		glLightfv(GL_LIGHT0, GL_POSITION, lightpos);

		GLfloat light_Ambient0[] = { 0.2f , 0.2f , 0.2f , 1};
		GLfloat light_Diffuse0[] = { 1,1,1,1};
		GLfloat light_Specular0[]= { 0.2f , 0.2f , 0.2f , 1};
		glLightfv(GL_LIGHT0,GL_AMBIENT,light_Ambient0);
		glLightfv(GL_LIGHT0,GL_DIFFUSE,light_Diffuse0);
		glLightfv(GL_LIGHT0,GL_SPECULAR,light_Specular0);

		GLfloat light_Ambient1[] = { 0,0,0,1};
		GLfloat light_Diffuse1[] = { 0.5f , 0.5f , 0.5f , 1};
		GLfloat light_Specular1[]= { 0,0,0,1};
		glLightfv(GL_LIGHT1,GL_AMBIENT,light_Ambient1);
		glLightfv(GL_LIGHT1,GL_DIFFUSE,light_Diffuse1);
		glLightfv(GL_LIGHT1,GL_SPECULAR,light_Specular1);
		
		GLfloat light_Ambient2[] = { 0,0,0, 1};
		GLfloat light_Diffuse2[] = { 0.5f , 0.5f , 0.5f , 1};
		GLfloat light_Specular2[]= { 0,0,0,1};
		glLightfv(GL_LIGHT2,GL_AMBIENT,light_Ambient2);
		glLightfv(GL_LIGHT2,GL_DIFFUSE,light_Diffuse2);
		glLightfv(GL_LIGHT2,GL_SPECULAR,light_Specular2);
		

		glTexEnvi( GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE ) ;
		
		glShadeModel( GL_SMOOTH ) ;
		glPixelStorei( GL_UNPACK_ALIGNMENT,4 ) ;

		glCullFace( GL_BACK ) ;
		
		glPolygonMode( GL_FRONT,GL_FILL ) ;
		glEnable( GL_CULL_FACE ) ;
		
		m_eyeParam = m_initEyeParam;

		//Modify the following state in the rendering sequence
		//glEnable( GL_NORMALIZE ) ;
		//glShadeModel( GL_FLAT ) ;
		//glEnable( GL_LIGHT0 );
		//glEnable( GL_LIGHT1 );
		//glEnable( GL_LIGHT2 );
		//glEnable( GL_LIGHTING ) ;
		//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		//glEnable( GL_TEXTURE_2D );
	}

};

////////////////////////////////////////////////////////////////////////////////////
//Projection functions//////////////////////////////////////////////////////////////
inline void t_projectEyeRayOntoPlane(const TVector3 &eyeP  , const TVector3 &ori, 
							  const TVector3 &P     , const TVector3 &N  ,
							        TVector3 &result)
{
	if(N * ori == 0){ fprintf( stderr, "error , system can't project stroke to plane 3\n" ); return result.Set(0,0,0); }
	double t = ( ( N * P ) - ( N * eyeP ) ) / ( N * ori );
	result.SetAdditionWithCoef( 1.0, eyeP, t, ori );
}

inline TVector3 t_projectOntoBuilboard( const CPoint &point, const TOGL *ogl)
{
	static const TVector3 P( 0, 0, 0);
	TVector3 rayP, rayDir, N, result; 

	ogl->GetCursorRay( point, rayP, rayDir);
	N.SetSubtract( ogl->getEyePoint(), ogl->getEyeFocus() );
	N.Normalize_Self();

	t_projectEyeRayOntoPlane( rayP, rayDir, P, N, result );
	return result;
}

inline TVector3 t_projectOntoPlane( const CPoint &point, const TVector3 &P, const TVector3 &N, const TOGL *ogl)
{
	TVector3 rayP, rayDir, result;
	ogl->GetCursorRay( point, rayP, rayDir);
	t_projectEyeRayOntoPlane( rayP, rayDir, P, N, result );
	return result;
}

//verts2D.data[2]にはdepth値が入る
inline void t_projectVerticesOn2D(const vector<TVector3> &verts, vector<TVector3> &verts2D, const TOGL &togl)
{
	verts2D.resize( verts.size() );

	if( !togl.IsDrawing() ) togl.MakeOpenGLCurrent();
	double model[16],proj[16] ;int vp[4] ;
	glGetDoublev(GL_MODELVIEW_MATRIX,model) ;
	glGetDoublev(GL_PROJECTION_MATRIX,proj) ;
	glGetIntegerv(GL_VIEWPORT, vp) ;
	for( int i=0; i < (int) verts.size(); ++i){
		gluProject(	verts  [i].data[0],  verts  [i].data[1], verts[i].data[2], model,proj,vp, 
			       &verts2D[i].data[0], &verts2D[i].data[1], &verts2D[i].data[2]) ;
		verts2D[i].data[1] = vp[3] - verts2D[i].data[1];
	}
	if( !togl.IsDrawing() ) wglMakeCurrent(NULL,NULL);
}
//2D screenへ射影
inline void t_projectVerticesOn2D(const int vSize, const TVector3 *verts, vector<TVector3> &verts2D, const TOGL &togl)
{
	verts2D.resize( vSize );

	if( !togl.IsDrawing() ) togl.MakeOpenGLCurrent();
	double model[16],proj[16] ;int vp[4] ;
	glGetDoublev(GL_MODELVIEW_MATRIX,model) ;
	glGetDoublev(GL_PROJECTION_MATRIX,proj) ;
	glGetIntegerv(GL_VIEWPORT, vp) ;
	for( int i=0; i < vSize; ++i){
		gluProject(	verts  [i].data[0],  verts  [i].data[1], verts[i].data[2], model,proj,vp, 
			       &verts2D[i].data[0], &verts2D[i].data[1], &verts2D[i].data[2]) ;
		verts2D[i].data[1] = vp[3] - verts2D[i].data[1];
	}
	if( !togl.IsDrawing() ) wglMakeCurrent(NULL,NULL);
}



inline void t_projectVerticesOntoScreen(const vector<TVector3> &verts, vector<TVector3> &projectedVerts,  const TOGL *ogl)
{
	if( !ogl->IsDrawing() ) ogl->MakeOpenGLCurrent();
	
	int vp[4];
	double model[16],proj[16] ;

	glGetIntegerv(GL_VIEWPORT,vp) ;
	glGetDoublev(GL_MODELVIEW_MATRIX,model) ;
	glGetDoublev(GL_PROJECTION_MATRIX,proj) ;
	for( int i = 0; i < (int) verts.size(); ++i)
	{
		gluProject(  verts         [i].data[0],  verts         [i].data[1],  verts         [i].data[2], model,proj,vp,
					&projectedVerts[i].data[0], &projectedVerts[i].data[1], &projectedVerts[i].data[2] );
	}

	if( !ogl->IsDrawing() ) wglMakeCurrent(NULL,NULL);
}

inline void t_projectVerticesOntoScreen(const int vertsNum, const TVector3* verts, vector<TVector3> &projectedVerts,  const TOGL *ogl)
{
	if( !ogl->IsDrawing() ) ogl->MakeOpenGLCurrent();
	
	int vp[4];
	double model[16],proj[16] ;

	glGetIntegerv(GL_VIEWPORT,vp) ;
	glGetDoublev(GL_MODELVIEW_MATRIX,model) ;
	glGetDoublev(GL_PROJECTION_MATRIX,proj) ;
	for( int i = 0; i < vertsNum; ++i)
	{
		gluProject(  verts         [i].data[0],  verts         [i].data[1],  verts         [i].data[2], model,proj,vp,
					&projectedVerts[i].data[0], &projectedVerts[i].data[1], &projectedVerts[i].data[2] );
	}

	if( !ogl->IsDrawing() ) wglMakeCurrent(NULL,NULL);
}



inline bool t_pickVerticesInRectOnScreen  (const vector<TVector3> &verts, const double minX, const double minY, 
																			const double maxX, const double maxY, 
																			vector<int> &insideVertices,        
																			const TOGL *ogl)
{
	insideVertices.clear();
	if( !ogl->IsDrawing() ) ogl->MakeOpenGLCurrent();
	{
		int vp[4] ;
		double model[16],proj[16], depth ;
		double x, y;
		glGetIntegerv(GL_VIEWPORT,vp) ;
		glGetDoublev(GL_MODELVIEW_MATRIX,model) ;
		glGetDoublev(GL_PROJECTION_MATRIX,proj) ;
		for( int i = 0; i < (int) verts.size(); ++i)
		{
			gluProject(  verts[i].data[0], verts[i].data[1], verts[i].data[2], model,proj,vp, &x, &y, &depth );
			y = vp[3] - y;
			//内外判定
			if( minX <= x && x <= maxX && minY <= y && y <= maxY ) insideVertices.push_back( i ); 
		}
	}
	if( !ogl->IsDrawing() ) wglMakeCurrent(NULL,NULL);

	return insideVertices.size() != 0;
}


inline bool t_pickVerticesInCircleOnScreen(const vector<TVector3> &verts, const double cX, const double cY, 
																			const double rX, const double rY, 
																			vector<int> &insideVertices,        
																			const TOGL *ogl )
{
	double rXrX = 1 / (rX * rX);
	double rYrY = 1 / (rY * rY);

	insideVertices.clear();
	if( !ogl->IsDrawing() ) ogl->MakeOpenGLCurrent();
	{
		int vp[4] ;
		double model[16],proj[16], depth ;
		double x, y;
		glGetIntegerv(GL_VIEWPORT,vp) ;
		glGetDoublev(GL_MODELVIEW_MATRIX,model) ;
		glGetDoublev(GL_PROJECTION_MATRIX,proj) ;
		for( int i = 0; i < (int) verts.size(); ++i)
		{
			gluProject(  verts[i].data[0], verts[i].data[1], verts[i].data[2], model,proj,vp, &x, &y, &depth );
			y = vp[3] - y;
			//内外判定
			if( (x-cX) * (x-cX) * rXrX + (y-cY) * (y-cY) * rYrY <= 1) insideVertices.push_back( i ); 
		}
	}
	if( !ogl->IsDrawing() ) wglMakeCurrent(NULL,NULL);

	return insideVertices.size() != 0;
}

inline void t_drawCubeFrame( double r, float w, double red, double green, double blue)
{
	glDisable( GL_LIGHTING );
	glColor3d(  red, green,  blue);
	glLineWidth ( w );

	glBegin( GL_LINE_STRIP );
	glVertex3d( -r, -r, -r); glVertex3d( +r, -r, -r); glVertex3d( +r, +r, -r); glVertex3d( -r, +r, -r); glVertex3d( -r, -r, -r);
	glVertex3d( -r, -r, +r); glVertex3d( +r, -r, +r); glVertex3d( +r, +r, +r); glVertex3d( -r, +r, +r); glVertex3d( -r, -r, +r);
	glEnd();
	glBegin( GL_LINES );
	glVertex3d( +r, -r, -r); glVertex3d( +r, -r, +r);
	glVertex3d( +r, +r, -r); glVertex3d( +r, +r, +r);
	glVertex3d( -r, +r, -r); glVertex3d( -r, +r, +r); 
	glEnd();
}

inline void t_DrawCylinder(double len, double r)
{
	
	//上面
	glBegin(GL_TRIANGLE_FAN);
	glNormal3d(0, 1, 0);
	for(int i = 0; i < 6; i++) 
	{
		double a = i * M_PI / 3;
		glVertex3d(r * sin(a), len, r * cos(a));
	}
	glEnd();
	//側面
	glBegin(GL_TRIANGLE_STRIP);
	for(int i = 0; i <= 6; i++) {
		double a = - i * M_PI / 3;
		glNormal3d(sin(a), 0, cos(a));
		glVertex3d(r * sin(a), 0,   r * cos(a));
		glVertex3d(r * sin(a), len, r * cos(a));
	}
	glEnd();
	//底面
	glBegin(GL_TRIANGLE_FAN);
	glNormal3d(0, -1, 0);
	for(int i = 0; i < 6; i++) {
		double a = - i * M_PI / 3;
		glVertex3d(r * sin(a), 0, r * cos(a));
	}
	glEnd();
}

inline void t_DrawCylinder_Spike(double len, double r)
{
	
	//上面
	glBegin(GL_TRIANGLE_FAN);
	glNormal3d(0, 1, 0);
	for(int i = 0; i < 12; i++) 
	{
		double a = i * M_PI / 6;
		double R = ( i%2 ==0 )?0.8*r:r;  
		glVertex3d(R * sin(a), len, R * cos(a));
	}
	glEnd();
	//側面
	glBegin(GL_TRIANGLE_STRIP);
	for(int i = 0; i <= 12; i++) {
		double a = - i * M_PI / 6;
		double R = ( i%2 ==0 )?0.8*r:r;  
		glNormal3d(sin(a), 0, cos(a));
		glVertex3d(R * sin(a), 0,   R * cos(a));
		glVertex3d(R * sin(a), len, R * cos(a));
	}
	glEnd();
	//底面
	glBegin(GL_TRIANGLE_FAN);
	glNormal3d(0, -1, 0);
	for(int i = 0; i < 12; i++) {
		double R = ( i%2 ==0 )?0.8*r:r;  
		double a = - i * M_PI / 6;
		glVertex3d(R * sin(a), 0, R * cos(a));
	}
	glEnd();
}



inline void t_drawCone(double len, double r) {
	//座標変換
	
	glBegin(GL_TRIANGLE_FAN);
	glNormal3d(0, 1, 0);
	glVertex3d(0, len, 0);
	for(int i = 0; i <= 6; i++) {
		double a = i * M_PI / 3;
		glNormal3d(sin(a), 0, cos(a));
		glVertex3d(r * sin(a), 0, r * cos(a));
	}
	glEnd();

	glBegin(GL_TRIANGLE_FAN);
	glNormal3d(0, -1, 0);
	for(int i = 0; i < 6; i++) {
		double a = - i * M_PI / 3;
		glVertex3d(r * sin(a), 0, r * cos(a));
	}
	glEnd();
}



inline void t_drawCircleOnScreen( const CPoint &point, double r)
{
	glMatrixMode( GL_MODELVIEW );
	glPushMatrix  ();
	glLoadIdentity();

	double modelMat[16],projMat[16];
	int vp[4];
	glGetDoublev(GL_MODELVIEW_MATRIX,modelMat) ;
	glGetDoublev(GL_PROJECTION_MATRIX,projMat) ;
	glGetIntegerv(GL_VIEWPORT,vp) ;

	double x,y,z;
	gluUnProject(point.x, vp[3]-point.y, 0.0001, modelMat, projMat, vp, &x,&y,&z);
	TVector3 center(x,y,z);
		
	glBegin( GL_LINE_STRIP );
	for( int i=0; i<30; ++i) glVertex3d( center.data[0] + r * cos( 2*M_PI * i / 29.0), 
			                             center.data[1] + r * sin( 2*M_PI * i / 29.0), center.data[2]); 
	glEnd();

	glPopMatrix();
}

#endif	// __TOGL_H_INCLUDED__