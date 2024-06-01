#include "StdAfx.h"
#include "TCore.h"
#include "MouseListener.h"

MouseListener::~MouseListener(void)
{
}

MouseListener::MouseListener(void)
{
	m_bLButton = m_bRButton = m_bMButton = false;
	 m_bPlaceOutCP = m_bPlaceInCP = m_bPlaceBoundCP = false;
}





void MouseListener::LButtonDown( const CPoint &p){
	m_bLButton =  true;
	TCore::getInst()->m_ogl.GetCursorRay( p, m_wp0, TVector3() );
	m_wp0.data[2] = 0;
	m_wp = m_wp0;

	if( isShiftKeyOn() ){
		m_bPlaceInCP    = true;
		m_placePs.clear();
		m_placePs.push_back( TVector2( m_wp[0], m_wp[1] ));
	}else {
		m_bPlaceBoundCP = true;
	}


	TCore::getInst()->m_ogl.RedrawWindow();
}
void MouseListener::LButtonUp( const CPoint &p){
	m_bLButton =  false;

	if( m_bPlaceBoundCP ){
		TVector3 pp = (m_wp-m_wp0).Normalize();
		TCore::getInst()->addBoundConst( TVector2( m_wp0[0], m_wp0[1]) , TVector2( pp[0],pp[1]) );
	}
	if( m_bPlaceInCP ){
		TCore::getInst()->addInOutConst( m_placePs, true );
	}
	m_bPlaceOutCP = m_bPlaceInCP = m_bPlaceBoundCP = false;
	m_placePs.clear();
	TCore::getInst()->m_ogl.RedrawWindow();
}



void MouseListener::RButtonDown( const CPoint &p){
	m_bRButton = true;
	
	TCore::getInst()->m_ogl.GetCursorRay( p, m_wp0, TVector3() );
	m_wp0.data[2] = 0;
	m_wp = m_wp0;
	if( isShiftKeyOn() ){
		m_bPlaceOutCP = true;
		m_placePs.clear();
		m_placePs.push_back( TVector2( m_wp[0], m_wp[1] ));
	}
	else  TCore::getInst()->m_ogl.ButtonDownForOrthoTrans(p);
	TCore::getInst()->m_ogl.RedrawWindow();
}

void MouseListener::RButtonUp( const CPoint &p){
	m_bRButton = false;
	if( m_bPlaceOutCP ){
		TCore::getInst()->addInOutConst( m_placePs, false );
	}
	m_placePs.clear();
	m_bPlaceOutCP = m_bPlaceInCP = m_bPlaceBoundCP = false;
	TCore::getInst()->m_ogl.ButtonUp();
	TCore::getInst()->m_ogl.RedrawWindow();
}





void MouseListener::MButtonUp( const CPoint &p){
	m_bMButton = false;
	TCore::getInst()->m_ogl.ButtonUp();
	TCore::getInst()->m_ogl.RedrawWindow();
}

void MouseListener::MButtonDown( const CPoint &p){
	m_bMButton = true;
	TCore::getInst()->m_ogl.ButtonDownForOrthoZoom(p);
	TCore::getInst()->m_ogl.RedrawWindow();
}

void MouseListener::LButtonDblClk( const CPoint &p){}
void MouseListener::RButtonDblClk( const CPoint &p){}
void MouseListener::MButtonDblClk( const CPoint &p){}

void MouseListener::MouseMove( const CPoint &p){

	if( m_bPlaceInCP || m_bPlaceOutCP )
	{
		TCore::getInst()->m_ogl.GetCursorRay( p, m_wp, TVector3() );
		m_wp.data[2] = 0;
		TVector2 tmp( m_wp.data );
		if( t_distance( tmp, m_placePs.back()) > 0.2 ) m_placePs.push_back( tmp ); 
	}
	else if(      m_bLButton ){
		TCore::getInst()->m_ogl.GetCursorRay( p, m_wp, TVector3() );
		m_wp.data[2] = 0;
	}
	else if( m_bRButton ) {
		TCore::getInst()->m_ogl.MouseMove( p );
	}
	else if( m_bMButton ){
		TCore::getInst()->m_ogl.MouseMove( p );
	}

	if( m_bLButton || m_bRButton || m_bMButton ) 
		TCore::getInst()->m_ogl.RedrawWindow();
}
void MouseListener::MouseWheel( const CPoint &p, int dt ){}



void MouseListener::drawTmpManip()
{
	glPushMatrix();
	if( m_bPlaceBoundCP ){

		glTranslated( 0,0,4);
		glDisable( GL_LIGHTING );
		glPointSize( 7 );
		glBegin( GL_POINTS );
		glColor3d( 1,0,0); 
		glVertex3dv( m_wp0.data);
		glEnd();

		TVector3 dir = m_wp - m_wp0;
		dir.Normalize_Self();

		glLineWidth( 3 );
		glBegin( GL_LINES );
		glColor3d( 1,0,0); glVertex3dv( m_wp0.data);
		glColor3d( 0,1,0); glVertex3dv( (m_wp0+0.05*dir).data);
		glEnd();

	}
	if( m_bPlaceInCP || m_bPlaceOutCP )
	{
		glTranslated( 0,0,4);
		glPointSize( 7 );
		if( m_bPlaceInCP  ) glColor3d( 1,0,0);
		else                glColor3d( 0,1,0);
		glBegin( GL_POINTS);
		for( int i=0; i<(int)m_placePs.size(); ++i)
			glVertex2dv( m_placePs[i].data );
		glEnd();
		glTranslated( 0,0,-4);
	}
	glPopMatrix();
}

