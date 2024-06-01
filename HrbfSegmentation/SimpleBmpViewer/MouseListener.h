#pragma once

#include "TMath.h"



class MouseListener
{
	TVector3 m_wp0, m_wp;// world points
	bool m_bLButton, m_bRButton, m_bMButton;
	bool m_bPlaceOutCP, m_bPlaceInCP, m_bPlaceBoundCP;
	vector< TVector2 > m_placePs;
	MouseListener(void);
public:
	~MouseListener(void);

	inline static MouseListener* getInst(){ static MouseListener p; return &p;}

	void LButtonUp( const CPoint &p);
	void RButtonUp( const CPoint &p);
	void MButtonUp( const CPoint &p);

	void LButtonDown( const CPoint &p);
	void RButtonDown( const CPoint &p);
	void MButtonDown( const CPoint &p);
	void LButtonDblClk( const CPoint &p);
	void RButtonDblClk( const CPoint &p);
	void MButtonDblClk( const CPoint &p);
	void MouseMove( const CPoint &p);
	void MouseWheel( const CPoint &p, int dt );


	void drawTmpManip();
};

