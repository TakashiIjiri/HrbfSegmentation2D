
// SimpleBmpViewerView.cpp : implementation of the CSimpleBmpViewerView class
//

#include "stdafx.h"
// SHARED_HANDLERS can be defined in an ATL project implementing preview, thumbnail
// and search filter handlers and allows sharing of document code with that project.
#ifndef SHARED_HANDLERS
#include "SimpleBmpViewer.h"
#endif

#include "SimpleBmpViewerDoc.h"
#include "SimpleBmpViewerView.h"
#include "TCore.h"
#include "MouseListener.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CSimpleBmpViewerView

IMPLEMENT_DYNCREATE(CSimpleBmpViewerView, CView)

BEGIN_MESSAGE_MAP(CSimpleBmpViewerView, CView)
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CSimpleBmpViewerView::OnFilePrintPreview)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_MBUTTONDOWN()
	ON_WM_MBUTTONDBLCLK()
	ON_WM_LBUTTONUP()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONDBLCLK()
	ON_WM_MBUTTONUP()
	ON_WM_SIZE()
	ON_WM_MOUSEHWHEEL()
	ON_WM_MOUSEMOVE()
	ON_WM_NCLBUTTONDBLCLK()
	ON_WM_NCLBUTTONDOWN()
	ON_WM_NCLBUTTONUP()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONDBLCLK()
	ON_WM_ERASEBKGND()
	ON_WM_KEYDOWN()
	ON_WM_MOUSEWHEEL()
END_MESSAGE_MAP()

// CSimpleBmpViewerView construction/destruction

CSimpleBmpViewerView::CSimpleBmpViewerView()
{
	// TODO: add construction code here
}

CSimpleBmpViewerView::~CSimpleBmpViewerView()
{
}

BOOL CSimpleBmpViewerView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CSimpleBmpViewerView drawing

void CSimpleBmpViewerView::OnDraw(CDC* /*pDC*/)
{
	CSimpleBmpViewerDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	glLineWidth( 3 );
	glColor3d( 1,1,0);
	glDisable( GL_LIGHTING );
	TCore::getInst()->m_ogl.OnDraw_Begin_forOrtho( -10, 30);
	//TCore::getInst()->m_ogl.OnDraw_Begin();
	TCore::getInst()->drawScene();
	MouseListener::getInst()->drawTmpManip();

	TCore::getInst()->m_ogl.OnDraw_End();
}

//‰æ‘œ‚Ì•\Ž¦
//‰æ‘œ‚ðpick
//graph cut
//RBF



// CSimpleBmpViewerView printing


void CSimpleBmpViewerView::OnFilePrintPreview()
{
#ifndef SHARED_HANDLERS
	AFXPrintPreview(this);
#endif
}

BOOL CSimpleBmpViewerView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CSimpleBmpViewerView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CSimpleBmpViewerView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}


void CSimpleBmpViewerView::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CSimpleBmpViewerView diagnostics

#ifdef _DEBUG
void CSimpleBmpViewerView::AssertValid() const
{
	CView::AssertValid();
}

void CSimpleBmpViewerView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CSimpleBmpViewerDoc* CSimpleBmpViewerView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CSimpleBmpViewerDoc)));
	return (CSimpleBmpViewerDoc*)m_pDocument;
}
#endif //_DEBUG


// CSimpleBmpViewerView message handlers


int CSimpleBmpViewerView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	TCore::getInst()->m_ogl.OnCreate( this );
	return 0;
}


void CSimpleBmpViewerView::OnDestroy(){
	CView::OnDestroy();
	TCore::getInst()->m_ogl.OnDestroy();
}
void CSimpleBmpViewerView::OnSize(UINT nType, int cx, int cy)
{
	CView::OnSize(nType, cx, cy);
	TCore::getInst()->m_ogl.OnSize( cx, cy );
}
void CSimpleBmpViewerView::OnMButtonDown  (UINT nFlags, CPoint point){ MouseListener::getInst()->MButtonDown  ( point ); CView::OnMButtonDown(   nFlags, point);}
void CSimpleBmpViewerView::OnMButtonDblClk(UINT nFlags, CPoint point){ MouseListener::getInst()->MButtonDblClk( point ); CView::OnMButtonDblClk( nFlags, point);}
void CSimpleBmpViewerView::OnMButtonUp(    UINT nFlags, CPoint point){ MouseListener::getInst()->MButtonUp(     point ); CView::OnMButtonUp(     nFlags, point);}
void CSimpleBmpViewerView::OnLButtonDown  (UINT nFlags, CPoint point){ MouseListener::getInst()->LButtonDown  ( point ); CView::OnLButtonDown   (nFlags, point);}
void CSimpleBmpViewerView::OnLButtonDblClk(UINT nFlags, CPoint point){ MouseListener::getInst()->LButtonDblClk( point ); CView::OnLButtonDblClk( nFlags, point);}
void CSimpleBmpViewerView::OnLButtonUp    (UINT nFlags, CPoint point){ MouseListener::getInst()->LButtonUp    ( point ); CView::OnLButtonUp    ( nFlags, point);}
void CSimpleBmpViewerView::OnRButtonDown(  UINT nFlags, CPoint point){ MouseListener::getInst()->RButtonDown  ( point ); CView::OnRButtonDown(   nFlags, point);}
void CSimpleBmpViewerView::OnRButtonDblClk(UINT nFlags, CPoint point){ MouseListener::getInst()->RButtonDblClk( point ); CView::OnRButtonDblClk( nFlags, point);}
void CSimpleBmpViewerView::OnRButtonUp    (UINT nFlags, CPoint point){ MouseListener::getInst()->RButtonUp    ( point ); }

void CSimpleBmpViewerView::OnMouseMove(UINT nFlags, CPoint point)
{
	MouseListener::getInst()->MouseMove( point );
	CView::OnMouseMove(nFlags, point);
}
void CSimpleBmpViewerView::OnMouseHWheel(UINT nFlags, short zDelta, CPoint pt)
{
	MouseListener::getInst()->MouseWheel( pt, zDelta );
	CView::OnMouseHWheel(nFlags, zDelta, pt);
}








void CSimpleBmpViewerView::OnNcLButtonDblClk(UINT nHitTest, CPoint point)
{
	CView::OnNcLButtonDblClk(nHitTest, point);
}


void CSimpleBmpViewerView::OnNcLButtonDown(UINT nHitTest, CPoint point)
{
	CView::OnNcLButtonDown(nHitTest, point);
}


void CSimpleBmpViewerView::OnNcLButtonUp(UINT nHitTest, CPoint point){
	CView::OnNcLButtonUp(nHitTest, point);
}


BOOL CSimpleBmpViewerView::OnEraseBkgnd(CDC* pDC)
{
	return TRUE;
	//return CView::OnEraseBkgnd(pDC);
}


void CSimpleBmpViewerView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	TCore::getInst()->OnKeyDown(nChar);

	CView::OnKeyDown(nChar, nRepCnt, nFlags);
}


BOOL CSimpleBmpViewerView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	if( isCtrKeyOn() ) TCore::getInst()->m_phai  += zDelta > 0 ? 3 : -3;
	else               TCore::getInst()->m_theta += zDelta > 0 ? 3 : -3;

	TCore::getInst()->m_ogl.RedrawWindow();
	return CView::OnMouseWheel(nFlags, zDelta, pt);
}
