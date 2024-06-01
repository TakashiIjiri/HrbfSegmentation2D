#pragma once
#include "TOGL.h"
#include "TOGLImage.h"

#include "DlgParams.h"
#include "TTriangleMesh.h"
#include "RBFManager.h"


class TBoundConst
{
public:
	TVector2 m_pos, m_dir;
	TBoundConst(const TVector2 pos, const TVector2 dir){
		m_pos = pos;
		m_dir = dir; m_dir.Normalize_Self();
	}
};


class TInOutConst
{
public:
	TVector2 m_pos;
	bool     m_tf;
	TInOutConst(const TVector2 pos, bool inOut){
		m_pos = pos   ;
		m_tf  = inOut;
	}
};


enum SOLVE_MODE{
	RBF_COLOR      ,  RBF_SCALE     ,  RBF_JOINT,  
	RBF_JOINT_GraJ ,  RBF_SCALE_GraS,  RBF_COLOR_GraC,
	HRBF_JOINT_GraJ, HRBF_SCALE_GraS, HRBF_COLOR_GraC
};


class TCore
{
	TVector2  m_imgSize ;
	double    m_pixPitch, m_colPitch;
	DlgParams m_dlg    ;

	bool m_bAutoDrawMode;
	vector< TBoundConst > m_CPs_bound;//UserÇ…éwíËÇ≥ÇÍÇΩêßñÒ
	vector< TInOutConst > m_CPs_inOut;
	vector< RBF_PosCP   > m_posCPs   ;//é¿ç€Ç…HRBFsÇ…óòópÇ≥ÇÍÇΩï®
	vector< RBF_GradCP  > m_gradCPs  ;//é¿ç€Ç…HRBFsÇ…óòópÇ≥ÇÍÇΩï®
	vector< TVector2    > m_UVs;

	//Images & image manifold & boundary lines and so on... 
	TOGL2DImage  m_imgOriginal, m_imgSmooth, m_imgVis;
	TOGL3DImage4 m_imgVolume;

	vector< pair<TVector2, TVector2> > m_boundLines;
	TTriangleMesh m_mesh, m_imgManif, m_zeroMesh, m_visPlane;

	TCore(void);
public:
	~TCore(void);
	inline static TCore* getInst(){ static TCore p; return &p;}
	TOGL m_ogl;

	void drawScene();
	void addBoundConst( TVector2 p0,  TVector2 dir);
	void addInOutConst(vector< TVector2 > &ps,  bool tf     );

	void OnKeyDown( char nChar );

	void updateSegmentation();
	void initDlg();
	//Mode switch///////////////////////
	SOLVE_MODE m_solverMode;//1RBF,2RBFÅ§s,3RBFÅ§j,4HRBFÅ§s,5HRBFÅ§j
	double m_phai, m_theta;

	bool   m_bPosConstF, m_bPosConstB, m_bPosConstM;
	bool   m_bGraConstF, m_bGraConstB, m_bGraConstM;
	/////////////////////////////
	double m_HRBF_colorDirCoef;

};

static inline int t_getPixIdx( const TVector2 &p, const int &W, const double &pitch )
{
	int x = (int)( p[0] / pitch );
	int y = (int)( p[1] / pitch );
	return x + W*y;	
}

inline bool isCtrKeyOn  (){ return GetKeyState( VK_CONTROL ) < 0 ; }
inline bool isTabKeyOn  (){ return GetKeyState( VK_TAB     ) < 0 ; }
inline bool isSpaceKeyOn(){ return GetKeyState( VK_SPACE   ) < 0 ; }
inline bool isShiftKeyOn(){ return GetKeyState( VK_SHIFT   ) < 0 ; }
inline bool isAltKeyOn  (){ return GetKeyState( VK_MENU    ) < 0 ; }

static inline void heuColor(double inVal, byte &r, byte &g, byte &b)
{
	inVal *= -1;
	inVal = max( -1, inVal);
	inVal = min(  1, inVal);
	int v = (int) ( inVal *128 + 128);
	//if( 127 <= v && v <= 129 ) {
	//	r=0; g=0; b=0;
	//	return;
	//}
	int texSize   = 256;  
	int texSize_2 = texSize / 2;
	int ri = 512 - 1024 * v / texSize;
	int gi = 512 - 512 * abs(v - texSize_2) / texSize_2;
	int bi = -512 + 1024 * v / texSize;
	r = (byte)max(min(ri, 255), 0);
	g = (byte)max(min(gi, 255), 0);
	b = (byte)max(min(bi, 255), 0);
}

