// DlgParams.cpp : 実装ファイル
//

#include "stdafx.h"
#include "SimpleBmpViewer.h"
#include "DlgParams.h"
#include "afxdialogex.h"
#include "TCore.h"
#include "RBFManager.h"

// DlgParams ダイアログ

IMPLEMENT_DYNAMIC(DlgParams, CDialogEx)

DlgParams::DlgParams(CWnd* pParent /*=NULL*/)
	: CDialogEx(DlgParams::IDD, pParent)
{

}

DlgParams::~DlgParams()
{
}

void DlgParams::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(DlgParams, CDialogEx)
	ON_BN_CLICKED(IDC_CHECK_BF_X			, &DlgParams::OnBnClickedCheckBfX)
	ON_BN_CLICKED(IDC_CHECK_BF_XXX			, &DlgParams::OnBnClickedCheckBfXxx)
	ON_BN_CLICKED(IDC_CHECK_BF_XLOGX		, &DlgParams::OnBnClickedCheckBfXlogx)
	ON_BN_CLICKED(IDC_CHECK_BF_X_ROUGH		, &DlgParams::OnBnClickedCheckBfXRough)
	ON_BN_CLICKED(IDC_CHECK_BF_XXX_ROUGH	, &DlgParams::OnBnClickedCheckBfXxxRough)
	ON_BN_CLICKED(IDC_CHECK_BF_XLOGX_ROUGH	, &DlgParams::OnBnClickedCheckBfXlogxRough)
	ON_BN_CLICKED(IDC_CHECK_POLYNOMINAL_0	, &DlgParams::OnBnClickedCheckPolynominal0)
	ON_BN_CLICKED(IDC_CHECK_POLYNOMINAL_1	, &DlgParams::OnBnClickedCheckPolynominal1)
	ON_BN_CLICKED(IDC_CHECK_POLYNOMINAL_2	, &DlgParams::OnBnClickedCheckPolynominal2)
	ON_BN_CLICKED(IDC_CHECK_APP_COEF0		, &DlgParams::OnBnClickedCheckAppCoef0)
	ON_BN_CLICKED(IDC_CHECK_APP_COEF001		, &DlgParams::OnBnClickedCheckAppCoef001)
	ON_BN_CLICKED(IDC_CHECK_APP_COEF01		, &DlgParams::OnBnClickedCheckAppCoef01)
	ON_BN_CLICKED(IDC_CHECK_APP_COEF1		, &DlgParams::OnBnClickedCheckAppCoef1)
	ON_BN_CLICKED(IDC_CHECK_POLYNOMINAL_NA    , &DlgParams::OnBnClickedCheckPolynominalNa)
	ON_BN_CLICKED(IDC_CHECK_POS_CONST_ONBOUND , &DlgParams::OnBnClickedCheckPosConstOnbound)
	ON_BN_CLICKED(IDC_CHECK_POS_CONST_BACK    , &DlgParams::OnBnClickedCheckPosConstBack)
	ON_BN_CLICKED(IDC_CHECK_POS_CONST_FORE    , &DlgParams::OnBnClickedCheckPosConstFore)
	ON_BN_CLICKED(IDC_CHECK_GRAD_CONST_BACK   , &DlgParams::OnBnClickedCheckGradConstBack)
	ON_BN_CLICKED(IDC_CHECK_GRAD_CONST_ONBOUND, &DlgParams::OnBnClickedCheckGradConstOnbound)
	ON_BN_CLICKED(IDC_CHECK_GRAD_CONST_FORE   , &DlgParams::OnBnClickedCheckGradConstFore)

	ON_BN_CLICKED(IDC_CHECK_RBF_COLOR		, &DlgParams::OnBnClickedCheckRbfColor)
	ON_BN_CLICKED(IDC_CHECK_RBF_COLOR_GRA_C	, &DlgParams::OnBnClickedCheckRbfColorGraC)
	ON_BN_CLICKED(IDC_CHECK_RBF_SCALE		, &DlgParams::OnBnClickedCheckRbfScale)
	ON_BN_CLICKED(IDC_CHECK_RBF_SCALE_GRA_S	, &DlgParams::OnBnClickedCheckRbfScaleGraS)
	ON_BN_CLICKED(IDC_CHECK_RBF_JOINT		, &DlgParams::OnBnClickedCheckRbfJoint)
	ON_BN_CLICKED(IDC_CHECK_RBF_JOINT_GRA_J	, &DlgParams::OnBnClickedCheckRbfJointGraJ)
	ON_BN_CLICKED(IDC_CHECK_HRBF_COLOR_GRA_C, &DlgParams::OnBnClickedCheckHrbfColorGraC)
	ON_BN_CLICKED(IDC_CHECK_HRBF_SCALE_GRA_S, &DlgParams::OnBnClickedCheckHrbfScaleGraS)
	ON_BN_CLICKED(IDC_CHECK_HRBF_JOINT_GRA_J, &DlgParams::OnBnClickedCheckHrbfJointGraJ)
	ON_BN_CLICKED(IDC_CHECK_RANGE_COEF_0, &DlgParams::OnBnClickedCheckRangeCoef0)
	ON_BN_CLICKED(IDC_CHECK_RANGE_COEF_025, &DlgParams::OnBnClickedCheckRangeCoef025)
	ON_BN_CLICKED(IDC_CHECK_RANGE_COEF_05, &DlgParams::OnBnClickedCheckRangeCoef05)
	ON_BN_CLICKED(IDC_CHECK_RANGE_COEF_075, &DlgParams::OnBnClickedCheckRangeCoef075)
	ON_BN_CLICKED(IDC_CHECK_RANGE_COEF_10, &DlgParams::OnBnClickedCheckRangeCoef10)
END_MESSAGE_MAP()


// DlgParams メッセージ ハンドラー



void DlgParams::OnBnClickedCheckPolynominalNa(){RBFManager::getInst()->m_polynomMode   =-1; TCore::getInst()->updateSegmentation();updateChecks(); }
void DlgParams::OnBnClickedCheckPolynominal0(){ RBFManager::getInst()->m_polynomMode   = 0; TCore::getInst()->updateSegmentation();updateChecks(); } 
void DlgParams::OnBnClickedCheckPolynominal1(){ RBFManager::getInst()->m_polynomMode   = 1; TCore::getInst()->updateSegmentation();updateChecks(); }
void DlgParams::OnBnClickedCheckPolynominal2(){ RBFManager::getInst()->m_polynomMode   = 2; TCore::getInst()->updateSegmentation();updateChecks(); }
void DlgParams::OnBnClickedCheckBfX         (){ RBFManager::getInst()->m_basisFuncMode = 1; TCore::getInst()->updateSegmentation();updateChecks();}
void DlgParams::OnBnClickedCheckBfXxx       (){ RBFManager::getInst()->m_basisFuncMode = 2; TCore::getInst()->updateSegmentation();updateChecks();}
void DlgParams::OnBnClickedCheckBfXlogx     (){ RBFManager::getInst()->m_basisFuncMode = 3; TCore::getInst()->updateSegmentation();updateChecks();}
void DlgParams::OnBnClickedCheckBfXRough    (){ RBFManager::getInst()->m_basisFuncMode = 4; TCore::getInst()->updateSegmentation();updateChecks();}
void DlgParams::OnBnClickedCheckBfXxxRough  (){ RBFManager::getInst()->m_basisFuncMode = 5; TCore::getInst()->updateSegmentation();updateChecks();}
void DlgParams::OnBnClickedCheckBfXlogxRough(){ RBFManager::getInst()->m_basisFuncMode = 6; TCore::getInst()->updateSegmentation();updateChecks();}

void DlgParams::OnBnClickedCheckAppCoef001(){ RBFManager::getInst()->m_approxCoef = 0.001; TCore::getInst()->updateSegmentation(); updateChecks();}
void DlgParams::OnBnClickedCheckAppCoef01 (){ RBFManager::getInst()->m_approxCoef = 0.01 ; TCore::getInst()->updateSegmentation(); updateChecks();}
void DlgParams::OnBnClickedCheckAppCoef1  (){ RBFManager::getInst()->m_approxCoef = 0.1  ; TCore::getInst()->updateSegmentation(); updateChecks();}
void DlgParams::OnBnClickedCheckAppCoef0  (){ RBFManager::getInst()->m_approxCoef = 0    ; TCore::getInst()->updateSegmentation(); updateChecks();}

void DlgParams::OnBnClickedCheckPosConstOnbound(){ TCore* t=TCore::getInst(); t->m_bPosConstM = !t->m_bPosConstM; t->updateSegmentation(); updateChecks(); }
void DlgParams::OnBnClickedCheckPosConstBack   (){ TCore* t=TCore::getInst(); t->m_bPosConstB = !t->m_bPosConstB; t->updateSegmentation(); updateChecks(); }
void DlgParams::OnBnClickedCheckPosConstFore   (){ TCore* t=TCore::getInst(); t->m_bPosConstF = !t->m_bPosConstF; t->updateSegmentation(); updateChecks(); }

void DlgParams::OnBnClickedCheckGradConstOnbound(){TCore* t=TCore::getInst(); t->m_bGraConstM = !t->m_bGraConstM; t->updateSegmentation(); updateChecks();}
void DlgParams::OnBnClickedCheckGradConstBack   (){TCore* t=TCore::getInst(); t->m_bGraConstB = !t->m_bGraConstB; t->updateSegmentation(); updateChecks();}
void DlgParams::OnBnClickedCheckGradConstFore   (){TCore* t=TCore::getInst(); t->m_bGraConstF = !t->m_bGraConstF; t->updateSegmentation(); updateChecks();}





void DlgParams::updateChecks()
{
	RBFManager *r = RBFManager::getInst();
	TCore      *t = TCore::getInst();


	
	((CButton*)GetDlgItem(IDC_CHECK_RBF_COLOR       ))->SetCheck( t->m_solverMode == RBF_COLOR       );
	((CButton*)GetDlgItem(IDC_CHECK_RBF_SCALE       ))->SetCheck( t->m_solverMode == RBF_SCALE       );
	((CButton*)GetDlgItem(IDC_CHECK_RBF_JOINT       ))->SetCheck( t->m_solverMode == RBF_JOINT       );
	((CButton*)GetDlgItem(IDC_CHECK_RBF_COLOR_GRA_C ))->SetCheck( t->m_solverMode == RBF_COLOR_GraC  );
	((CButton*)GetDlgItem(IDC_CHECK_RBF_SCALE_GRA_S ))->SetCheck( t->m_solverMode == RBF_SCALE_GraS  );
	((CButton*)GetDlgItem(IDC_CHECK_RBF_JOINT_GRA_J ))->SetCheck( t->m_solverMode == RBF_JOINT_GraJ  );
	
	((CButton*)GetDlgItem(IDC_CHECK_HRBF_COLOR_GRA_C ))->SetCheck( t->m_solverMode == HRBF_COLOR_GraC  );
	((CButton*)GetDlgItem(IDC_CHECK_HRBF_SCALE_GRA_S ))->SetCheck( t->m_solverMode == HRBF_SCALE_GraS  );
	((CButton*)GetDlgItem(IDC_CHECK_HRBF_JOINT_GRA_J ))->SetCheck( t->m_solverMode == HRBF_JOINT_GraJ  );

	((CButton*)GetDlgItem(IDC_CHECK_BF_X		  ))->SetCheck( r->m_basisFuncMode==1 );
	((CButton*)GetDlgItem(IDC_CHECK_BF_X		  ))->SetCheck( r->m_basisFuncMode==1 );
	((CButton*)GetDlgItem(IDC_CHECK_BF_XXX		  ))->SetCheck( r->m_basisFuncMode==2 );
	((CButton*)GetDlgItem(IDC_CHECK_BF_XLOGX      ))->SetCheck( r->m_basisFuncMode==3 );
	((CButton*)GetDlgItem(IDC_CHECK_BF_X_ROUGH    ))->SetCheck( r->m_basisFuncMode==4 );
	((CButton*)GetDlgItem(IDC_CHECK_BF_XXX_ROUGH  ))->SetCheck( r->m_basisFuncMode==5 );
	((CButton*)GetDlgItem(IDC_CHECK_BF_XLOGX_ROUGH))->SetCheck( r->m_basisFuncMode==6 );


	((CButton*)GetDlgItem(IDC_CHECK_POLYNOMINAL_NA))->SetCheck( r->m_polynomMode==-1 );
	((CButton*)GetDlgItem(IDC_CHECK_POLYNOMINAL_0 ))->SetCheck( r->m_polynomMode==0  );
	((CButton*)GetDlgItem(IDC_CHECK_POLYNOMINAL_1 ))->SetCheck( r->m_polynomMode==1  );
	((CButton*)GetDlgItem(IDC_CHECK_POLYNOMINAL_2 ))->SetCheck( r->m_polynomMode==2  );

	((CButton*)GetDlgItem(IDC_CHECK_APP_COEF001    ))->SetCheck( r->m_approxCoef==0.001);
	((CButton*)GetDlgItem(IDC_CHECK_APP_COEF01     ))->SetCheck( r->m_approxCoef==0.01 );
	((CButton*)GetDlgItem(IDC_CHECK_APP_COEF1      ))->SetCheck( r->m_approxCoef==0.1  );
	((CButton*)GetDlgItem(IDC_CHECK_APP_COEF0      ))->SetCheck( r->m_approxCoef==0    );

	((CButton*)GetDlgItem(IDC_CHECK_POS_CONST_ONBOUND))->SetCheck( t->m_bPosConstM );
	((CButton*)GetDlgItem(IDC_CHECK_POS_CONST_BACK   ))->SetCheck( t->m_bPosConstB );
	((CButton*)GetDlgItem(IDC_CHECK_POS_CONST_FORE   ))->SetCheck( t->m_bPosConstF );

	((CButton*)GetDlgItem(IDC_CHECK_GRAD_CONST_ONBOUND))->SetCheck( t->m_bGraConstM );
	((CButton*)GetDlgItem(IDC_CHECK_GRAD_CONST_BACK   ))->SetCheck( t->m_bGraConstB );
	((CButton*)GetDlgItem(IDC_CHECK_GRAD_CONST_FORE   ))->SetCheck( t->m_bGraConstF );

	((CButton*)GetDlgItem(IDC_CHECK_RANGE_COEF_0      ))->SetCheck( t->m_HRBF_colorDirCoef==0    );
	((CButton*)GetDlgItem(IDC_CHECK_RANGE_COEF_025    ))->SetCheck( t->m_HRBF_colorDirCoef==0.1    );
	((CButton*)GetDlgItem(IDC_CHECK_RANGE_COEF_05     ))->SetCheck( t->m_HRBF_colorDirCoef==0.25    );
	((CButton*)GetDlgItem(IDC_CHECK_RANGE_COEF_075    ))->SetCheck( t->m_HRBF_colorDirCoef==0.5    );
	((CButton*)GetDlgItem(IDC_CHECK_RANGE_COEF_10     ))->SetCheck( t->m_HRBF_colorDirCoef==1.0    );


	this->RedrawWindow();
}

BOOL DlgParams::OnInitDialog()
{
	CDialogEx::OnInitDialog();
	updateChecks();

	return TRUE;  // return TRUE unless you set the focus to a control
	// 例外 : OCX プロパティ ページは必ず FALSE を返します。
}



void DlgParams::OnBnClickedCheckRbfColor     (){ TCore::getInst()->m_solverMode = RBF_COLOR; TCore::getInst()->updateSegmentation(); updateChecks();}
void DlgParams::OnBnClickedCheckRbfScale     (){ TCore::getInst()->m_solverMode = RBF_SCALE; TCore::getInst()->updateSegmentation(); updateChecks();}
void DlgParams::OnBnClickedCheckRbfJoint	 (){ TCore::getInst()->m_solverMode = RBF_JOINT; TCore::getInst()->updateSegmentation(); updateChecks();}

void DlgParams::OnBnClickedCheckRbfColorGraC (){ TCore::getInst()->m_solverMode = RBF_COLOR_GraC; TCore::getInst()->updateSegmentation(); updateChecks();}
void DlgParams::OnBnClickedCheckRbfScaleGraS (){ TCore::getInst()->m_solverMode = RBF_SCALE_GraS; TCore::getInst()->updateSegmentation(); updateChecks();}


void DlgParams::OnBnClickedCheckRbfJointGraJ (){ TCore::getInst()->m_solverMode = RBF_JOINT_GraJ; TCore::getInst()->updateSegmentation(); updateChecks();}

void DlgParams::OnBnClickedCheckHrbfColorGraC(){ TCore::getInst()->m_solverMode = HRBF_COLOR_GraC; TCore::getInst()->updateSegmentation(); updateChecks();}
void DlgParams::OnBnClickedCheckHrbfScaleGraS(){ TCore::getInst()->m_solverMode = HRBF_SCALE_GraS; TCore::getInst()->updateSegmentation(); updateChecks();}


void DlgParams::OnBnClickedCheckHrbfJointGraJ(){ TCore::getInst()->m_solverMode = HRBF_JOINT_GraJ; TCore::getInst()->updateSegmentation(); updateChecks();}


void DlgParams::OnBnClickedCheckRangeCoef0()
{
	TCore::getInst()->m_HRBF_colorDirCoef = 0;
	 TCore::getInst()->updateSegmentation(); updateChecks();
}


void DlgParams::OnBnClickedCheckRangeCoef025()
{
	TCore::getInst()->m_HRBF_colorDirCoef = 0.1;
	 TCore::getInst()->updateSegmentation(); updateChecks();

}


void DlgParams::OnBnClickedCheckRangeCoef05()
{
	TCore::getInst()->m_HRBF_colorDirCoef = 0.25;
	 TCore::getInst()->updateSegmentation(); updateChecks();
}


void DlgParams::OnBnClickedCheckRangeCoef075()
{
	TCore::getInst()->m_HRBF_colorDirCoef = 0.5;
	 TCore::getInst()->updateSegmentation(); updateChecks();
}


void DlgParams::OnBnClickedCheckRangeCoef10()
{
	TCore::getInst()->m_HRBF_colorDirCoef =1.0;
	 TCore::getInst()->updateSegmentation(); updateChecks();
}
