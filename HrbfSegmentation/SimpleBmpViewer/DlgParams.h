#pragma once


// DlgParams ダイアログ
#include "resource.h"




class DlgParams : public CDialogEx
{
	DECLARE_DYNAMIC(DlgParams)

public:
	DlgParams(CWnd* pParent = NULL);   // 標準コンストラクター
	virtual ~DlgParams();

// ダイアログ データ
	enum { IDD = IDD_DIALOG_PARAM };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV サポート

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedCheckBfX();
	afx_msg void OnBnClickedCheckBfXxx();
	afx_msg void OnBnClickedCheckBfXlogx();
	afx_msg void OnBnClickedCheckBfXRough();
	afx_msg void OnBnClickedCheckBfXxxRough();
	afx_msg void OnBnClickedCheckBfXlogxRough();
	afx_msg void OnBnClickedCheckPolynominal0();
	afx_msg void OnBnClickedCheckPolynominal1();
	afx_msg void OnBnClickedCheckPolynominal2();
	afx_msg void OnBnClickedCheckAppCoef0();
	afx_msg void OnBnClickedCheckAppCoef001();
	afx_msg void OnBnClickedCheckAppCoef01();
	afx_msg void OnBnClickedCheckAppCoef1();

	void updateChecks();
	virtual BOOL OnInitDialog();
	afx_msg void OnBnClickedCheckPolynominalNa();
	afx_msg void OnBnClickedCheckPosConstOnbound();
	afx_msg void OnBnClickedCheckPosConstBack();
	afx_msg void OnBnClickedCheckPosConstFore();
	afx_msg void OnBnClickedCheckGradConstBack();
	afx_msg void OnBnClickedCheckGradConstOnbound();
	afx_msg void OnBnClickedCheckGradConstFore();

	afx_msg void OnBnClickedCheckRbfColor();
	afx_msg void OnBnClickedCheckRbfColorGraC();
	afx_msg void OnBnClickedCheckRbfScale();
	afx_msg void OnBnClickedCheckRbfScaleGraS();
	afx_msg void OnBnClickedCheckRbfJoint();
	afx_msg void OnBnClickedCheckRbfJointGraJ();
	afx_msg void OnBnClickedCheckHrbfColorGraC();
	afx_msg void OnBnClickedCheckHrbfScaleGraS();
	afx_msg void OnBnClickedCheckHrbfJointGraJ();
	afx_msg void OnBnClickedCheckRangeCoef0();
	afx_msg void OnBnClickedCheckRangeCoef025();
	afx_msg void OnBnClickedCheckRangeCoef05();
	afx_msg void OnBnClickedCheckRangeCoef075();
	afx_msg void OnBnClickedCheckRangeCoef10();
};
