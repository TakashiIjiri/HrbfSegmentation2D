// dllmain.cpp : DLL �̏��������[�`�����`���܂��B
//

#include "stdafx.h"
#include <afxwin.h>
#include <afxdllx.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

static AFX_EXTENSION_MODULE TGcutDllDLL = { NULL, NULL };

extern "C" int APIENTRY
DllMain(HINSTANCE hInstance, DWORD dwReason, LPVOID lpReserved)
{
	// lpReserved ���g���ꍇ�͂������폜���Ă�������
	UNREFERENCED_PARAMETER(lpReserved);

	if (dwReason == DLL_PROCESS_ATTACH)
	{
		TRACE0("TGcutDll.DLL Initializing!\n");
		
		// �g�� DLL �� 1 �񂾂����������܂��B
		if (!AfxInitExtensionModule(TGcutDllDLL, hInstance))
			return 0;

		// ���� DLL �����\�[�X �`�F�[���֑}�����܂��B
		// ���� : ���̊g�� DLL ���ÖٓI�ɁAMFC �A�v���P�[�V�����ł͂Ȃ�
		//  ActiveX �R���g���[���Ȃǂ� MFC �W�� DLL �ɂ���ă����N����Ă���ꍇ�A
		//  �ȉ��̍s�� DllMain ����폜����
		//  ����폜���āA���̊g�� DLL ����G�N�X�|�[�g
		//  �z�u���Ă��������B���������āA���̊g�� DLL ���g���W�� DLL �́A
		//  ���̊֐��𖾎��I�ɌĂяo���āA
		//  �����������邽�߂ɖ����I�ɂ��̊֐����Ăяo���܂��B
		//  ����ȊO�̏ꍇ�́ACDynLinkLibrary �I�u�W�F�N�g��
		//  �W�� DLL �̃��\�[�X �`�F�[���փA�^�b�`���ꂸ�A
		//  ���̌��ʏd��Ȗ��ƂȂ�܂��B

		new CDynLinkLibrary(TGcutDllDLL);

	}
	else if (dwReason == DLL_PROCESS_DETACH)
	{
		TRACE0("TGcutDll.DLL Terminating!\n");

		// �f�X�g���N�^���Ăяo�����O�Ƀ��C�u�������I�����܂�
		AfxTermExtensionModule(TGcutDllDLL);
	}
	return 1;   // OK
}