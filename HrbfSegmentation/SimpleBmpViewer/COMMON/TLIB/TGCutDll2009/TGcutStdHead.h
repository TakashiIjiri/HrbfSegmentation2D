#ifndef __TGCUTSTDHEAD_H__
#define __TGCUTSTDHEAD_H__


//please define
//TGCUT_STATIC_LINK 
//if you want use static lib befor including this file.

/*

#define TGCUT_STATIC_LINK
#inclide "TGcutStdHead.h"

*/


//�@���C�u�����t�@�C����ǂݍ��݂܂��B
#ifdef _WIN64

#ifdef TGCUT_STATIC_LINK 
#pragma comment(lib, "TGcutDlls64.lib")
#define DLL_EXPORT_IMPORT_KEY   
#else
#pragma comment(lib, "TGcutDll64.lib")
#define DLL_EXPORT_IMPORT_KEY  __declspec(dllimport)
#endif



#else


#ifdef TGCUT_STATIC_LINK 
#pragma comment(lib, "TGcutDlls.lib")
#define DLL_EXPORT_IMPORT_KEY   
#else
#pragma comment(lib, "TGcutDll.lib")
#define DLL_EXPORT_IMPORT_KEY  __declspec(dllimport)
#endif


#endif


//�@�w�b�_�[�t�@�C����ǂݍ��݂܂��B
#include "graph.h"
#include "TTest.h"


#endif	//__STDHEAD_H__




//�������Ȃ� static lib�쐬��, static lib���p��
//#define DLL_EXPORT_IMPORT_KEY   

//export .dll�쐬��
//#define DLL_EXPORT_IMPORT_KEY    __declspec(dllexport)//class��export����

//import  .dll ���p��
//#define DLL_EXPORT  __declspec(dllimport)
