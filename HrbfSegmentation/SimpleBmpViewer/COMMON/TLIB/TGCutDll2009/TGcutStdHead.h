#ifndef __TGCUTSTDHEAD_H__
#define __TGCUTSTDHEAD_H__


//please define
//TGCUT_STATIC_LINK 
//if you want use static lib befor including this file.

/*

#define TGCUT_STATIC_LINK
#inclide "TGcutStdHead.h"

*/


//　ライブラリファイルを読み込みます。
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


//　ヘッダーファイルを読み込みます。
#include "graph.h"
#include "TTest.h"


#endif	//__STDHEAD_H__




//何もしない static lib作成時, static lib利用時
//#define DLL_EXPORT_IMPORT_KEY   

//export .dll作成時
//#define DLL_EXPORT_IMPORT_KEY    __declspec(dllexport)//classをexportする

//import  .dll 利用時
//#define DLL_EXPORT  __declspec(dllimport)
