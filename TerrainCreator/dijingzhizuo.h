// dijingzhizuo.h : PROJECT_NAME 应用程序的主头文件
//

#pragma once
#include "NumEdit.h"

#ifndef __AFXWIN_H__
	#error "在包含此文件之前包含“stdafx.h”以生成 PCH 文件"
#endif

#include "resource.h"		// 主符号

#define WM_OPERENV  		WM_USER+1 
#define WM_LATLON   		WM_USER+2 
#define WM_OPERENV_UPDATE  	WM_USER+3
#define WM_GENEAPT          WM_USER+4
#define WM_3DBUILD          WM_USER+5
#define WM_LATLON_UPDATE  	WM_USER+6

// CdijingzhizuoApp:
// 有关此类的实现，请参阅 dijingzhizuo.cpp
//

class CdijingzhizuoApp : public CWinApp
{
public:
	CdijingzhizuoApp();


// 重写
	public:
	virtual BOOL InitInstance();

// 实现

	DECLARE_MESSAGE_MAP()
};

extern CdijingzhizuoApp theApp;