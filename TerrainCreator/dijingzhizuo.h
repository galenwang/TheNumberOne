// dijingzhizuo.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once
#include "NumEdit.h"

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������

#define WM_OPERENV  		WM_USER+1 
#define WM_LATLON   		WM_USER+2 
#define WM_OPERENV_UPDATE  	WM_USER+3
#define WM_GENEAPT          WM_USER+4
#define WM_3DBUILD          WM_USER+5
#define WM_LATLON_UPDATE  	WM_USER+6

// CdijingzhizuoApp:
// �йش����ʵ�֣������ dijingzhizuo.cpp
//

class CdijingzhizuoApp : public CWinApp
{
public:
	CdijingzhizuoApp();


// ��д
	public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern CdijingzhizuoApp theApp;