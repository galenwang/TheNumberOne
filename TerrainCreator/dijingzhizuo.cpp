// dijingzhizuo.cpp : ����Ӧ�ó��������Ϊ��
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "dijingzhizuoDlg.h"
#include "Splash.h" 
#include "time.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CdijingzhizuoApp

BEGIN_MESSAGE_MAP(CdijingzhizuoApp, CWinApp)
	ON_COMMAND(ID_HELP, &CWinApp::OnHelp)
END_MESSAGE_MAP()


// CdijingzhizuoApp ����

CdijingzhizuoApp::CdijingzhizuoApp()
{
	// TODO: �ڴ˴���ӹ�����룬
	// ��������Ҫ�ĳ�ʼ�������� InitInstance ��
}


// Ψһ��һ�� CdijingzhizuoApp ����

CdijingzhizuoApp theApp;


// CdijingzhizuoApp ��ʼ��

BOOL CdijingzhizuoApp::InitInstance()
{
	CCommandLineInfo cmdInfo;
	ParseCommandLine(cmdInfo);
	CSplashWnd::EnableSplashScreen(cmdInfo.m_bShowSplash);
	

	
	
	// ���һ�������� Windows XP �ϵ�Ӧ�ó����嵥ָ��Ҫ
	// ʹ�� ComCtl32.dll �汾 6 ����߰汾�����ÿ��ӻ���ʽ��
	//����Ҫ InitCommonControlsEx()�����򣬽��޷��������ڡ�
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// ��������Ϊ��������Ҫ��Ӧ�ó�����ʹ�õ�
	// �����ؼ��ࡣ
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();

	AfxEnableControlContainer();

	// ��׼��ʼ��
	// ���δʹ����Щ���ܲ�ϣ����С
	// ���տ�ִ���ļ��Ĵ�С����Ӧ�Ƴ�����
	// ����Ҫ���ض���ʼ������
	// �������ڴ洢���õ�ע�����
	// TODO: Ӧ�ʵ��޸ĸ��ַ�����
	// �����޸�Ϊ��˾����֯��
	SetRegistryKey(_T("Ӧ�ó��������ɵı���Ӧ�ó���"));

	CdijingzhizuoDlg dlg;
	m_pMainWnd = &dlg;
	INT_PTR nResponse = dlg.DoModal();
	if (nResponse == IDOK)
	{
		// TODO: �ڴ˷��ô����ʱ��
		//  ��ȷ�������رնԻ���Ĵ���
	}
	else if (nResponse == IDCANCEL)
	{
		// TODO: �ڴ˷��ô����ʱ��
		//  ��ȡ�������رնԻ���Ĵ���
	}

	// ���ڶԻ����ѹرգ����Խ����� FALSE �Ա��˳�Ӧ�ó���
	//  ����������Ӧ�ó������Ϣ�á�
	return FALSE;
}


void MessageToDialog01(char *msg)
{
	/* �ҵ���ǰ�Ի���ָ�� */
	CdijingzhizuoDlg* curDialog = (CdijingzhizuoDlg*)theApp.m_pMainWnd;
	if(curDialog == NULL)
		return ;
	CDialog *cur = curDialog->m_tabPages[0];
	if(cur == NULL)
		return;

	/* ��Ϣǰ�����ʱ��ͷ */
	struct tm *local;  
	time_t t;  
	t = time(NULL);  
	local = localtime(&t); 
	
	/* ƴװ��������Ϣ */
	char fullmsg[128];
	sprintf(fullmsg, "%02d:%02d:%02d: %s", local->tm_hour, local->tm_min, local->tm_sec, msg);
	cur->SendMessage(WM_OPERENV,0, (LPARAM)fullmsg);
}



void MessageToDialog02(char *msg)
{
	/* �ҵ���ǰ�Ի���ָ�� */
	CdijingzhizuoDlg* curDialog = (CdijingzhizuoDlg*)theApp.m_pMainWnd;
	if(curDialog == NULL)
		return ;
	CDialog *cur = curDialog->m_tabPages[1];
	if(cur == NULL)
		return;

	/* ��Ϣǰ�����ʱ��ͷ */
	struct tm *local;  
	time_t t;  
	t = time(NULL);  
	local = localtime(&t); 
	
	/* ƴװ��������Ϣ */
	char fullmsg[128];
	sprintf(fullmsg, "%02d:%02d:%02d: %s", local->tm_hour, local->tm_min, local->tm_sec, msg);
	cur->SendMessage(WM_LATLON,0,(LPARAM)fullmsg);
}


void UpdateDialog01()
{
	/* �ҵ���ǰ�Ի���ָ�� */
	CdijingzhizuoDlg* curDialog = (CdijingzhizuoDlg*)theApp.m_pMainWnd;
	if(curDialog == NULL)
		return ;
	CDialog *cur = curDialog->m_tabPages[0];
	if(cur == NULL)
		return;

	/* ����һ������  */
	cur->SendMessage(WM_OPERENV_UPDATE,0,0);
	
}

