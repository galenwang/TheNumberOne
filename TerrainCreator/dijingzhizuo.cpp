// dijingzhizuo.cpp : 定义应用程序的类行为。
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


// CdijingzhizuoApp 构造

CdijingzhizuoApp::CdijingzhizuoApp()
{
	// TODO: 在此处添加构造代码，
	// 将所有重要的初始化放置在 InitInstance 中
}


// 唯一的一个 CdijingzhizuoApp 对象

CdijingzhizuoApp theApp;


// CdijingzhizuoApp 初始化

BOOL CdijingzhizuoApp::InitInstance()
{
	CCommandLineInfo cmdInfo;
	ParseCommandLine(cmdInfo);
	CSplashWnd::EnableSplashScreen(cmdInfo.m_bShowSplash);
	

	
	
	// 如果一个运行在 Windows XP 上的应用程序清单指定要
	// 使用 ComCtl32.dll 版本 6 或更高版本来启用可视化方式，
	//则需要 InitCommonControlsEx()。否则，将无法创建窗口。
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// 将它设置为包括所有要在应用程序中使用的
	// 公共控件类。
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();

	AfxEnableControlContainer();

	// 标准初始化
	// 如果未使用这些功能并希望减小
	// 最终可执行文件的大小，则应移除下列
	// 不需要的特定初始化例程
	// 更改用于存储设置的注册表项
	// TODO: 应适当修改该字符串，
	// 例如修改为公司或组织名
	SetRegistryKey(_T("应用程序向导生成的本地应用程序"));

	CdijingzhizuoDlg dlg;
	m_pMainWnd = &dlg;
	INT_PTR nResponse = dlg.DoModal();
	if (nResponse == IDOK)
	{
		// TODO: 在此放置处理何时用
		//  “确定”来关闭对话框的代码
	}
	else if (nResponse == IDCANCEL)
	{
		// TODO: 在此放置处理何时用
		//  “取消”来关闭对话框的代码
	}

	// 由于对话框已关闭，所以将返回 FALSE 以便退出应用程序，
	//  而不是启动应用程序的消息泵。
	return FALSE;
}


void MessageToDialog01(char *msg)
{
	/* 找到当前对话框指针 */
	CdijingzhizuoDlg* curDialog = (CdijingzhizuoDlg*)theApp.m_pMainWnd;
	if(curDialog == NULL)
		return ;
	CDialog *cur = curDialog->m_tabPages[0];
	if(cur == NULL)
		return;

	/* 消息前面加上时间头 */
	struct tm *local;  
	time_t t;  
	t = time(NULL);  
	local = localtime(&t); 
	
	/* 拼装完整的消息 */
	char fullmsg[128];
	sprintf(fullmsg, "%02d:%02d:%02d: %s", local->tm_hour, local->tm_min, local->tm_sec, msg);
	cur->SendMessage(WM_OPERENV,0, (LPARAM)fullmsg);
}



void MessageToDialog02(char *msg)
{
	/* 找到当前对话框指针 */
	CdijingzhizuoDlg* curDialog = (CdijingzhizuoDlg*)theApp.m_pMainWnd;
	if(curDialog == NULL)
		return ;
	CDialog *cur = curDialog->m_tabPages[1];
	if(cur == NULL)
		return;

	/* 消息前面加上时间头 */
	struct tm *local;  
	time_t t;  
	t = time(NULL);  
	local = localtime(&t); 
	
	/* 拼装完整的消息 */
	char fullmsg[128];
	sprintf(fullmsg, "%02d:%02d:%02d: %s", local->tm_hour, local->tm_min, local->tm_sec, msg);
	cur->SendMessage(WM_LATLON,0,(LPARAM)fullmsg);
}


void UpdateDialog01()
{
	/* 找到当前对话框指针 */
	CdijingzhizuoDlg* curDialog = (CdijingzhizuoDlg*)theApp.m_pMainWnd;
	if(curDialog == NULL)
		return ;
	CDialog *cur = curDialog->m_tabPages[0];
	if(cur == NULL)
		return;

	/* 更新一下数据  */
	cur->SendMessage(WM_OPERENV_UPDATE,0,0);
	
}

