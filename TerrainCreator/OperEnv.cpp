// OperEnv.cpp : 实现文件
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "OperEnv.h"
#include "ProgressDialog.h"
#include "interface.h"
#include "travelmap.h"


// COperEnv 对话框
extern TerrainInput currTi;
extern int needUpdate;
extern char baseTextDir[_MAX_PATH];
extern char baseHgtDir[_MAX_PATH];
extern char fsData[_MAX_PATH];
extern char workDirectory[_MAX_PATH];

extern int hTxtNum;
extern int hHgtNum;
extern unsigned int numHTI;
extern std::vector<HighNode> highTxtnode;
extern std::vector<HighNode> highHgtnode;
extern std::vector<HighTerrainInfo> allHTI;
extern int RefreshTheConfigFile();
extern int SetGlobalVars();
extern void UpdateDialog01();
extern void DoInitializationCheck();

DWORD WINAPI ThreadBatchCT (PVOID pParam);
HANDLE hThread;
int InfoFlag;
extern int MessageTo;
extern void MessageToDialog(char *msg);
extern string getCurrLODString();
extern int getCurrLodInt(const char *lodstr);

DWORD WINAPI ThreadBatchCT (PVOID pParam)
{
	int *flag=(int*)pParam;

	if(!BatchCreateTerrain())
	{
		*flag=2;
		return 0;
	}
	*flag=1;
	return 1;
}

CCharHandle::CCharHandle()
{
	//pResult=NULL;
}
CCharHandle::~CCharHandle()
{
	vector<char*>::iterator iter;
	for(iter=resV.begin();iter!=resV.end();iter++)
		delete *iter;
}

IMPLEMENT_DYNAMIC(COperEnv, CBkDialogST)

COperEnv::COperEnv(CWnd* pParent /*=NULL*/)
	: CBkDialogST(COperEnv::IDD, pParent)
	, m_Lod(_T(""))
	, m_listbox1(_T(""))
{


}

COperEnv::~COperEnv()
{
}



void COperEnv::DoDataExchange(CDataExchange* pDX)
{
	CBkDialogST::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_TREE_TERRAIN, m_tree);
	DDX_Text(pDX, IDC_EDIT_LOD, m_Lod);
 	DDX_Control(pDX, IDC_LIST_TERRAIN, m_listbox1_ctrl);
	DDX_LBString(pDX, IDC_LIST_TERRAIN, m_listbox1);
}


BEGIN_MESSAGE_MAP(COperEnv, CBkDialogST)
	ON_NOTIFY(TVN_SELCHANGED, IDC_TREE_TERRAIN, &COperEnv::OnTvnSelchangedTree)
	ON_NOTIFY(TVN_ENDLABELEDIT, IDC_TREE_TERRAIN, &COperEnv::OnTvnEndlabeleditTree)
	ON_NOTIFY(NM_DBLCLK, IDC_TREE_TERRAIN, &COperEnv::OnDblClkTree)
	ON_BN_CLICKED(IDC_REFRESH_TERRAIN, 	&COperEnv::OnBnClickedRefresh_Button)
	ON_BN_CLICKED(IDC_RESET_TERRAIN, 		&COperEnv::OnBnClickedReset_Button)
	ON_BN_CLICKED(IDC_PiLiang_Button, &COperEnv::OnBnClickedPiLiang_Button)
	ON_MESSAGE	 (WM_OPERENV,  &COperEnv::OnListSet)
	ON_MESSAGE	 (WM_OPERENV_UPDATE,  &COperEnv::OnUpdateEdit)
	ON_WM_CTLCOLOR()
	ON_WM_PAINT()
	ON_WM_TIMER()
END_MESSAGE_MAP()


// COperEnv 消息处理程序
void COperEnv::OnBnClickedRefresh_Button()
{
	/* 提示消息。*/
	MessageTo = 0;
	MessageToDialog("按下 \"地景生成\" 页面的 \"刷新\" 按钮。");

	m_tree.DeleteAllItems();
	RefreshTheConfigFile();
	initTree();
}

void COperEnv::OnBnClickedReset_Button()
{
	/* 提示消息。*/
	MessageTo = 0;
	MessageToDialog("按下 \"地景生成\" 页面的 \"重置\" 按钮。");

	hTxtNum = 0;	highTxtnode.clear();
	hHgtNum = 0;	highHgtnode.clear();
	numHTI = 0;		allHTI.clear();
	
	OnBnClickedRefresh_Button();
}


void COperEnv::OnBnClickedPiLiang_Button()
{
	UpdateParameters();
	
	/* 提示消息。*/
	MessageTo = 0;
	MessageToDialog("按下 \"地景生成\" 页面的 \"生成\" 按钮。");
	
	GetDlgItem(IDC_PiLiang_Button)->EnableWindow(0);

	InfoFlag = 0;
	hThread = CreateThread(NULL, 0, ThreadBatchCT,(void*)&InfoFlag,0,NULL);

	GetDlgItem(IDC_PiLiang_Button)->EnableWindow(1);
}

BOOL COperEnv::OnInitDialog()
{
	CDialog::OnInitDialog();

	typedef BOOL (WINAPI *lpfnSetLayeredWindowAttributes)(HWND hWnd, COLORREF crKey, BYTE bAlpha, DWORD dwFlags); 
	lpfnSetLayeredWindowAttributes SetLayeredWindowAttributes; 

	m_Refresh_Button.SubclassDlgItem(IDC_REFRESH_TERRAIN, this);
	m_Refresh_Button.SetIcon(IDI_Go);
	m_Refresh_Button.DrawTransparent();

	m_Reset_Button.SubclassDlgItem(IDC_RESET_TERRAIN, this);
	m_Reset_Button.SetIcon(IDI_Go);
	m_Reset_Button.DrawTransparent();
	
	m_PiLiang_Button.SubclassDlgItem(IDC_PiLiang_Button, this);
	m_PiLiang_Button.SetIcon(IDI_Go);
	m_PiLiang_Button.DrawTransparent();

	initTree();

	m_hAccel=::LoadAccelerators(::AfxGetResourceHandle(),MAKEINTRESOURCE(IDR_ACCELERATOR1));

	/* 下面一段代码不知道干什么用的。*/
	COLORREF maskColor=RGB(0,0,0);

#if 1	//_DEBUG
	HMODULE hUser32 = GetModuleHandle("user32.dll"); //加载动态链接库
#else
	HMODULE hUser32 = GetModuleHandle(L"user32.dll"); //加载动态链接库
#endif
	SetLayeredWindowAttributes = (lpfnSetLayeredWindowAttributes)GetProcAddress(hUser32,"SetLayeredWindowAttributes"); 

	//取得SetLayeredWindowAttributes函数指针 
	//为窗口加入WS_EX_LAYERED扩展属性
	SetWindowLong(this->GetSafeHwnd(), GWL_EXSTYLE, GetWindowLong(GetSafeHwnd(), GWL_EXSTYLE)^WS_EX_LAYERED); 

	//调用SetLayeredWinowAttributes函数
	SetLayeredWindowAttributes(this->GetSafeHwnd(), maskColor, 192, 2); 

	FreeLibrary(hUser32);   //释放动态链接库

	// TODO:  在此添加额外的初始化
	m_Lod.SetString(getCurrLODString().c_str());
	UpdateData(false);

	// 加一个定时器
	this->SetTimer(1,500,NULL);//

	return TRUE;  // return TRUE unless you set the focus to a control
}


LRESULT COperEnv::OnListSet(WPARAM wParam, LPARAM lParam)
{
	char* temp=(char*)lParam;
	CString sTemp(temp);

	// Select the next item of the currently selected one.
	m_listbox1_ctrl.AddString(sTemp);
	int nCount = m_listbox1_ctrl.GetCount();
	int curSel = m_listbox1_ctrl.SetCurSel(nCount - 1);
	
	UpdateData(TRUE);
	return 0;
}


HBRUSH COperEnv::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
	HBRUSH hbr = CDialog::OnCtlColor(pDC, pWnd, nCtlColor);

	// TODO:  在此更改 DC 的任何属性
	switch(nCtlColor)
	{
	case CTLCOLOR_STATIC:
		pDC ->SetBkMode(TRANSPARENT);
		pDC->SetTextColor(RGB(0,0,0)); 
		return (HBRUSH)::GetStockObject(HOLLOW_BRUSH); 
	}
	return hbr;
}

void COperEnv::OnPaint()
{
	CPaintDC dc(this); // device context for painting
	BOOL first = true;

	CBitmap m_bmpBK;
	m_bmpBK.LoadBitmap(IDB_BITMAP1);

	CRect rect,rc;
	GetClientRect(&rect);//获得目标尺寸，即窗口客户区的坐标
	GetWindowRect(&rc);

	BITMAP bitMap;//位图结构体
	m_bmpBK.GetBitmap(&bitMap);//获得原图片尺寸

	CDC dcMem; //目标DC
	dcMem.CreateCompatibleDC(&dc); //创建与dc兼容的内存DC
	dcMem.SelectObject(&m_bmpBK);//将位图对象m_bmpBK选入内存DC
	dc.StretchBlt(0,0,rect.Width(),rect.Height(),&dcMem,192,111,rect.Width(),rect.Height(),SRCCOPY);//bitMap.bmWidth-200,bitMap.bmHeight-250,SRCCOPY);
}

LRESULT COperEnv::OnUpdateEdit(WPARAM wParam, LPARAM lParam)
{
	UpdateParameters();
	return 0;
}


void COperEnv::UpdateParameters()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(true);
	currTi.lod = getCurrLodInt(m_Lod.GetString());

}


void COperEnv::initTree()
{
	HTREEITEM hroot,hchild;
	char buf[256];
	int i,j;

	//父节点内容
	for(i=0;i<3;i++)
	{
		switch(i)
		{
		case 0:
			{
				sprintf(buf, "用户高精度纹理图片TIF (%d)", hTxtNum);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<hTxtNum;j++)
				{
					sprintf(buf, "%s    %s", highTxtnode[j].fileName,highTxtnode[j].pathName);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 1:
			{
				sprintf(buf, "用户高精度高程图片TIF (%d)", hHgtNum);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<hHgtNum;j++)
				{
					sprintf(buf, "%s    %s", highHgtnode[j].fileName, highHgtnode[j].pathName);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 2:
			{
				sprintf(buf, "用户高精度地景制作情况 (%d)", numHTI);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<numHTI;j++)
				{
					sprintf(buf, "%s    %s    %d  %2d", allHTI[j].texFname, allHTI[j].hgtFname, allHTI[j].tileID, allHTI[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		default:
			break;
		}
	}

	m_tree.SetBkColor(RGB(247,247,255));
	m_tree.SetTextColor(RGB(0,0,255));
}


void COperEnv::OnTvnSelchangedTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMTREEVIEW pNMTreeView = reinterpret_cast<LPNMTREEVIEW>(pNMHDR);
	
	CString str=m_tree.GetItemText(pNMTreeView->itemNew.hItem);

	*pResult = 0;
}

void COperEnv::OnTvnEndlabeleditTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMTVDISPINFO pTVDispInfo = reinterpret_cast<LPNMTVDISPINFO>(pNMHDR);
	
	HTREEITEM selItem,parentItem;
	selItem = m_tree.GetSelectedItem();
	if(m_tree.ItemHasChildren(selItem))
		return;

	//获取父子索引
	int indexTree,subIndex;
	subIndex = m_tree.GetItemData(selItem);
	parentItem = m_tree.GetParentItem(selItem);
	if(parentItem == NULL) return;
	indexTree = m_tree.GetItemData(parentItem);

	if((indexTree == 0) || (indexTree == 1))
	{
		MessageBox("该内容不允许更改！");
		return;
	}

	CString strText;
	m_tree.GetEditControl()->GetWindowText(strText.GetBuffer(200),200);
	
	char *strSel = (LPSTR)(LPCTSTR)strText;
	int strLen =  strlen(strSel);
	if(strLen > 0)
	{
		if((strSel[strLen-1] != '0') && (strSel[strLen-1] != '1'))
		{
			MessageBox("错误更改！请输入0或者1");
			return;
		}
	}
	else
	{
		MessageBox("内容错误！");
		return;
	}

	if(indexTree == 2)
	{
		allHTI[subIndex].isOK = atoi(&strSel[strLen-1]);
	}

	char *strItem = (LPSTR)(LPCTSTR)m_strItem;
	int lenItem = strlen(strItem);
	strItem[lenItem-1] = strSel[strLen-1];
	strText.Format("%s", strItem);
	m_tree.SetItemText(m_tree.GetSelectedItem(),strText);
	*pResult = 0;
}

void COperEnv::OnDblClkTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	// TODO: 在此添加控件通知处理程序代码
	HTREEITEM selItem,parentItem;

	selItem = m_tree.GetSelectedItem();
	if(m_tree.ItemHasChildren(selItem))
	{
		return;
	}

	//获取父节点索引
	int indexTree;
	parentItem = m_tree.GetParentItem(selItem);
	if(parentItem == NULL) return;
	indexTree = m_tree.GetItemData(parentItem);

	if((indexTree == 0) || (indexTree == 1))
	{
		MessageBox("该内容不允许更改！");
		return;
	}

	m_strItem = m_tree.GetItemText(selItem);
	m_tree.ModifyStyle(NULL,TVS_EDITLABELS);
	m_tree.EditLabel(m_tree.GetSelectedItem());
}

void COperEnv::OnTimer(UINT_PTR nIDEvent)
{
	// 执行初始化检查
	DoInitializationCheck();

	CDialog::OnTimer(nIDEvent);
	this->KillTimer(1); 
}

