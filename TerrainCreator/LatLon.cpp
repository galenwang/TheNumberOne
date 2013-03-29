// LatLon.cpp : 实现文件
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "LatLon.h"
#include "ProgressDialog.h"
#include "interface.h"
#include "travelmap.h"

extern float  mWoodCoverage;
extern float  mRoadWidthInShp;
extern float  mLightPotDistance;
extern bool CurrentLampType;

extern int shpID;
extern int bldNum;
extern std::vector<CultNode> totalShp;
extern std::vector<LocateBuilding> totalUserBld;
extern std::vector<HighCultureInfo> allHCI;
extern int RefreshTheConfigFile();
extern int SetGlobalVars();
extern void UpdateDialog01();
extern int MessageTo;
extern void MessageToDialog(char *msg);

tagLatLon LatLon;

DWORD WINAPI ThreadBatchCC (PVOID pParam);
HANDLE hLLThread;
int InfoLLFlag;

DWORD WINAPI ThreadBatchCC (PVOID pParam)
{
	int *flag=(int*)pParam;

	if(!BeginCreateCulture())
	{
		*flag=2;
		return 0;
	}
	*flag=1;
	return 1;
}


// CLatLon 对话框

IMPLEMENT_DYNAMIC(CLatLon, CBkDialogST)

CLatLon::CLatLon(CWnd* pParent /*=NULL*/)
	: CBkDialogST(CLatLon::IDD, pParent)
	, m_listbox1(_T(""))
	, roadWide(20.0)
	, lightDistance(15.0)
	, woodCoverage(4000)
{

}

CLatLon::~CLatLon()
{
}

void CLatLon::DoDataExchange(CDataExchange* pDX)
{
	CBkDialogST::DoDataExchange(pDX);
	DDX_Control(pDX, 	IDC_TREE_CULTURE, m_tree);
	DDX_Control(pDX, 	IDC_LIST_CULTURE, m_listbox1_ctrl);
	DDX_LBString(pDX, IDC_LIST_CULTURE, m_listbox1);
	DDX_Text(pDX, IDC_ROADWIDE, roadWide);
	DDV_MinMaxFloat(pDX, roadWide, 5, 1000);
	DDX_Text(pDX, IDC_LIGHTDISTANCE, lightDistance);
	DDV_MinMaxFloat(pDX, lightDistance, 0, 100);
	DDX_Text(pDX, IDC_WOODCOVERAGE, woodCoverage);
	DDV_MinMaxFloat(pDX, woodCoverage, 400, 2000000);
	DDX_Check(pDX, IDC_CHECK_LAMP, (int &)m_lampType);
}


BEGIN_MESSAGE_MAP(CLatLon, CBkDialogST)
	ON_NOTIFY(TVN_SELCHANGED, IDC_TREE_CULTURE, &CLatLon::OnTvnSelchangedTree)
	ON_NOTIFY(TVN_ENDLABELEDIT, IDC_TREE_CULTURE, &CLatLon::OnTvnEndlabeleditTree)
	ON_NOTIFY(NM_DBLCLK, IDC_TREE_CULTURE, &CLatLon::OnDblClkTree)
	ON_BN_CLICKED(IDC_REFRESH_CULTURE, 	&CLatLon::OnBnClickedRefresh_Button)
	ON_BN_CLICKED(IDC_RESET_CULTURE, 		&CLatLon::OnBnClickedReset_Button)
	ON_BN_CLICKED(IDC_CREATE_BUTTON, &CLatLon::OnBnClickedCreateButton)
	ON_BN_CLICKED(IDC_CHECK_LAMP, &CLatLon::OnBnClickedLamp)

	ON_MESSAGE	 (WM_LATLON,   &CLatLon::OnListSet)
	ON_MESSAGE	 (WM_LATLON_UPDATE,  &CLatLon::OnUpdateEdit)

	ON_WM_CTLCOLOR()
	ON_WM_PAINT()

END_MESSAGE_MAP()


// CLatLon 消息处理程序
void CLatLon::OnBnClickedRefresh_Button()
{
	/* 提示消息。*/
	MessageTo = 1;
	MessageToDialog("按下 \"文化信息\" 页面的 \"刷新\" 按钮。");

	m_tree.DeleteAllItems();
	RefreshTheConfigFile();
	initTree();
}

void CLatLon::OnBnClickedReset_Button()
{
	/* 提示消息。*/
	MessageTo = 1;
	MessageToDialog("按下 \"文化信息\" 页面的 \"重置\" 按钮。");

	shpID = 0;	totalShp.clear();
	bldNum = 0;	totalUserBld.clear();
	allHCI.clear();
	
	OnBnClickedRefresh_Button();
}




BOOL CLatLon::OnInitDialog()
{
	CDialog::OnInitDialog();

	typedef BOOL (WINAPI *lpfnSetLayeredWindowAttributes)(HWND hWnd, COLORREF crKey, BYTE bAlpha, DWORD dwFlags); 
	lpfnSetLayeredWindowAttributes SetLayeredWindowAttributes; 

	m_Refresh_Button.SubclassDlgItem(IDC_REFRESH_CULTURE, this);
	m_Refresh_Button.SetIcon(IDI_Go);
	m_Refresh_Button.DrawTransparent();

	m_Reset_Button.SubclassDlgItem(IDC_RESET_CULTURE, this);
	m_Reset_Button.SetIcon(IDI_Go);
	m_Reset_Button.DrawTransparent();
	
	m_Create_Button.SubclassDlgItem(IDC_CREATE_BUTTON, this);
	m_Create_Button.SetIcon(IDI_Go);
	m_Create_Button.DrawTransparent();

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
	UpdateData(false);
	// TODO:  在此添加额外的初始化

	return TRUE;  // return TRUE unless you set the focus to a control
	// 异常: OCX 属性页应返回 FALSE
}

LRESULT CLatLon::OnListSet(WPARAM wParam, LPARAM lParam)
{
	//((CComboBox*)GetDlgItem(IDC_LIST1))->AddString(m_sCloudType[i]);
    char* temp=(char*)lParam;
	CString sTemp(temp);

	// Select the next item of the currently selected one.
	m_listbox1_ctrl.AddString(sTemp);
	int nCount = m_listbox1_ctrl.GetCount();
	int curSel = m_listbox1_ctrl.SetCurSel(nCount - 1);
	
	UpdateData(TRUE);
	return 0;
}

HBRUSH CLatLon::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
	HBRUSH hbr = CDialog::OnCtlColor(pDC, pWnd, nCtlColor);

	// TODO:  在此更改 DC 的任何属性
	switch(nCtlColor)
	{
	case CTLCOLOR_STATIC:
		pDC ->SetBkMode(TRANSPARENT);
		pDC->SetTextColor(RGB(0,0,0)); 
		return (HBRUSH)::GetStockObject(HOLLOW_BRUSH); 
/*	case CTLCOLOR_EDIT:
		//	pDC ->SetBkColor(RGB(0,0,0));
		pDC->SetTextColor(RGB(0,0,0)); 
		hbr   =   ::CreateSolidBrush(RGB(188,197,230));   
		//pDC-> SetBkMode(TRANSPARENT); 
	return (HBRUSH)::GetStockObject(HOLLOW_BRUSH); */
	}

	// TODO:  如果默认的不是所需画笔，则返回另一个画笔
	return hbr;
}

void CLatLon::OnPaint()
{
	CPaintDC dc(this); // device context for painting
	// TODO: 在此处添加消息处理程序代码
	// 不为绘图消息调用 CDialog::OnPaint()
	CBitmap m_bmpBK;
	m_bmpBK.LoadBitmap(IDB_BITMAP1);

	CRect rect,rc;
	GetClientRect(&rect);//获得目标尺寸，即窗口客户区的坐标
	GetWindowRect(&rc);

	//pParent-> ScreenToClient(  & rect   );

	BITMAP bitMap;//位图结构体
	m_bmpBK.GetBitmap(&bitMap);//获得原图片尺寸

	CDC dcMem; //目标DC
	dcMem.CreateCompatibleDC(&dc); //创建与dc兼容的内存DC
	dcMem.SelectObject(&m_bmpBK);//将位图对象m_bmpBK选入内存DC
	dc.StretchBlt(0,0,rect.Width(),rect.Height(),&dcMem,192,111,rect.Width(),rect.Height(),SRCCOPY);//bitMap.bmWidth-200,bitMap.bmHeight-250,SRCCOPY);

}

void CLatLon::OnBnClickedCreateButton()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData();

	/* 提示消息。*/
	MessageTo = 1;
	MessageToDialog("按下 \"文化信息\" 页面的 \"生成\" 按钮。");

	CurrentLampType					 = m_lampType;		
	mRoadWidthInShp				   = roadWide;
	mLightPotDistance			   = lightDistance;
	mWoodCoverage	  			   = woodCoverage;

	GetDlgItem(IDC_CREATE_BUTTON)->EnableWindow(0);
	InfoLLFlag = 0;
	hLLThread = CreateThread(NULL, 0, ThreadBatchCC,(void*)&InfoLLFlag,0,NULL);
	GetDlgItem(IDC_CREATE_BUTTON)->EnableWindow(1);

}


void CLatLon::OnBnClickedLamp()
{
	if(((CButton*)GetDlgItem(IDC_CHECK_LAMP))-> GetCheck())
	{
		GetDlgItem(IDC_LIGHTDISTANCE)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC13)->EnableWindow(TRUE);
	}
	else
	{
		GetDlgItem(IDC_LIGHTDISTANCE)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC13)->EnableWindow(FALSE);
	}
}

LRESULT CLatLon::OnUpdateEdit(WPARAM wParam, LPARAM lParam)
{
	UpdateParameters();
	return 0;
}


void CLatLon::UpdateParameters()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(true);
}


void CLatLon::initTree()
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
				sprintf(buf, "用户文化信息原始文件SHP (%d)", shpID);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<shpID;j++)
				{
					sprintf(buf, "%s    %s   %2d", totalShp[j].fileName, totalShp[j].pathName,totalShp[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 1:
			{
				sprintf(buf, "用户机场定位情况FLT (%d)", bldNum);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<bldNum;j++)
				{
					sprintf(buf, "%s    %s   %2d", totalUserBld[j].fltFname, totalUserBld[j].fltPname,totalUserBld[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 2:
			{
				sprintf(buf, "用户文化信息完成情况 (%d)", allHCI.size());
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<allHCI.size();j++)
				{
					sprintf(buf, "%s    %9d   %2d", allHCI[j].texFname, allHCI[j].tileID, allHCI[j].isOK);
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


void CLatLon::OnTvnSelchangedTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMTREEVIEW pNMTreeView = reinterpret_cast<LPNMTREEVIEW>(pNMHDR);
	
	CString str=m_tree.GetItemText(pNMTreeView->itemNew.hItem);

	*pResult = 0;
}

void CLatLon::OnTvnEndlabeleditTree(NMHDR *pNMHDR, LRESULT *pResult)
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

	char *strItem = (LPSTR)(LPCTSTR)m_strItem;
	int lenItem = strlen(strItem);
	strItem[lenItem-1] = strSel[strLen-1];
	strText.Format("%s", strItem);
	m_tree.SetItemText(m_tree.GetSelectedItem(),strText);
	*pResult = 0;
}

void CLatLon::OnDblClkTree(NMHDR *pNMHDR, LRESULT *pResult)
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

	m_strItem = m_tree.GetItemText(selItem);
	m_tree.ModifyStyle(NULL,TVS_EDITLABELS);
	m_tree.EditLabel(m_tree.GetSelectedItem());
}

