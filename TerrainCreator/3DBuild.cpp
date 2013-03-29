// 3DBuild.cpp : 实现文件
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "3DBuild.h"
#include "travelmap.h"


// C3DBuild 对话框
extern bool rebuildTerrainForBuilding;
#define USETHREAD
DWORD WINAPI ThreadCreate3D (PVOID pParam);
HANDLE h3dThread;
int G3dInfoFlag;

DWORD WINAPI ThreadCreate3D (PVOID pParam)
{
	int *flag=(int*)pParam;

	if(!LocateUserAirport())
	{
		*flag=2;
		return 0;
	}
	*flag=1;
	return 1;
}


IMPLEMENT_DYNAMIC(C3DBuild, CDialog)

C3DBuild::C3DBuild(CWnd* pParent /*=NULL*/)
	: CDialog(C3DBuild::IDD, pParent)
	, m_listbox1(_T(""))
	, m_rebuild(false)
{

}

C3DBuild::~C3DBuild()
{
}

void C3DBuild::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_Batch_Button, m_Batch_Button);
	DDX_Control(pDX, IDC_LIST1, m_listbox1_ctrl);
	DDX_LBString(pDX, IDC_LIST1, m_listbox1);
	DDX_Check(pDX, IDC_CHECK1, m_rebuild);
}


BEGIN_MESSAGE_MAP(C3DBuild, CDialog)
	ON_BN_CLICKED(IDC_Batch_Button, &C3DBuild::OnBnClickedBatch)
	ON_WM_CTLCOLOR()
	ON_WM_PAINT()
	ON_WM_ERASEBKGND()
	ON_MESSAGE	 (WM_3DBUILD,  &C3DBuild::OnListSet)
END_MESSAGE_MAP()


// C3DBuild 消息处理程序
void C3DBuild::OnBnClickedBatch()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData();
	rebuildTerrainForBuilding       = m_rebuild;

#ifdef USETHREAD
	G3dInfoFlag = 0;
	h3dThread = CreateThread(NULL, 0, ThreadCreate3D, (void*)&G3dInfoFlag,0,NULL);
#else
	LocateUserAirport();
#endif
}


HBRUSH C3DBuild::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
	HBRUSH hbr = CDialog::OnCtlColor(pDC, pWnd, nCtlColor);

	// TODO:  如果默认的不是所需画笔，则返回另一个画笔
	switch(nCtlColor)
	{
	case CTLCOLOR_STATIC:
		pDC ->SetBkMode(TRANSPARENT);
		pDC->SetTextColor(RGB(0,0,0)); 
		return (HBRUSH)::GetStockObject(HOLLOW_BRUSH); 
	}


	return hbr;
}

LRESULT C3DBuild::OnListSet(WPARAM wParam, LPARAM lParam)
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

void C3DBuild::OnPaint()
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

BOOL C3DBuild::OnEraseBkgnd(CDC* pDC)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值

	return true;//CDialog::OnEraseBkgnd(pDC);
}

BOOL C3DBuild::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  在此添加额外的初始化
	//m_Batch_Button.SubclassDlgItem(IDC_Batch_Button, this);
	m_Batch_Button.SetIcon(IDI_Go);
	m_Batch_Button.DrawTransparent();
	return TRUE;  // return TRUE unless you set the focus to a control
	// 异常: OCX 属性页应返回 FALSE
}

