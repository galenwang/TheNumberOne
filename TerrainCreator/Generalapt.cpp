// Generalapt.cpp : 实现文件
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "Generalapt.h"
#include "OperEnv.h"
#include "travelmap.h"


// CGeneralapt 对话框
extern GeneralAirport currGA;
extern bool rebuildTerrainForGeneapt;
extern bool CurrentCreateGeneralAirport;
#define USETHREAD
DWORD WINAPI ThreadCreateGA (PVOID pParam);
HANDLE hGaThread;
int GaInfoFlag;

DWORD WINAPI ThreadCreateGA (PVOID pParam)
{
	int *flag=(int*)pParam;

	if(!CreateGeneralAirport())
	{
		*flag=2;
		return 0;
	}
	*flag=1;
	return 1;
}

IMPLEMENT_DYNAMIC(CGeneralapt, CBkDialogST)

CGeneralapt::CGeneralapt(CWnd* pParent /*=NULL*/)
	: CBkDialogST(CGeneralapt::IDD, pParent)
	, shRunwayID(5)
	, dLatitude(43.45)
	, dLongtitude(105.23)
	, fAltitude(0.0)
	, cICAO(_T("A"))
	, m_listbox1(_T(""))	
	, m_createApt(false)
{

}

CGeneralapt::~CGeneralapt()
{
}

void CGeneralapt::DoDataExchange(CDataExchange* pDX)
{
	CBkDialogST::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_shRunwayID, shRunwayID);
	DDX_Text(pDX, IDC_dLatitude, dLatitude);
	DDX_Text(pDX, IDC_dLongtitude, dLongtitude);
	DDX_Text(pDX, IDC_fAltitude, fAltitude);
	//DDX_Text(pDX, IDC_fHeading, fHeading);
	//DDX_Text(pDX, IDC_fVLAAngles, fVLAAngles);
	//DDX_Text(pDX, IDC_fRunwaySlope, fRunwaySlope);
	DDX_Control(pDX, IDC_shRunwayWidth, shRunwayWidth_ctrl);
	DDX_Control(pDX, IDC_shRunwayLength, shRunwayLength_ctrl);
	//DDX_Control(pDX, IDC_shRunwayCondition, shRunwayCondition_ctrl);
	DDX_Control(pDX, IDC_shBuildingPosition, shBuildingPosition_ctrl);	
	//DDX_Control(pDX, IDC_shBuildingTypes, shBuildingTypes_ctrl);
	//DDX_Control(pDX, IDC_shAppLightTypes, shAppLightTypes_ctrl);
	//DDX_Control(pDX, IDC_shVLATypes, shVLATypes_ctrl);
	DDX_Text(pDX, IDC_cICAO, cICAO);
	DDX_Control(pDX, IDC_cICAO, cICAO_ctrl);
	DDX_Control(pDX, IDC_OK, m_OK);
	DDX_Control(pDX, IDC_LIST1, m_listbox1_ctrl);
	DDX_LBString(pDX, IDC_LIST1, m_listbox1);
	DDX_Check(pDX, IDC_CHECK_AIRPORT, (int &)m_createApt);
}


BEGIN_MESSAGE_MAP(CGeneralapt, CBkDialogST)
	ON_BN_CLICKED(IDC_OK, &CGeneralapt::OnBnClickedOK)

	ON_EN_KILLFOCUS(IDC_cICAO, &CGeneralapt::OnEnKillfocuscICAO)
	ON_EN_SETFOCUS(IDC_cICAO, &CGeneralapt::OnEnSetfocuscICAO)
	ON_WM_CTLCOLOR()
	ON_WM_ERASEBKGND()
	ON_WM_PAINT()

	ON_MESSAGE	 (WM_GENEAPT,  &CGeneralapt::OnListSet)
	ON_BN_CLICKED(IDC_CHECK_AIRPORT, &CGeneralapt::OnBnClickedCheckAirport)
END_MESSAGE_MAP()


// CGeneralapt 消息处理程序

void CGeneralapt::OnBnClickedOK()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData();

	currGA.shGenericAirportState	= 1										;
	currGA.shRunwayID             	= shRunwayID                            ;
	currGA.shRunwayWidth          	= shRunwayWidth_ctrl.GetCurSel()        ;
	//currGA.shRunwayCondition      	= shRunwayCondition_ctrl.GetCurSel()    ;
	currGA.shRunwayLength         	= shRunwayLength_ctrl.GetCurSel()       ;
	currGA.shBuildingPosition     	= shBuildingPosition_ctrl.GetCurSel()   ;
	//currGA.shBuildingTypes        	= shBuildingTypes_ctrl.GetCurSel()      ;
	//currGA.shAppLightTypes        	= shAppLightTypes_ctrl.GetCurSel()      ;
	//currGA.shVLATypes             	= shVLATypes_ctrl.GetCurSel()           ;
	currGA.dLatitude              	= dLatitude                             ;
	currGA.dLongtitude            	= dLongtitude                           ;
	currGA.fAltitude              	= fAltitude                             ;
	//currGA.fHeading               	= fHeading                              ;
	//currGA.fVLAAngles				= fVLAAngles		                    ;
	//currGA.fRunwaySlope           	= fRunwaySlope                          ;
	sprintf(currGA.cICAO, "GEN%s", (LPCSTR)(cICAO.GetString()))                                   ;
	//currGA.cLocation[64]			= cLocation
	rebuildTerrainForGeneapt       = false;
	CurrentCreateGeneralAirport    = m_createApt;
	

#ifdef USETHREAD
	GaInfoFlag = 0;
	hGaThread = CreateThread(NULL, 0, ThreadCreateGA,(void*)&GaInfoFlag,0,NULL);
#else
	CreateGeneralAirport();
#endif

}


BOOL CGeneralapt::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  在此添加额外的初始化
	m_cEdit[0].SubclassDlgItem(IDC_shRunwayID,this);
	m_cEdit[1].SubclassDlgItem(IDC_dLatitude,this);
	m_cEdit[2].SubclassDlgItem(IDC_dLongtitude,this);
	m_cEdit[3].SubclassDlgItem(IDC_fAltitude,this);
	//m_cEdit[4].SubclassDlgItem(IDC_fHeading,this);
	//m_cEdit[5].SubclassDlgItem(IDC_fVLAAngles,this);
	//m_cEdit[6].SubclassDlgItem(IDC_fRunwaySlope,this);
	cICAO_ctrl.SetLimitText(1);
	
	shRunwayWidth_ctrl.SetCurSel(0);       
	//shRunwayCondition_ctrl.SetCurSel(0);
	shRunwayLength_ctrl.SetCurSel(1) ;   
	shBuildingPosition_ctrl.SetCurSel(0);
	//shBuildingTypes_ctrl.SetCurSel(0);
//	shAppLightTypes_ctrl.SetCurSel(0);
	//shVLATypes_ctrl.SetCurSel(0);   

	//m_OK.SubclassDlgItem(IDC_OK, this);
	m_OK.SetIcon(IDI_Go);
	m_OK.DrawTransparent();

 
	//通用机场初始化
	((CButton*)GetDlgItem(IDC_CHECK_AIRPORT))->SetCheck(FALSE);
	GetDlgItem(IDC_shRunwayID)->EnableWindow(FALSE);
	GetDlgItem(IDC_shRunwayWidth)->EnableWindow(FALSE);
	GetDlgItem(IDC_dLatitude)->EnableWindow(FALSE);	
	GetDlgItem(IDC_shRunwayLength)->EnableWindow(FALSE);
	GetDlgItem(IDC_dLongtitude)->EnableWindow(FALSE);	
	GetDlgItem(IDC_shBuildingPosition)->EnableWindow(FALSE);
	GetDlgItem(IDC_fAltitude)->EnableWindow(FALSE);	
	GetDlgItem(IDC_cICAO)->EnableWindow(FALSE);

	GetDlgItem(IDC_STATIC1)->EnableWindow(FALSE);
	GetDlgItem(IDC_STATIC2)->EnableWindow(FALSE);
	GetDlgItem(IDC_STATIC5)->EnableWindow(FALSE);
	GetDlgItem(IDC_STATIC3)->EnableWindow(FALSE);
	GetDlgItem(IDC_STATIC6)->EnableWindow(FALSE);
	GetDlgItem(IDC_STATIC4)->EnableWindow(FALSE);
	GetDlgItem(IDC_STATIC7)->EnableWindow(FALSE);
	GetDlgItem(IDC_STATIC8)->EnableWindow(FALSE);
	GetDlgItem(IDC_STATIC9)->EnableWindow(FALSE);


	return TRUE;  // return TRUE unless you set the focus to a control
	// 异常: OCX 属性页应返回 FALSE
}

void CGeneralapt::OnEnKillfocuscICAO()
{
	// TODO: 在此添加控件通知处理程序代码
	flag=0;
}

void CGeneralapt::OnEnSetfocuscICAO()
{
	// TODO: 在此添加控件通知处理程序代码
	flag=1;
}

HBRUSH CGeneralapt::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
	HBRUSH hbr = CDialog::OnCtlColor(pDC, pWnd, nCtlColor);

	// TODO:  在此更改 DC 的任何属性
	switch(nCtlColor)
	{
	case CTLCOLOR_STATIC:
		pDC ->SetBkMode(TRANSPARENT);
		pDC->SetTextColor(RGB(0,0,0)); 
		return (HBRUSH)::GetStockObject(HOLLOW_BRUSH); 
	/*case CTLCOLOR_EDIT:
		//	pDC ->SetBkColor(RGB(0,0,0));
		pDC->SetTextColor(RGB(0,0,0)); 
	//	hbr   =   ::CreateSolidBrush(RGB(188,197,230));   
		//pDC-> SetBkMode(TRANSPARENT); 
		return (HBRUSH)::GetStockObject(HOLLOW_BRUSH);*/
	}


	return hbr;
}

BOOL CGeneralapt::OnEraseBkgnd(CDC* pDC)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值

	return TRUE; //CDialog::OnEraseBkgnd(pDC);
}

void CGeneralapt::OnPaint()
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


LRESULT CGeneralapt::OnListSet(WPARAM wParam, LPARAM lParam)
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

void CGeneralapt::OnBnClickedCheckAirport()
{
	// TODO: 在此添加控件通知处理程序代码
	if(((CButton*)GetDlgItem(IDC_CHECK_AIRPORT))-> GetCheck())
	{
		GetDlgItem(IDC_shRunwayID)->EnableWindow(TRUE);
		GetDlgItem(IDC_shRunwayWidth)->EnableWindow(TRUE);
		GetDlgItem(IDC_dLatitude)->EnableWindow(TRUE);	
		GetDlgItem(IDC_shRunwayLength)->EnableWindow(TRUE);
		GetDlgItem(IDC_dLongtitude)->EnableWindow(TRUE);	
		GetDlgItem(IDC_shBuildingPosition)->EnableWindow(TRUE);
		GetDlgItem(IDC_fAltitude)->EnableWindow(TRUE);	
		GetDlgItem(IDC_cICAO)->EnableWindow(TRUE);

		GetDlgItem(IDC_STATIC1)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC2)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC5)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC3)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC6)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC4)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC7)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC8)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC9)->EnableWindow(TRUE);
	}
	else
	{
		GetDlgItem(IDC_shRunwayID)->EnableWindow(FALSE);
		GetDlgItem(IDC_shRunwayWidth)->EnableWindow(FALSE);
		GetDlgItem(IDC_dLatitude)->EnableWindow(FALSE);	
		GetDlgItem(IDC_shRunwayLength)->EnableWindow(FALSE);
		GetDlgItem(IDC_dLongtitude)->EnableWindow(FALSE);	
		GetDlgItem(IDC_shBuildingPosition)->EnableWindow(FALSE);
		GetDlgItem(IDC_fAltitude)->EnableWindow(FALSE);	
		GetDlgItem(IDC_cICAO)->EnableWindow(FALSE);

		GetDlgItem(IDC_STATIC1)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC2)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC5)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC3)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC6)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC4)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC7)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC8)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC9)->EnableWindow(FALSE);
	}
}

