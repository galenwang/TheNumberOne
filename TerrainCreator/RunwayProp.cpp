// RunwayProp.cpp : 实现文件
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "RunwayProp.h"
#include "interface.h"
#include "OperEnv.h"

// CRunwayProp 对话框
extern char workDirectory[];

IMPLEMENT_DYNAMIC(CRunwayProp, CBkDialogST)

CRunwayProp::CRunwayProp(CWnd* pParent /*=NULL*/)
	: CBkDialogST(CRunwayProp::IDD, pParent)
	, m_Lenght(0)
	, m_Width(0)
	, m_Heading(0)
	, m_Longitude(0)
	, m_Latitude(0)
	, m_DispThr1(0)
	, m_DispThr2(0)
	, m_Stopway1(0)
	, m_Stopway2(0)
	, m_RunwayLight1(_T(""))
	, m_RunwayLight2(_T(""))
	, m_ApproachLight1(_T(""))
	, m_ApproachLight2(_T(""))
	, m_GlideslopeLight1(_T(""))
	, m_GlideslopeLight2(_T(""))
	, m_SurfaceType(_T(""))
	, m_Marking(_T(""))
	, m_Should(_T(""))
	, m_Roughness(0)
	, m_Signage(FALSE)
	, m_RunwayAltitude(0)
{

}

CRunwayProp::~CRunwayProp()
{
}

void CRunwayProp::DoDataExchange(CDataExchange* pDX)
{
	CBkDialogST::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_Lenght_EDIT, m_Lenght);
	DDX_Text(pDX, IDC_Width_EDIT, m_Width);
	DDX_Text(pDX, IDC_Heading_EDIT, m_Heading);
	DDV_MinMaxFloat(pDX, m_Heading, 0, 360);
	DDX_Text(pDX, IDC_Longitude_EDIT, m_Longitude);
	DDX_Text(pDX, IDC_Latitude_EDIT, m_Latitude);
	DDX_Text(pDX, IDC_DispThr1_EDIT, m_DispThr1);
	DDX_Text(pDX, IDC_DispThr2_EDIT, m_DispThr2);
	DDX_Text(pDX, IDC_Stopway1_EDIT, m_Stopway1);
	DDX_Text(pDX, IDC_Stopway2_EDIT, m_Stopway2);
	DDX_CBString(pDX, IDC_RunwayLight1_COMB, m_RunwayLight1);
	DDX_CBString(pDX, IDC_RunwayLight2_COMB, m_RunwayLight2);
	DDX_CBString(pDX, IDC_ApproachLight1_COMB, m_ApproachLight1);
	DDX_CBString(pDX, IDC_ApproachLight2_COMB, m_ApproachLight2);
	DDX_CBString(pDX, IDC_GlideslopeLight1_COMB, m_GlideslopeLight1);
	DDX_CBString(pDX, IDC_GlideslopeLight2_COMB, m_GlideslopeLight2);
	DDX_CBString(pDX, IDC_SurfaceType_COMB, m_SurfaceType);
	DDX_CBString(pDX, IDC_COMBO9, m_Marking);
	DDX_CBString(pDX, IDC_Should_COMB, m_Should);
	DDX_Text(pDX, IDC_EDIT10, m_Roughness);
	DDV_MinMaxFloat(pDX, m_Roughness, 0, 1);
	DDX_Check(pDX, IDC_CHECK1, m_Signage);
	DDV_MinMaxFloat(pDX, m_Longitude, -180, 180);
	DDV_MinMaxFloat(pDX, m_Latitude, -90, 90);

	DDX_Control(pDX, IDC_RunwayLight1_COMB, m_RunwayLight1_ctrl);
	DDX_CBString(pDX, IDC_RunwayLight1_COMB, m_RunwayLight1);
	DDX_Control(pDX, IDC_ApproachLight1_COMB, m_ApproachLight1_ctrl);
	DDX_Control(pDX, IDC_GlideslopeLight1_COMB, m_GlideslopeLight1_ctrl);
	DDX_Control(pDX, IDC_RunwayLight2_COMB, m_RunwayLight2_ctrl);
	DDX_Control(pDX, IDC_ApproachLight2_COMB, m_ApproachLight2_ctrl);
	DDX_Control(pDX, IDC_GlideslopeLight2_COMB, m_GlideslopeLight2_ctrl);
	DDX_Control(pDX, IDC_SurfaceType_COMB, m_SurfaceType_ctrl);
	DDX_Text(pDX, IDC_RunwayAlt_EDIT, m_RunwayAltitude);
	DDX_Control(pDX, IDC_Should_COMB, m_Should_ctrl);

//	DDX_Control(pDX, IDC_CHECK1, m_CHECK1);
	DDX_Control(pDX, IDC_OK, m_OK);
	DDX_Control(pDX, IDC_cancel, m_cancel);
}


BEGIN_MESSAGE_MAP(CRunwayProp, CBkDialogST)
	ON_BN_CLICKED(IDC_OK, &CRunwayProp::OnBnClickedOK)
	ON_WM_CTLCOLOR()
	ON_WM_ERASEBKGND()
	ON_WM_PAINT()
END_MESSAGE_MAP()


// CRunwayProp 消息处理程序

void CRunwayProp::OnBnClickedOK()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(true);
	tagRunwayProp RunwayProp;

	RunwayProp.Lenght=m_Lenght;
	RunwayProp.Width=m_Width;
	RunwayProp.Heading=m_Heading;
	RunwayProp.Longitude=m_Longitude;
	RunwayProp.Latitude=m_Latitude;
	RunwayProp.DispThr1=m_DispThr1;
	RunwayProp.DispThr2=m_DispThr2;
	RunwayProp.Stopway1=m_Stopway1;
	RunwayProp.Stopway2=m_Stopway2;
#if 1	//_DEBUG	
	strcpy_s(RunwayProp.RunwayLight1,(LPCSTR)m_RunwayLight1.GetBuffer());
	strcpy_s(RunwayProp.RunwayLight2,(LPCSTR)m_RunwayLight2.GetBuffer());
	strcpy_s(RunwayProp.ApproachLight1,(LPCSTR)m_ApproachLight1.GetBuffer());
	strcpy_s(RunwayProp.ApproachLight2,(LPCSTR)m_ApproachLight2.GetBuffer());
	strcpy_s(RunwayProp.GlideslopeLight2,(LPCSTR)m_GlideslopeLight2.GetBuffer());
	strcpy_s(RunwayProp.GlideslopeLight1,(LPCSTR)m_GlideslopeLight1.GetBuffer());
	strcpy_s(RunwayProp.SurfaceType,(LPCSTR)m_SurfaceType.GetBuffer());
	strcpy_s(RunwayProp.Marking,(LPCSTR)m_Marking.GetBuffer());
	strcpy_s(RunwayProp.Should,(LPCSTR)m_Should.GetBuffer());
#else
	CCharHandle InterStr;
	strcpy_s(RunwayProp.RunwayLight1,InterStr.UnicodeToAnsi(m_RunwayLight1.GetBuffer()));
	strcpy_s(RunwayProp.RunwayLight2,InterStr.UnicodeToAnsi(m_RunwayLight2.GetBuffer()));
	strcpy_s(RunwayProp.ApproachLight1,InterStr.UnicodeToAnsi(m_ApproachLight1.GetBuffer()));
	strcpy_s(RunwayProp.ApproachLight2,InterStr.UnicodeToAnsi(m_ApproachLight2.GetBuffer()));
	strcpy_s(RunwayProp.GlideslopeLight2,InterStr.UnicodeToAnsi(m_GlideslopeLight2.GetBuffer()));
	strcpy_s(RunwayProp.GlideslopeLight1,InterStr.UnicodeToAnsi(m_GlideslopeLight1.GetBuffer()));
	strcpy_s(RunwayProp.SurfaceType,InterStr.UnicodeToAnsi(m_SurfaceType.GetBuffer()));
	strcpy_s(RunwayProp.Marking,InterStr.UnicodeToAnsi(m_Marking.GetBuffer()));
	strcpy_s(RunwayProp.Should,InterStr.UnicodeToAnsi(m_Should.GetBuffer()));
#endif
	
	RunwayProp.Roughness=m_Roughness;
	RunwayProp.Signage=m_Signage; 
	
#if 0			//专门输出给通用机场用；	
	short 	shGenericAirportState;//		通用机场状态		[0，1]
	short 	shRunwayID;//              	通用机场跑道编号		[01，36]
	short 	shRunwayWidth;//           	通用机场跑道宽度		[0，2]
	short 	shRunwayCondition;//       	通用机场跑道条件		[0，3]
	short 	shRunwayLength;//          	通用机场跑道长度		[0，7]
	short 	shBuildingPosition;//      	通用机场建筑物位置		[0，1]
	short 	shBuildingTypes;//         	通用机场建筑物类型		[0，1]
	short 	shAppLightTypes;//         	通用机场跑道进近灯光类型		[0，2]
	short 	shVLATypes;//              	通用机场跑道助降灯光类型		[0，2]
	double	dLatitude;//               	通用机场跑道纬度	deg	[-90.0，90.0]
	double	dLongtitude;//             	通用机场跑道经度	deg	[-180.0，180.0]  
	float 	fAltitude;//               	通用机场跑道高度	ft	[-500.0，50000.0]
	float 	fHeading;//                	通用机场跑道航向角	deg	[0.0，360.0]     
	float 	fVLAAngles;//				通用机场跑道助降灯光角度	deg	[2.4，4.0] 
	float 	fRunwaySlope;//            	通用机场跑道倾斜角度	deg	[-0.5，0.5]  


	shGenericAirportState=1;
	fAltitude= m_RunwayAltitude;   
	if (m_Width>=0 && m_Width<=30)
	{
		shRunwayWidth=0;
		
	} 
	else
	if (m_Width>30 && m_Width<=45)
	{
		shRunwayWidth=1;
	} 
	else
		if(m_Width>45 )
	{
		shRunwayWidth=2;
	} 
		
	   //0:Dry/Normal 1:Wet 2:Snow 3:Ice
		//Asphalt;Concrete;Turf;Dirt;Gravel;Water;Dry lakebed;Other;
	switch(m_SurfaceType_ctrl.GetCurSel())
	{
	case 0:
	case 1:
	case 2:
	case 3:
	case 4:
	case 6:
		shRunwayCondition=0;
		break;
	case 5:
		shRunwayCondition=1;
		break;
	case 7:
		shRunwayCondition=2;
		break;

	}
	//0:0m 1:1500m 2:2000m 3:2500m 4:3000m 5:3500m 6:4000m 7:4500m
	if (m_Lenght==0)
		shRunwayLength=0;
	else if (m_Lenght>0&&m_Lenght<=1500)
			shRunwayLength=1;
	else if (m_Lenght>1500&&m_Lenght<=2000)
		shRunwayLength=2;
	else if (m_Lenght>2000&&m_Lenght<=2500)
		shRunwayLength=3;
	else if (m_Lenght>2500&&m_Lenght<=3000)
		shRunwayLength=4;
	else if (m_Lenght>3000&&m_Lenght<=3500)
		shRunwayLength=5;
	else if (m_Lenght>3500&&m_Lenght<=4000)
		shRunwayLength=6;
	else if (m_Lenght>4000)
		shRunwayLength=7;
	
	
 
	shBuildingPosition=1;//
	shBuildingTypes=1;//   
	shAppLightTypes=1;//   
	shVLATypes=2;//        
	dLatitude=m_Latitude;//         
	dLongtitude=m_Longitude;//       
	fHeading=m_Heading;//          
	fVLAAngles=3.000;//		
	fRunwaySlope=0.3;//
	shRunwayID = (int)(m_Heading / 10.0) + 1; 
	


#else

/*
1   1572 1 1 ZLXY XIANYANG
10  34.446863  108.751861 05x  48.75   9849    0.0000  197.0197  148 344354 02 0 3 0.25 0 0300.0300
*/
	
	char rwyStr[2][256], lights[7];
	int  i;
	
	memset(rwyStr, 0, sizeof(rwyStr));
	sprintf(rwyStr[0], "1   %d  1  1  GENZ GENARALAIRPORT", m_RunwayAltitude);
	sprintf(rwyStr[1], "10  %f  %f  ", m_Latitude, m_Longitude);
	
	int runwayID = (int)(m_Heading / 10.0) + 1; 
	sprintf(rwyStr[1], "%s%02dx  %f  %d  ", rwyStr[1], runwayID, m_Heading, m_Width);
	sprintf(rwyStr[1], "%s%04d.%04d  04d.%04d  %f  ", rwyStr[1], m_DispThr1, m_DispThr2, m_Stopway1, m_Stopway2, m_Width);
	
	
	//////////////////////////////////////////////////////////////////////////
	//跑道灯_1
	//	None;Edge lighting;REIL;Center lighting;TZ lighting;
	//	rwylt1>= 2   Has edge lighting  
	//	rwylt1>= 3   Has REIL lighting 
	//	rwylt1>= 4		Has centerline lighting 
	//	rwylt1>= 5	 Has touchdown zone lighting
	//	rwylt1 == 6 /	taxiway lit, assume centerline too 
	//	rwylt1 == 6  taxiway blue lit 

	switch(m_RunwayLight1_ctrl.GetCurSel())
	{
	case 0:
		i=0;
		break;
	case 1:
	case 2:
	case 3:
	case 4:
		i=m_RunwayLight1_ctrl.GetCurSel()+1;	
		break;
	default:
		i=0;
		break;
	}
	lights[1] = 30 + i;


	//////////////////////////////////////////////////////////////////////////
	//滑行道——1
	//None;VASI;PAPI;SSL PAPI;
	//vasi1 == 2 /* Has VASI
	//vasi1 == 3 /* Has PAPI 
	switch(m_GlideslopeLight1_ctrl.GetCurSel())
	{
	case 1:
		i=2;
		break;
	case 2:
		i=3;
		break;
	default:
		i=1;
		break;
	}
	lights[0] = 30 + i;

	//////////////////////////////////////////////////////////////////////////
	//进场灯——1
	//None;ALS;ALSF I;ALSF II;CAL;     CAL II;LDIN;MALS;MALSF;NSTD;    MALSR;MIL OVRN;ODALS;RAIL;SALS;   SALSF;SSALF;SSALR;SSALS;
	//	app1 == -1 /* MALS not supported by data base 
	//	app1 == 2 /* SSALS 
	//	app1 == 3 /* SALSF 
	//	app1 == 4 /* ALSF-I
	//	app1 == 5 /* ALSF-II 
	//	app1 == 6 /* ODALS Omni-directional approach light system
	//	app1 == 7 || app1 == 8 /* Calvert 1, 2, and 3
	switch(m_ApproachLight1_ctrl.GetCurSel())
	{
	case 18:
		i=2;
		break;
	case 15:
		i=3;
		break;
	case 2:
		i=4;
		break;
	case 3:
		i=5;
		break;
	case 13:
		i=6;
		break;
	case 4:
	case 5:
		i=7;
	default:
		i=0;
		break;
	}
	lights[2] = 30 + i;
	
	//////////////////////////////////////////////////////////////////////////
	///跑道灯——2
	//	None;Edge lighting;REIL;Center lighting;TZ lighting;
	//	rwylt1>= 2 /* Has edge lighting 
	//	rwylt1>= 3 /* Has REIL lighting 
	//	rwylt1>= 4 /* Has centerline lighting
	//	rwylt1>= 5 /* Has touchdown zone lighting
	//	rwylt1 == 6 /* taxiway lit, assume centerline too
	//	rwylt1 == 6 /* taxiway blue lit 
	switch(m_RunwayLight2_ctrl.GetCurSel())
	{
	case 0:
		i=0;
		break;
	case 1:
	case 2:
	case 3:
	case 4:
		i=m_RunwayLight2_ctrl.GetCurSel()+1;	
		break;
	}
	lights[4] = 30 + i;

	//////////////////////////////////////////////////////////////////////////
	//滑行道——2
	//None;VASI;PAPI;SSL PAPI;
	//vasi1 == 2 /* Has VASI
	//vasi1 == 3 /* Has PAPI 
	switch(m_GlideslopeLight2_ctrl.GetCurSel())
	{
	case 1:
		i=2;
		break;
	case 2:
		i=3;
		break;
	default:
		i=1;
		break;
	}
	lights[3] = 30 + i;
	
	//////////////////////////////////////////////////////////////////////////
	//进场灯——2
	//None;ALS;ALSF I;ALSF II;CAL;     CAL II;LDIN;MALS;MALSF;NSTD;    MALSR;MIL OVRN;ODALS;RAIL;SALS;   SALSF;SSALF;SSALR;SSALS;
	//app1 == -1 /* MALS not supported by data base
	//app1 == 2 /* SSALS 
	//app1 == 3 /* SALSF 
	//app1 == 4 /* ALSF-I 
	//app1 == 5 /* ALSF-II 
	//app1 == 6 /* ODALS Omni-directional approach light system
	//app1 == 7 || app1 == 8 /* Calvert 1, 2, and 3
	switch(m_ApproachLight2_ctrl.GetCurSel())
	{
	case 18:
		i=2;
		break;
	case 15:
		i=3;
		break;
	case 2:
		i=4;
		break;
	case 3:
		i=5;
		break;
	case 13:
		i=6;
		break;
	default:
		i=0;
		break;
	}
	lights[5] = 30 + i;
	lights[6] = 0;
	sprintf(rwyStr[1], "%s%s  ", rwyStr[1], lights);
	

	//////////////////////////////////////////////////////////////////////////
	//跑道材质
	//Asphalt;Concrete;Turf;Dirt;Gravel;Water;Dry lakebed;Other;
	//1 /* Asphalt 
	//2 /* Concrete 
	//3 /* Turf/Grass 
	//4 /* Dirt 
	//5 /* Gravel 
	//12 /* Dry Lakebed
	//13 /* Water runway (buoy's?) 
	switch(m_SurfaceType_ctrl.GetCurSel())
	{

	case 0:
	case 1:
	case 2:
	case 3:
	case 4:
		i=m_SurfaceType_ctrl.GetCurSel()+1;
		break;
	case 5:
		i=13;
		break;
	case 6:
		i=12;
		break;
	default:
		i=-1;
		break;
	}
	sprintf(rwyStr[1], "%s%02d  ", rwyStr[1], i);

	/////////////////////////////////////////////////////////////////////
	//路肩
	switch(m_Should_ctrl.GetCurSel())
	{
	case 0:
	case 1:
	case 2:
		i=m_Should_ctrl.GetCurSel()+1;
	}
	sprintf(rwyStr[1], "%s%d  ", rwyStr[1], i);

	//////////////////////////////////////////////////////////////////////////
	//标记
	//None;Visual;Non-precision;Precision;Helipad;
	//0 /* No known markings, lets assume Visual 
	//1 /* Visual
	//2 /* Non-precision 
	//3 /* Precision
	switch(m_SurfaceType_ctrl.GetCurSel())
	{
	case 0:
	case 1:
	case 2:
	case 3:
		i=m_SurfaceType_ctrl.GetCurSel()+1;
		break;
	case 4:
		i=0;
		break;	
	}
	sprintf(rwyStr[1], "%s%d  ", rwyStr[1], i);
	sprintf(rwyStr[1], "%s%f  0.00  0300.0300", rwyStr[1], m_Roughness);

	/* 写dat文件 */
	FILE *fi;
	char fn[256];
	sprintf(fn, "%s\\generalapt\\GENZ.dat", workDirectory);
	fi = fopen(fn, "rt");
	if(fi == NULL)
		return;
	fprintf(fi, "%s\n", rwyStr[0]);
	fprintf(fi, "%s", 	rwyStr[1]);
	fclose(fi);
	
	/* 编写命令行调用genapt*/


#endif

	
}

BOOL CRunwayProp::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  在此添加额外的初始化
	m_cEdit[0].SubclassDlgItem(IDC_Lenght_EDIT,this);
	m_cEdit[1].SubclassDlgItem(IDC_Width_EDIT,this);
	m_cEdit[2].SubclassDlgItem(IDC_Heading_EDIT,this);
	m_cEdit[3].SubclassDlgItem(IDC_Longitude_EDIT,this);
	m_cEdit[4].SubclassDlgItem(IDC_Latitude_EDIT,this);
	m_cEdit[5].SubclassDlgItem(IDC_DispThr1_EDIT,this);
	m_cEdit[6].SubclassDlgItem(IDC_Stopway1_EDIT,this);
	m_cEdit[7].SubclassDlgItem(IDC_DispThr2_EDIT,this);
	m_cEdit[8].SubclassDlgItem(IDC_Stopway2_EDIT,this);
	m_cEdit[9].SubclassDlgItem(IDC_EDIT10,this);

//	m_CHECK1.SetIcon(IDI_LedOn, IDI_LedOff);
//	m_CHECK1.DrawTransparent();

	m_OK.SetThemeHelper(&m_Theme);				    
	m_OK.SetIcon(IDI_ok, (int)BTNST_AUTO_GRAY);
	m_OK.OffsetColor(CButtonST::BTNST_COLOR_BK_IN, 30);
	m_OK.SetTooltipText(_T("Close the application"));
	//m_OK.DrawTransparent();

	m_cancel.SetThemeHelper(&m_Theme);
	m_cancel.SetIcon(IDI_Cancel, (int)BTNST_AUTO_GRAY);
	m_cancel.OffsetColor(CButtonST::BTNST_COLOR_BK_IN, 30);
	m_cancel.SetTooltipText(_T("Step progress bars"));
	//m_cancel.DrawTransparent();

	if (m_Theme.IsAppThemed() == FALSE)
	{
		m_OK.DrawTransparent();
		m_cancel.DrawTransparent();
	} 


	return TRUE;  // return TRUE unless you set the focus to a control
	// 异常: OCX 属性页应返回 FALSE
}

HBRUSH CRunwayProp::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
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
		return (HBRUSH)::GetStockObject(HOLLOW_BRUSH);*/
	}


	return hbr;
}

BOOL CRunwayProp::OnEraseBkgnd(CDC* pDC)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值

	return TRUE;// CDialog::OnEraseBkgnd(pDC);
}

void CRunwayProp::OnPaint()
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
