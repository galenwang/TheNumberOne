// dijingzhizuoDlg.cpp : ʵ���ļ�
//

#include <stdio.h>
#include "stdafx.h"
#include "dijingzhizuo.h"
#include "dijingzhizuoDlg.h"
#include "Splash.h" 
#include "travelmap.h" 

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

extern int ifFoundBasePath;
extern FILE *fLog;

// ����Ӧ�ó��򡰹��ڡ��˵���� CAboutDlg �Ի���

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// �Ի�������
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

// ʵ��
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	
}



void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CdijingzhizuoDlg �Ի���




CdijingzhizuoDlg::CdijingzhizuoDlg(CWnd* pParent /*=NULL*/)
	: CBkDialogST(CdijingzhizuoDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}
CdijingzhizuoDlg::~CdijingzhizuoDlg()
{
	/* �ر�Log�ļ���*/
	fclose(fLog);
	
	for (int i=0;i<3;i++)
	{
		delete m_tabPages[i];
	}
	delete pImageList;
	 
}

void CdijingzhizuoDlg::DoDataExchange(CDataExchange* pDX)
{
	CBkDialogST::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_TAB1, m_Tab);
	DDX_Control(pDX, IDC_STATIC_1, m_openenv_pic);
	DDX_Control(pDX, IDC_STATIC_6, m_runway_pic);
//	DDX_Control(pDX, IDC_STATIC_2, m_lonlat_pic);
	DDX_Control(pDX, IDC_STATIC_4, m_generalapt_pic);
//	DDX_Control(pDX, IDC_STATIC_5, m_3dBuild_pic);
}

BEGIN_MESSAGE_MAP(CdijingzhizuoDlg, CBkDialogST)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_NOTIFY(TCN_SELCHANGE, IDC_TAB1, &CdijingzhizuoDlg::OnTcnSelchangeTab1)
	ON_WM_CREATE()
	ON_WM_TIMER()
//	ON_STN_CLICKED(IDC_STATIC_1, &CdijingzhizuoDlg::OnStnClickedStatic1)
END_MESSAGE_MAP()


// CdijingzhizuoDlg ��Ϣ�������

BOOL CdijingzhizuoDlg::OnInitDialog()
{
	CBkDialogST::OnInitDialog();
	SetBitmap(IDB_BITMAP1);

	// ��������...���˵�����ӵ�ϵͳ�˵��С�

	// IDM_ABOUTBOX ������ϵͳ���Χ�ڡ�
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// ���ô˶Ի����ͼ�ꡣ��Ӧ�ó��������ڲ��ǶԻ���ʱ����ܽ��Զ�
	//  ִ�д˲���
	SetIcon(m_hIcon, TRUE);			// ���ô�ͼ��
	SetIcon(m_hIcon, FALSE);		// ����Сͼ��

	// TODO: �ڴ���Ӷ���ĳ�ʼ������
	m_Tab.SetItemSize(CSize (16,160));
	m_Tab.SetMinTabWidth(16);
	pImageList = new CImageList();
	pImageList->Create(32, 32, ILC_COLOR32, 3, 3);
	CBitmap bmp;


	m_Tab.SetImageList(pImageList);   //m_tabΪCTabCtrl���� 

	// 	�ڶ������趨����ͼ�� 
	UINT   nMask   =  TCIF_IMAGE ;  //  and  TCIF_TEXT ;
	int   nImageIndex; 
	int   nItemIndex=0; 
	CString   sTabLabel; 
	
	nImageIndex=0; 
	sTabLabel= "��ǩ1 "; 
	m_Tab.InsertItem(nMask,   nItemIndex,   sTabLabel,   
	nImageIndex,   0L); 
	
	nImageIndex++; 
	nItemIndex++; 
	sTabLabel= "��ǩ2 "; 
	m_Tab.InsertItem(nMask,   nItemIndex,   sTabLabel,   
	nImageIndex,   0L); 
	
	nImageIndex++; 
	nItemIndex++; 
	sTabLabel= "��ǩ3 "; 
	m_Tab.InsertItem(nMask,   nItemIndex,   sTabLabel,   
	nImageIndex,   0L); 
	
	m_tabPages[0]=new COperEnv;
	m_tabPages[1]=new CLatLon;
	m_tabPages[2]=new CGeneralapt;
	m_nNumberOfPages=3;
	
	m_tabPages[0]->Create(IDD_OPERENV, &m_Tab);
	m_tabPages[1]->Create(IDD_LATLON, &m_Tab);
	m_tabPages[2]->Create(IDD_GENERALAPT, &m_Tab);
	
	m_tabPages[0]->ShowWindow(SW_SHOW);
	m_tabPages[1]->ShowWindow(SW_HIDE);
	m_tabPages[2]->ShowWindow(SW_HIDE);
	
	CRect rect;
	m_Tab.GetClientRect(rect);
	rect.InflateRect(-160,20,2,2);
	
	m_tabPages[0]->MoveWindow(&rect) ;
	m_tabPages[1]->MoveWindow(&rect) ;
	m_tabPages[2]->MoveWindow(&rect) ;
	m_tabCurrent=0;
	
	bm1.LoadBitmap(IDB_BITMAP3);
	bm2.LoadBitmap(IDB_BITMAP4);
	
	bm5.LoadBitmap(IDB_BITMAP11);
	bm6.LoadBitmap(IDB_BITMAP12);
	
	bm9.LoadBitmap(IDB_BITMAP9);
	bm10.LoadBitmap(IDB_BITMAP10);
	
	return TRUE;  // ���ǽ��������õ��ؼ������򷵻� TRUE
}

void CdijingzhizuoDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// �����Ի��������С����ť������Ҫ����Ĵ���
//  �����Ƹ�ͼ�ꡣ����ʹ���ĵ�/��ͼģ�͵� MFC Ӧ�ó���
//  �⽫�ɿ���Զ���ɡ�

void CdijingzhizuoDlg::OnPaint()
{
	CRect rect,rc;
	GetWindowRect(&rc);

	GetClientRect(&rect);
	if (IsIconic())
	{
		CPaintDC dc(this); // ���ڻ��Ƶ��豸������

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// ʹͼ���ڹ����������о���
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// ����ͼ��
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
	this->KillTimer(1); 

}

//���û��϶���С������ʱϵͳ���ô˺���ȡ�ù��
//��ʾ��
HCURSOR CdijingzhizuoDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


void CdijingzhizuoDlg::OnTcnSelchangeTab1(NMHDR *pNMHDR, LRESULT *pResult)
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	m_tabPages[m_tabCurrent]->ShowWindow(SW_HIDE);
	m_tabCurrent=m_Tab.GetCurSel();
	m_tabPages[m_tabCurrent]->ShowWindow(SW_SHOW);
	*pResult = 0;

	switch (m_tabCurrent)
	{
	case 0:
		m_openenv_pic.SetBitmap(HBITMAP(bm1));
		m_runway_pic.SetBitmap(HBITMAP(bm9));
		m_generalapt_pic.SetBitmap(HBITMAP(bm5));
		break;
	case 1:
		m_openenv_pic.SetBitmap(HBITMAP(bm2));
		m_runway_pic.SetBitmap(HBITMAP(bm9));
		m_generalapt_pic.SetBitmap(HBITMAP(bm6));
		break;
	case 2:
		m_openenv_pic.SetBitmap(HBITMAP(bm2));
		m_runway_pic.SetBitmap(HBITMAP(bm10));
		m_generalapt_pic.SetBitmap(HBITMAP(bm5));
		break;	
	default:
		break;
	}
}


int CdijingzhizuoDlg::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CDialog::OnCreate(lpCreateStruct) == -1)
		return -1;

	// TODO:  �ڴ������ר�õĴ�������
	CSplashWnd::ShowSplashScreen(this); 
	this->MoveWindow(0,0,0,0); 
	this->SetTimer(1,500,NULL);//

	//Galen��Ӵ���
	if(!initCreateTerrain())
	{
		if(!ifFoundBasePath)
		{
			MessageBox("����·��û����ȫ�ҵ���");
		}else
			exit(1);
	}

	return 0;
}

void CdijingzhizuoDlg::OnTimer(UINT_PTR nIDEvent)
{
	// TODO: �ڴ������Ϣ�����������/�����Ĭ��ֵ
	
	int cxIcon = GetSystemMetrics(SM_CXSCREEN );
	int cyIcon = GetSystemMetrics(SM_CYSCREEN );
	CRect rect;
	int x = (cxIcon - 794  + 1) / 2;
	int y = (cyIcon - 650  + 1) / 2;

	this->MoveWindow(x,y,794,650); 
	CRect rc;
	GetWindowRect(&rc);

	CDialog::OnTimer(nIDEvent);
}

