// ProgressDialog.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "ProgressDialog.h"
#include "interface.h"
#include "OperEnv.h"
#include "dijingzhizuoDlg.h"



// CProgressDialog �Ի���

IMPLEMENT_DYNAMIC(CProgressDialog, CDialog)

CProgressDialog::CProgressDialog(CWnd* pParent /*=NULL*/)
	: CDialog(CProgressDialog::IDD, pParent)
{
	
theParent=NULL;
	
}

CProgressDialog::~CProgressDialog()
{
	
}

void CProgressDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CProgressDialog, CDialog)
	ON_WM_PAINT()
	ON_BN_CLICKED(IDOK, &CProgressDialog::OnBnClickedOk)
	ON_BN_CLICKED(IDC_BUTTON1, &CProgressDialog::OnBnClickedButton1)
END_MESSAGE_MAP()


// CProgressDialog ��Ϣ�������

void CProgressDialog::OnPaint()
{
	CPaintDC pdc(this); // device context for painting
	// TODO: �ڴ˴������Ϣ����������
	CRect rc;
	this->GetClientRect(&rc);
	      
	CBrush brush,brush1;
	if(mark)
	{
		brush.CreateSolidBrush(RGB(255,0,0));
		pdc.SelectObject(brush); 
		for (int i=0 ;i<m_XSum;i++)
			for(int j=0;j<m_YSum;j++)
			{
				pdc.Rectangle((rc.left+i*12),(rc.top+j*12),(rc.left+12*(i+1)),(rc.top+12*(j+1))); //���ƾ���
			}

	}else
	{
		brush.CreateSolidBrush(RGB(0,0,255)); 
		pdc.SelectObject(brush); 

		pdc.Rectangle((rc.left+m_X*12),(rc.top+m_Y*12),(rc.left+12*(m_X+1)),(rc.top+12*(m_Y+1))); //���ƾ���
		Invalidate(false);
		
	}
}

void CProgressDialog::OnBnClickedOk()
{
	// TODO: �ڴ���ӿؼ�֪ͨ��������
	

	OnOK();
}
char test[]="add a list";
void CProgressDialog::OnBnClickedButton1()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	mark=false;
	//OnPaint();
	SendMessage(WM_PAINT);
	m_X++;
	m_Y++;
	COperEnv ctrl;

	if(theParent)
		theParent->SendMessage(WM_OPERENV,0,(LPARAM)test);
//	CdijingzhizuoDlg *dlg;//=(CdijingzhizuoDlg*)m_tabPages[0];

//OperEnv* EnvClouds = (COperEnv*) (((CdijingzhizuoDlg*) GetParent()->GetParent())->m_tabPages[0]);
//	HWND hWnd = ::FindWindow("COperEnv",NULL);

 	
// 	ctrl.m_listbox1_ctrl(WM_OPERENV,0,(LPARAM)"add a list");
 	
 
 	
	
	//COperEnv * ctrl=( COperEnv * ) CdijingzhizuoDlg->m_tabPages[0];
//	CEnvClouds* EnvClouds = (CEnvClouds*) (((CEnvWeather*) GetParent()->GetParent())->m_tabPages[0]);

}

BOOL CProgressDialog::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  �ڴ���Ӷ���ĳ�ʼ��
	
	 m_XSum=20;
	 m_YSum=10;
	 m_X=0;
	 m_Y=0;

	return TRUE;  // return TRUE unless you set the focus to a control
	// �쳣: OCX ����ҳӦ���� FALSE
}
void CProgressDialog::SetParentPointer(COperEnv* pParent)
{
	theParent=pParent;
}