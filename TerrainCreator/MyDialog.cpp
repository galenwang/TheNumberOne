// MyDialog.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "MyDialog.h"


// CMyDialog �Ի���

IMPLEMENT_DYNAMIC(CMyDialog, CDialog)

CMyDialog::CMyDialog(CWnd* pParent /*=NULL*/)
	//: CDialog(CMyDialog::IDD, pParent)
{

}


CMyDialog::~CMyDialog()
{
}

void CMyDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CMyDialog, CDialog)
	ON_WM_PAINT()
	ON_WM_ERASEBKGND()
END_MESSAGE_MAP()


// CMyDialog ��Ϣ�������

void CMyDialog::OnPaint()
{
	CPaintDC dc(this); // device context for painting
	// TODO: �ڴ˴������Ϣ����������
	// ��Ϊ��ͼ��Ϣ���� CDialog::OnPaint()
}

BOOL CMyDialog::OnEraseBkgnd(CDC* pDC)
{
	// TODO: �ڴ������Ϣ�����������/�����Ĭ��ֵ

	return TRUE;//CDialog::OnEraseBkgnd(pDC);
}
