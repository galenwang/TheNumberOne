// dijingzhizuoDlg.h : ͷ�ļ�
//

#pragma once
#include "afxcmn.h"
#include "OperEnv.h"
//#include "RunwayProp.h"
#include "LatLon.h"
#include "Generalapt.h"
#include "BkDialogST.h"
#include "MyTabCtrl.h"
#include "afxwin.h"
#include "3DBuild.h"
#include "ManageProfile.h"

// CdijingzhizuoDlg �Ի���
class CdijingzhizuoDlg : public CBkDialogST
{
// ����
public:
	CdijingzhizuoDlg(CWnd* pParent = NULL);	// ��׼���캯��
	~CdijingzhizuoDlg();

// �Ի�������
	enum { IDD = IDD_DIJINGZHIZUO_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV ֧��
	
	

// ʵ��
protected:
	HICON m_hIcon;

	// ���ɵ���Ϣӳ�亯��
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnTcnSelchangeTab1(NMHDR *pNMHDR, LRESULT *pResult);
	CMyTabCtrl m_Tab;
	CDialog *m_tabPages[3]  ;
	int m_tabCurrent;
	int m_nNumberOfPages;
	int initialize;
	CImageList* pImageList;
protected:

public:
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	CStatic m_openenv_pic;
	CBitmap bm1,bm2,bm3,bm4,bm5,bm6,bm7,bm8,bm9,bm10;
//	afx_msg void OnStnClickedStatic1();
	CStatic m_runway_pic;
	CStatic m_lonlat_pic;
	CStatic m_generalapt_pic;
	CStatic m_3dBuild_pic;
};
