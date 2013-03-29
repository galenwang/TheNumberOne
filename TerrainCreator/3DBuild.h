#pragma once
#include "afxwin.h"
#include "BtnST.h"
#include "BkDialogST.h"
#include "StaticST.h"
// C3DBuild �Ի���

class C3DBuild : public CDialog
{
	DECLARE_DYNAMIC(C3DBuild)

public:
	C3DBuild(CWnd* pParent = NULL);   // ��׼���캯��
	virtual ~C3DBuild();

// �Ի�������
	enum { IDD = IDD_3DBUILD };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

	DECLARE_MESSAGE_MAP()
public:
	afx_msg HBRUSH OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor);
	afx_msg void OnPaint();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	virtual BOOL OnInitDialog();
	CButtonST m_Batch_Button;
	afx_msg void OnBnClickedBatch();

	//ListBox Controller
	CListBox m_listbox1_ctrl;
	afx_msg LRESULT OnListSet(WPARAM wParam, LPARAM lParam);
	CString m_listbox1;
	
	BOOL m_rebuild;
};
