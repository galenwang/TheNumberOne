#pragma once
#include "afxwin.h"
#include "BtnST.h"
#include "BkDialogST.h"
#include "StaticST.h"

// CManageProfile 对话框

class CManageProfile : public CDialog
{
	DECLARE_DYNAMIC(CManageProfile)

public:
	HACCEL m_hAccel;
	CManageProfile(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~CManageProfile();

// 对话框数据
	enum { IDD = IDD_PROFILE };
	CTreeCtrl	m_tree;
	CString	m_strItem;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()

public:
	afx_msg HBRUSH OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor);
	afx_msg void OnPaint();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	virtual BOOL OnInitDialog();
	CButtonST m_Edit_Button,m_Save_Button,m_Delete_Button,m_Refresh_Button;

	void initTree();
	afx_msg void OnBnClickedEdit();
	afx_msg void OnTvnSelchangedTree(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnTvnEndlabeleditTree(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnDblClkTree(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnBnClickedSave();
	afx_msg void OnBnClickedDelete();
	afx_msg void OnBnClickedRefresh();
};
