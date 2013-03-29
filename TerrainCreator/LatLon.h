#pragma once
#include "afxwin.h"
#include "BtnST.h"
#include "BkDialogST.h"
#include "StaticST.h"
#include "XPStyleButtonST.h"
#include "ThemeHelperST.h"

// CLatLon 对话框

class CLatLon : public CBkDialogST
{
	DECLARE_DYNAMIC(CLatLon)

public:
	HACCEL m_hAccel;
	CLatLon(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~CLatLon();

// 对话框数据
	enum { IDD = IDD_LATLON };
	CTreeCtrl	m_tree;
	CString	m_strItem;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	float roadWide;
	float woodCoverage;
	float lightDistance;
	BOOL m_lampType;

	CButtonST m_Create_Button;
	CButtonST m_Refresh_Button;
	CButtonST m_Reset_Button;
	afx_msg void OnBnClickedCreateButton();
	afx_msg void OnBnClickedRefresh_Button();
	afx_msg void OnBnClickedReset_Button();
	afx_msg void OnBnClickedLamp();

	virtual BOOL OnInitDialog();

	void initTree();
	afx_msg void OnTvnSelchangedTree(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnTvnEndlabeleditTree(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnDblClkTree(NMHDR *pNMHDR, LRESULT *pResult);

	CListBox m_listbox1_ctrl;
	afx_msg LRESULT OnListSet(WPARAM wParam, LPARAM lParam);
	CString m_listbox1;	
	afx_msg LRESULT OnUpdateEdit(WPARAM wParam, LPARAM lParam);

	afx_msg HBRUSH OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor);
	afx_msg void OnPaint();
	void UpdateParameters();

};
