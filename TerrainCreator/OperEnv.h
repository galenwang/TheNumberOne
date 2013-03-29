#pragma once
#include <vector>
#include "afxwin.h"
#include "BtnST.h"
#include "BkDialogST.h"
#include "StaticST.h"
#include "XPStyleButtonST.h"
#include "ThemeHelperST.h"
using namespace std;

class CCharHandle
{
public:
	CCharHandle();
	~CCharHandle();
public:
	inline char* UnicodeToAnsi( const wchar_t* szStr )
	{
		int nLen = WideCharToMultiByte( CP_ACP, 0, szStr, -1, NULL, 0, NULL, NULL );
		if (nLen == 0)
		{
		   return NULL;
		}
		char* pResult = new char[nLen];
		WideCharToMultiByte( CP_ACP, 0, szStr, -1, pResult, nLen, NULL, NULL );
		resV.push_back(pResult);
		return pResult;
	}
private:
     vector<char*> resV;
};


// COperEnv 对话框

class COperEnv : public CBkDialogST
{
	DECLARE_DYNAMIC(COperEnv)

public:
	HACCEL m_hAccel;
	COperEnv(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~COperEnv();

// 对话框数据
	enum { IDD = IDD_OPERENV };
	CTreeCtrl	m_tree;
	CString	m_strItem;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	
	CString m_Lod;
	CButtonST m_PiLiang_Button;
	CButtonST m_Refresh_Button;
	CButtonST m_Reset_Button;
	afx_msg void OnBnClickedPiLiang_Button();
	afx_msg void OnBnClickedRefresh_Button();
	afx_msg void OnBnClickedReset_Button();
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
	afx_msg void OnTimer(UINT_PTR nIDEvent);
};
