#pragma once
#include "afxwin.h"
#include "BtnST.h"
#include "BkDialogST.h"
#include "StaticST.h"

// CGeneralapt 对话框

class CGeneralapt : public CBkDialogST
{
	DECLARE_DYNAMIC(CGeneralapt)

public:
	CGeneralapt(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~CGeneralapt();

// 对话框数据
	enum { IDD = IDD_GENERALAPT };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	short shRunwayID;
	double dLatitude;
	double dLongtitude;
	float fAltitude;
	float fHeading;
	float fVLAAngles;
	float fRunwaySlope;
	CComboBox shRunwayWidth_ctrl;
	CComboBox shRunwayLength_ctrl;
	CComboBox shRunwayCondition_ctrl;
	CComboBox shBuildingTypes_ctrl;
	CComboBox shAppLightTypes_ctrl;
	CComboBox shBuildingPosition_ctrl;
	CComboBox shVLATypes_ctrl;
	CString cICAO;
	CString cLocation;
	afx_msg void OnBnClickedOK();
	virtual BOOL OnInitDialog();
	CEdit cICAO_ctrl;

	//ListBox Controller
	CListBox m_listbox1_ctrl;
	afx_msg LRESULT OnListSet(WPARAM wParam, LPARAM lParam);
	CString m_listbox1;
	
	afx_msg void OnEnKillfocuscICAO();
	int flag;
	afx_msg void OnEnSetfocuscICAO();
	CNumEdit m_cEdit[7];
	afx_msg HBRUSH OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	CButtonST m_OK;
	CButtonST m_BUTTON2;
	afx_msg void OnPaint();
	bool m_createApt;
	afx_msg void OnBnClickedCheckAirport();
};
