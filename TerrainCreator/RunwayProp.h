#pragma once
#include "afxwin.h"
#include "BtnST.h"
#include "BkDialogST.h"
#include "StaticST.h"
#include "XPStyleButtonST.h"
#include "ThemeHelperST.h"

// CRunwayProp 对话框

class CRunwayProp : public CBkDialogST
{
	DECLARE_DYNAMIC(CRunwayProp)

public:
	CRunwayProp(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~CRunwayProp();

// 对话框数据
	enum { IDD = IDD_RUNWAYPROP };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	int m_Lenght;
	int m_Width;
	float m_Heading;
	float m_Longitude;
	float m_Latitude;
	int m_DispThr1;
	int m_DispThr2;
	int m_Stopway1;
	int m_Stopway2;
	CString m_RunwayLight1;
	CString m_RunwayLight2;
	CString m_ApproachLight1;
	CString m_ApproachLight2;
	CString m_GlideslopeLight1;
	CString m_GlideslopeLight2;
	CString m_SurfaceType;
	CString m_Marking;
	CString m_Should;
	float m_Roughness;
	BOOL m_Signage;
	CNumEdit m_cEdit[10];
	afx_msg void OnBnClickedOK();
	virtual BOOL OnInitDialog();
	CComboBox m_RunwayLight1_ctrl;
	CComboBox m_ApproachLight1_ctrl;
	CComboBox m_GlideslopeLight1_ctrl;
	CComboBox m_RunwayLight2_ctrl;
	CComboBox m_ApproachLight2_ctrl;
	CComboBox m_GlideslopeLight2_ctrl;
	CComboBox m_SurfaceType_ctrl;
	int  m_RunwayAltitude;
	CComboBox m_Should_ctrl;
	afx_msg HBRUSH OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnPaint();
	//CButtonST m_CHECK1;
	CThemeHelperST m_Theme;
	CXPStyleButtonST m_OK;
	CXPStyleButtonST m_cancel;

};
