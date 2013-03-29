//
//	Class:		CStaticST
//
//	Compiler:	Visual C++
//	Tested on:	Visual C++ 6.0
//
//	Version:	See GetVersionC() or GetVersionI()
//
//	Created:	11/November/2002
//	Updated:	29/November/2002
//
//	Author:		Davide Calabro'		davide_calabro@yahoo.com
//
#ifndef _STATICST_H_
#define _STATICST_H_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

// Return values
#ifndef	STATICST_OK
#define	STATICST_OK						0
#endif
#ifndef	STATICST_INVALIDCLICKTYPE
#define	STATICST_INVALIDCLICKTYPE		1
#endif
#ifndef	STATICST_BADPARAM
#define	STATICST_BADPARAM				2
#endif
#ifndef	STATICST_INVALIDRESOURCE
#define	STATICST_INVALIDRESOURCE		3
#endif

class CStaticST : public CStatic
{
public:
	CStaticST();
	virtual ~CStaticST();

	enum	{	STATICST_CLICK	= 0,		// A single click on the control
				STATICST_DBLCLICK,			// A double click on the control

				STATICST_MAX_CLICKS
			};

	DWORD SetIcon(int nIcon, BOOL bRepaint = TRUE, HINSTANCE hInstance = NULL);
	DWORD SetIcon(HICON hIcon, BOOL bRepaint = TRUE);

	void SetTooltipText(int nId, BOOL bActivate = TRUE);
	void SetTooltipText(LPCTSTR lpszText, BOOL bActivate = TRUE);
	void ActivateTooltip(BOOL bActivate = TRUE);

	void SetColors(COLORREF crTextColor, COLORREF crBkColor, BOOL bRepaint = TRUE);
	void SetTextColor(COLORREF crTextColor, BOOL bRepaint = TRUE);
	void SetBkColor(COLORREF crBkColor, BOOL bRepaint = TRUE);

	void SetBold(BOOL bBold, BOOL bRepaint = TRUE);

	DWORD SetMessageClick(BYTE byClickType, HWND hWnd, UINT nMsg);

	DWORD DrawTransparent(BOOL bRepaint = FALSE);
	DWORD SetBk(CDC* pDC);

	static short GetVersionI()		{return 10;}
	static LPCTSTR GetVersionC()	{return (LPCTSTR)_T("1.0");}


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CStaticST)
	public:
	virtual BOOL PreTranslateMessage(MSG* pMsg);
	//}}AFX_VIRTUAL

protected:
	virtual DWORD OnDrawBackground(CDC* pDC, LPCRECT pRect);

	//{{AFX_MSG(CStaticST)
	afx_msg void OnPaint();
	afx_msg void OnNcPaint();
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDblClk(UINT nFlags, CPoint point);
	afx_msg HBRUSH CtlColor(CDC* pDC, UINT nCtlColor);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	//}}AFX_MSG

	CFont			m_Font;

	BOOL			m_bDrawTransparent;	// Draw transparent?
	COLORREF		m_crTextColor;
	COLORREF		m_crBkColor;

	HICON			m_hIcon;
	BYTE			m_cyIcon;
	BYTE			m_cxIcon;

private:
	void FreeResources(BOOL bCheckForNULL = TRUE);
	void PaintBk(CDC* pDC);
	void InitToolTip();
	void DrawItem(CDC* pDC);
	void DrawTheIcon(CDC* pDC, CRect* rpCtrl, BOOL bText);
	void DrawTheText(CDC* pDC, CRect* rpCtrl, CString& sText, DWORD dwStyle, DWORD dwExStyle);

	CToolTipCtrl	m_ToolTip;

	CDC				m_dcBk;
	CBitmap			m_bmpBk;
	CBitmap*		m_pbmpOldBk;

#pragma pack(1)
	typedef struct _STRUCT_CLICKS
	{
		HWND		hWnd;			// Window that will receive the message
		UINT		nMsg;			// Message to send
	} STRUCT_CLICKS;
#pragma pack()

	STRUCT_CLICKS	m_csClicks[STATICST_MAX_CLICKS];

	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif
